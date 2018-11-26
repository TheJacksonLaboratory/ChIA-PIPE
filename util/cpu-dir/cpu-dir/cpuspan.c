// TODO:
// 1. filter by linker type
// 2. summary collectively for a list of files
// 3. handle the case of sorted by name, and sorted by coordinates
// 4. R-script to generate the graphics!!!

#include "zlib.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h> // for multi-threading
#include <errno.h>
#include <math.h> // for floor
#include "utils.h"
#include "kseq.h"
#include "kstring.h"
#include "bwa.h"
KSEQ_DECLARE(gzFile)

#include "sam_header.h"
#include "sam.h"
#include "cputagcommon.h"
#include "cputagbam.h"
#include "cpuspan.h"

extern double G_t_real;
extern const char *CPU_version;

#define CPU_SPAN_OUTPUT_SELF_LIGATION	0
#define CPU_SPAN_OUTPUT_INTRA_LIGATION	1
#define CPU_SPAN_OUTPUT_INTER_LIGATION	2
#define CPU_SPAN_OUTPUT_COUNT			(CPU_SPAN_OUTPUT_INTER_LIGATION+1)
#define CPU_SPAN_OUTPUT_NONE			CPU_SPAN_OUTPUT_COUNT

#define NUM_PAIR_ORIENTATIONS	4
#define PF_STATES_PASS	0
#define PF_STATES_FAIL	1
#define NUM_PASS_FAIL_STATES	(PF_STATES_FAIL+1)

// we will not investiage the combination of <pair_orientation> x <linker combination>
#define CPU_SPAN_LINKER_AA	0
#define CPU_SPAN_LINKER_BB	1
#define CPU_SPAN_LINKER_AB	2
#define CPU_SPAN_LINKER_BA	3
#define CPU_SPAN_LINKER_NONE	4
#define NUM_LINKER_COMBINATIONS	5
#define NUM_OF_FIELDS_PER_SPAN_BIN (NUM_PAIR_ORIENTATIONS+NUM_LINKER_COMBINATIONS)
#define LINKER_OFFSET_IN_BIN(o,l) ((o)*NUM_OF_FIELDS_PER_SPAN_BIN+NUM_PAIR_ORIENTATIONS+l)
#define PAIR_ORIENTATION_OFFSET_IN_BIN(o,p) ((o)*NUM_OF_FIELDS_PER_SPAN_BIN+p)

typedef struct {
	int chunk_size;
	int n_threads;
	
	int lower;
	int upper;
	int n_bins;
	int selfLigation;
	
	int toGroup;
	
	kstring_t outputPrefix;
	
	pairqc_opt_t pairqc_opt;
} span_opt_t;

span_opt_t *span_opt_init()
{
	span_opt_t *o;
	o = calloc(1, sizeof(span_opt_t));
	
	o->chunk_size = 1000000;
	o->n_threads = 1;
	
	o->pairqc_opt.minFractionalMAPQ = 0.95;
	//o->pairqc_opt.minMAPQ = 19;
	o->pairqc_opt.minMAPQ = -1;
	o->pairqc_opt.maxX1 = 10; // equavalent to mapq 19
	
	o->lower = 1;
	o->upper = 100000;
	o->n_bins = 1000;
	o->selfLigation = G_SELF_LIGATION;
	
	o->toGroup = 0;
	
	memset(&(o->outputPrefix), 0, sizeof(kstring_t));
	return o;
}

void span_opt_terminate(span_opt_t *o)
{
	free(o->outputPrefix.s);
}

typedef struct {
	long n_differentChrom[NUM_PASS_FAIL_STATES];
	long n_sameChrom[NUM_PASS_FAIL_STATES];
	long n_sameChromLow[NUM_PASS_FAIL_STATES];
	long n_sameChromHigh[NUM_PASS_FAIL_STATES];
	
	double lower;
	double upper;
	double interval;
	
	int nBins;
	long bin_under[NUM_OF_FIELDS_PER_SPAN_BIN];
	long bin_over[NUM_OF_FIELDS_PER_SPAN_BIN];
	long *bins;
} span_histogram_t;

static inline void init_Span_Histogram (int lower, int upper, int nBins, span_histogram_t *h)
{
	int32_t x;
	h->n_differentChrom[PF_STATES_PASS] = h->n_differentChrom[PF_STATES_FAIL] = 0;
	h->n_sameChrom[PF_STATES_PASS] = h->n_sameChrom[PF_STATES_FAIL] = 0;
	h->n_sameChromLow[PF_STATES_PASS] = h->n_sameChromLow[PF_STATES_FAIL] = 0;
	h->n_sameChromHigh[PF_STATES_PASS] = h->n_sameChromHigh[PF_STATES_FAIL] = 0;
	
	h->lower = lower;
	h->upper = upper;
	h->interval = (h->upper-h->lower+1)/(1.0f*nBins);
	h->nBins = nBins;
	for(x=0; x<NUM_OF_FIELDS_PER_SPAN_BIN; ++x) { h->bin_under[x] = h->bin_over[x] = 0; }
	h->bins = calloc(NUM_OF_FIELDS_PER_SPAN_BIN*nBins, sizeof(long)); // consider pairing information
	if (!h->bins) {
		fprintf(stderr, "[E::%s] Failed to allocate %d bins\n", __func__, 4*nBins);
	}
}

static inline int binOffset (span_histogram_t *h, int span)
{
	int offset;
	
	if (span<0) span*=-1;
	
	offset = (int) floor(((span-h->lower)/h->interval)+0.5);
	if (offset<0) offset = -1;
	if (offset>=h->nBins) offset = h->nBins;
	
	return offset;
}

static inline int binPairOffset (span_histogram_t *h, int span, int pairOrientation)
{
	int offset = binOffset(h, span);
	
	if (-1!=offset && h->nBins!=offset) {
		offset *= NUM_OF_FIELDS_PER_SPAN_BIN;
		offset += pairOrientation;
	}
	
	return offset;
}

static inline void binInterval (span_histogram_t *h, int bin, int *lower, int *upper)
{
	if (0>bin) {
		int value = (0 - 0.5) * h->interval + h->lower - 1;
		if (value<0) {
			// already underflow with bin zero, and thus we will have no bin_under
			*lower = *upper = -1;
		} else {
			//value>=0
			*upper = value;
			*lower = 0;
		}
	} else if (bin>=h->nBins) {
		int value = (h->nBins - 0.5) * h->interval + h->lower;
		*lower = value;
		*upper = -1;
	} else {
		int value = (bin - 0.5) * h->interval + h->lower;
		*lower = value;
		value = (bin + 1 - 0.5) * h->interval + h->lower - 1;
		*upper = value;
	}
}

int32_t getPairOrientation (bam1_core_t *bamCore) {
	int pairOrientation = 0;
	if (bamCore->isize<0) {
		// negative span, i.e. pos > mpos
		// but we prefer to see the condition as pos < mpos
		if (BAM_FREVERSE==(BAM_FREVERSE&bamCore->flag)) {
			pairOrientation |= 0x01;
		}
		if (BAM_FMREVERSE==(BAM_FMREVERSE&bamCore->flag)) {
			pairOrientation |= 0x02;
		}
	} else {
		// positive span, i.e. pos < mpos
		if (BAM_FREVERSE==(BAM_FREVERSE&bamCore->flag)) {
			pairOrientation |= 0x02;
		}
		if (BAM_FMREVERSE==(BAM_FMREVERSE&bamCore->flag)) {
			pairOrientation |= 0x01;
		}
	}
	return pairOrientation;
}

static inline void count_Histogram (span_histogram_t *h, int span, int pairOrientation, int finalLinker)
{
	// TODO: consider pairing
	int offset = binOffset (h, span);
	if (-1==offset) {
		h->bin_under[PAIR_ORIENTATION_OFFSET_IN_BIN(0,pairOrientation)]++;
		h->bin_under[LINKER_OFFSET_IN_BIN(0,finalLinker)]++;
	}
	else if (h->nBins==offset) {
		h->bin_over[PAIR_ORIENTATION_OFFSET_IN_BIN(0,pairOrientation)]++;
		h->bin_over[LINKER_OFFSET_IN_BIN(0,finalLinker)]++;
	}
	else {
		h->bins[PAIR_ORIENTATION_OFFSET_IN_BIN(offset,pairOrientation)]++;
		h->bins[LINKER_OFFSET_IN_BIN(offset,finalLinker)]++;
	}
}

static inline void terminate_Span_Histogram (span_histogram_t *h)
{
	free(h->bins);
}

typedef struct {
	const span_opt_t *opt;
	
	int64_t n_processed;
	
	bam1_t *bamRecs;
	cpu_readid_t *readids;
	
	int nTagGroups;
	tag_group_t *tagGroups;
	
	span_histogram_t *spanHistograms;
	
	int (*pairQCCheck)(const bam1_t *,const bam1_t *,const pairqc_opt_t *);
} worker_t;

static inline int init_CPSpan_Outputs (const char *outPrefix, char *out_mode, bam_header_t *header, samfile_t **outfds, int selfLigation) {
	// CPU_SPAN_OUTPUT_SELF_LIGATION
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.self.%d.bam", outPrefix, selfLigation);
		outfds[CPU_SPAN_OUTPUT_SELF_LIGATION] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_SPAN_OUTPUT_SELF_LIGATION]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_SPAN_OUTPUT_SELF_LIGATION+1;
		}
		free(filename.s);
	}
	// CPU_SPAN_OUTPUT_INTRA_LIGATION
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.intra.%d.bam", outPrefix, selfLigation);
		outfds[CPU_SPAN_OUTPUT_INTRA_LIGATION] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_SPAN_OUTPUT_INTRA_LIGATION]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_SPAN_OUTPUT_INTRA_LIGATION+1;
		}
		free(filename.s);
	}
	// CPU_SPAN_OUTPUT_INTER_LIGATION
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.inter.bam", outPrefix);
		outfds[CPU_SPAN_OUTPUT_INTER_LIGATION] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_SPAN_OUTPUT_INTER_LIGATION]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_SPAN_OUTPUT_INTER_LIGATION+1;
		}
		free(filename.s);
	}
	
	return 0;
}

static inline int terminate_CPSpan_Outputs (samfile_t **outfds) {
	// TODO: always synchronized with #defined
	int nFailed = 0;
	if (outfds[CPU_SPAN_OUTPUT_SELF_LIGATION]) {
		samclose(outfds[CPU_SPAN_OUTPUT_SELF_LIGATION]);
	}
	if (outfds[CPU_SPAN_OUTPUT_INTRA_LIGATION]) {
		samclose(outfds[CPU_SPAN_OUTPUT_INTRA_LIGATION]);
	}
	if (outfds[CPU_SPAN_OUTPUT_INTER_LIGATION]) {
		samclose(outfds[CPU_SPAN_OUTPUT_INTER_LIGATION]);
	}
	
	return nFailed;
}

static void readid_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	// read id parsing does not care about if reads are paried-end
	int nRet = parseReadId(bam1_qname(&(w->bamRecs[i])), &(w->readids[i]), (BAM_FREAD2==(BAM_FREAD2 & w->bamRecs[i].core.flag))?1:0);
	if (0!=nRet) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to parse read id \"%s\". Make sure that read is prepared by CPU.\n", __func__, bam1_qname(&(w->bamRecs[i])));
	}
}

static void span_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	// we are going to assume that these are paired-end data
	// in case, the file is already coordinate-sorted, we will just use read/1 for the proper insert size
	// we will also check that the tag is set by CPU
	
	int j, k;
	for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
		// first in pair, paired, and not secondary!
		bam1_t *bam = &(w->bamRecs[k]);
		if ((BAM_FPAIRED|BAM_FREAD1)==(bam->core.flag & (BAM_FPAIRED|BAM_FREAD1|BAM_FUNMAP|BAM_FMUNMAP|BAM_FSECONDARY))) {
			uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
			if (0!=yt) {
				char *a_yt = bam_aux2Z(yt);
				if (0==strcmp(a_yt, G_PAIRING_TAG_UU)) {
					bam1_t *mbam = w->tagGroups[i].bamRecs[1];
					int mapqState = w->pairQCCheck(bam, mbam, &(w->opt->pairqc_opt));
					
					if (bam->core.tid!=bam->core.mtid) {
						w->spanHistograms[tid].n_differentChrom[mapqState]++;
						w->tagGroups[i].outputClass = (PF_STATES_PASS==mapqState) ? CPU_SPAN_OUTPUT_INTER_LIGATION : CPU_SPAN_OUTPUT_NONE;
					} else {
						int32_t nPairOrientation = getPairOrientation(&bam->core);
						
						// let's count the span, exclude different chromosome case
						count_Histogram (&(w->spanHistograms[tid]), bam->core.isize, nPairOrientation, w->tagGroups[i].finalLinker);
						
						w->spanHistograms[tid].n_sameChrom[mapqState]++;
						int32_t isize = (bam->core.isize>=0) ? bam->core.isize : -1*bam->core.isize;
						if (isize<w->opt->selfLigation) {
							w->spanHistograms[tid].n_sameChromLow[mapqState]++;
							w->tagGroups[i].outputClass = (PF_STATES_PASS==mapqState) ? CPU_SPAN_OUTPUT_SELF_LIGATION : CPU_SPAN_OUTPUT_NONE;
						} else {
							w->spanHistograms[tid].n_sameChromHigh[mapqState]++;
							w->tagGroups[i].outputClass = (PF_STATES_PASS==mapqState) ? CPU_SPAN_OUTPUT_INTRA_LIGATION : CPU_SPAN_OUTPUT_NONE;
						}
					}
					break;
				}
			} else {
				// data has not been processed by CPU!!!
				// TODO: we might have to a lenient mode for passing through
				//       but in this case, we likely need to compute the span and determine if we could even use the data
				//if (bwa_verbose >= 1)
				fprintf(stderr, "[E::%s] %s tag unavailable for \"%s\". Make sure that mapper produce %s tag like CPU..\n", __func__, G_PAIRING_TAG, bam1_qname(bam), G_PAIRING_TAG);
			}
		}
	}
	
}

static inline int finalLinkerToSpanLinker (int16_t finalLinker)
{
	if (LINKER_AA==(finalLinker&LINKER_AA)) {
		return CPU_SPAN_LINKER_AA;
	} else if (LINKER_BB==(finalLinker&LINKER_BB)) {
		return CPU_SPAN_LINKER_BB;
	} else if (LINKER_AB==(finalLinker&LINKER_AB)) {
		return CPU_SPAN_LINKER_AB;
	} else if (LINKER_BA==(finalLinker&LINKER_BA)) {
		return CPU_SPAN_LINKER_BA;
	} else {
		return CPU_SPAN_LINKER_NONE;
		/*fprintf(stderr, "[E::%s] Unrecognized final linker type (%d,%x) in read id\n", __func__, finalLinker, finalLinker);
		return -1;*/
	}
}

void span_process_alns(const span_opt_t *opt, int64_t n_processed, int n, bam1_t *bamRecs, cpu_readid_t *readids, int *pNumTagGroups, tag_group_t *tagGroups, CPUBamBuffer_t* bamsBuffer, span_histogram_t *tspans, span_histogram_t *spans, int (*pairQCCheck)(const bam1_t *,const bam1_t *,const pairqc_opt_t *))
{
	int i, j;
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	if (n<=0) return;
	
	w.opt = opt;
	w.n_processed = n_processed; // TODO:
	// TODO: w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	w.bamRecs = bamRecs;
	w.readids = readids;
	w.spanHistograms = tspans;
	w.pairQCCheck = pairQCCheck;
	
	// parse the read id
	kt_for(opt->n_threads, readid_worker, &w, n);
	
	// tag grouping
	// TODO: parallelize later
	int numTagGroups = 0;
	long pairId = readids[0].pairId;
	tagGroups[numTagGroups].readIndex = 0; tagGroups[numTagGroups].n_reads = 1; tagGroups[numTagGroups].bamRecs[0] = &(bamRecs[0]);
	tagGroups[numTagGroups].finalLinker = finalLinkerToSpanLinker(readids[0].finalLinkerType);
	int currTagGroup = numTagGroups; numTagGroups++;
	for(i=1; i<n; ++i) {
		if (pairId == readids[i].pairId) {
			tagGroups[currTagGroup].n_reads++;
			if (tagGroups[currTagGroup].n_reads>2) {
				fprintf(stderr, "[W::%s] Pair id %ld contains %d records..\n", __func__, pairId, tagGroups[currTagGroup].n_reads);
			} else if (2==tagGroups[currTagGroup].n_reads) {
				tagGroups[currTagGroup].bamRecs[1] = &(bamRecs[i]);
			}
		} else {
			pairId = readids[i].pairId;
			//CPU_SPAN_LINKER_AA
			tagGroups[numTagGroups].finalLinker =finalLinkerToSpanLinker(readids[i].finalLinkerType);
			tagGroups[numTagGroups].readIndex = i; tagGroups[numTagGroups].n_reads = 1;
			tagGroups[numTagGroups].bamRecs[0] = &(bamRecs[i]);
			currTagGroup++; numTagGroups++;
		}
	}
	
	if (!bamsBuffer->lastBlock) {
		// this is not the last block, and thus we cannot guarantee the completeness of the last tagGroup!
		// let's cache the last block for the next iteration
		numTagGroups--;
		add_CPUBamBuffer(&(bamRecs[tagGroups[numTagGroups].readIndex]), tagGroups[numTagGroups].n_reads, bamsBuffer);
	}
	*pNumTagGroups = numTagGroups;
	
	// pairing selection
	w.nTagGroups = numTagGroups;
	w.tagGroups = tagGroups;
	kt_for(opt->n_threads, span_worker, &w, numTagGroups);
	
	// accumulate all the span information
	for(i=0; i<opt->n_threads; ++i) {
		for(j=0; j<NUM_OF_FIELDS_PER_SPAN_BIN*spans->nBins;++j) { // there are 4 possible pairings
			spans->bins[j] += tspans[i].bins[j];
		}
		for(j=0; j<NUM_OF_FIELDS_PER_SPAN_BIN; ++j) {
			spans->bin_under[j] += tspans[i].bin_under[j];
			spans->bin_over[j] += tspans[i].bin_over[j];
		}
		
		for(j=0; j<NUM_PASS_FAIL_STATES; ++j) {
			spans->n_differentChrom[j] += tspans[i].n_differentChrom[j];
			spans->n_sameChrom[j] += tspans[i].n_sameChrom[j];
			spans->n_sameChromLow[j] += tspans[i].n_sameChromLow[j];
			spans->n_sameChromHigh[j] += tspans[i].n_sameChromHigh[j];
		}
	}
}

int main_span(int argc, char *argv[])
{
	int c, n, i, j, k;
	span_opt_t *opt;
	
	bam1_t *bamRecs;
	CPUBamBuffer_t bamsBuffer;
	
	double t_diff;
	double t_timepointIO, t_timepointProcessing;
	double t_diffIO, t_diffProcessing;
	
	int64_t n_processed = 0;
	int64_t n_pairProcessed = 0;
	
	int is_bamin = 1;
	char in_mode[5];
	char *fn_list = 0;
	samfile_t *in = 0;
	
	int ret = 0;
	int is_SpanBin = 0;
	
	int compress_level = -1;
	int is_header = 0;
	int is_bamout = 1;
	char out_mode[5];
	samfile_t *outfds[CPU_SPAN_OUTPUT_COUNT] = {0};
	
	int (*pairQCCheck)(const bam1_t *,const bam1_t *,const pairqc_opt_t *) = getPairQCStateByUniqueLoci;
	
	init_CPUBamBuffer(&bamsBuffer);
	
	opt = span_opt_init();
	strcpy(in_mode, "r");
	strcpy(out_mode, "w");
	while ((c = getopt(argc, argv, "giSs:t:O:l:u:b:q:a:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'g') opt->toGroup = 1;
		else if (c == 'S') is_bamin = 0;
		else if (c == 'l') opt->lower = atoi(optarg);
		else if (c == 'u') opt->upper = atoi(optarg);
		else if (c == 'b') opt->n_bins = atoi(optarg);
		else if (c == 'O') { ks_set(optarg, &(opt->outputPrefix)); }
		else if (c == 'i') is_SpanBin = 1;
		else if (c == 's') opt->selfLigation = atoi(optarg), opt->selfLigation = opt->selfLigation > 0 ? opt->selfLigation : G_SELF_LIGATION;
		else if (c == 'q') {
			if (strstr(optarg, ".")) {
				opt->pairqc_opt.minFractionalMAPQ = atof(optarg);
				opt->pairqc_opt.minMAPQ = -1;
				pairQCCheck = getPairQCStateByUniqueLoci;
			}
			else {
				opt->pairqc_opt.minMAPQ = atoi(optarg);
				opt->pairqc_opt.minFractionalMAPQ = -1.0;
				pairQCCheck = getPairQCStateByMapQ;
			}
		}
		else if (c == 'a') opt->pairqc_opt.maxX1 = atoi(optarg);
		else {
			span_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	
	if (is_SpanBin) {
		span_histogram_t spans; init_Span_Histogram (opt->lower, opt->upper, opt->n_bins, &spans);
		fprintf(stdout, "# span-bin distribution\n");
		fprintf(stdout, "#span\tbin\n");
		for(i=opt->lower; i<=opt->upper; ++i) {
			int nBin = binOffset(&spans, i);
			fprintf(stdout, "%d\t%d\n", i, nBin);
		}
		terminate_Span_Histogram (&spans);
		
		span_opt_terminate(opt);
		free(opt);
		return 0;
	}
	
	if (opt->n_threads < 1) opt->n_threads = 1;
	
	if (is_bamin) strcat(in_mode, "b");
	if (is_bamout) strcat(out_mode, "b");
	if (is_header) strcat(out_mode, "h");
	if (compress_level >= 0) {
		char tmp[2];
		tmp[0] = compress_level + '0'; tmp[1] = '\0';
		strcat(out_mode, tmp);
	}
	
	
	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu span [options] <in.sam/.bam>\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -l INT     lower bound [%d]\n", opt->lower);
		fprintf(stderr, "       -u INT     upper bound [%d]\n", opt->upper);
		fprintf(stderr, "       -b INT     number of bins [%d]\n", opt->n_bins);
		fprintf(stderr, "       -s INT     self-ligation distance [%d bp]\n", opt->selfLigation);
		fprintf(stderr, "       -S         input sam format\n");
		fprintf(stderr, "       -g         output different categories read to separate file\n");
		fprintf(stderr, "       -q INT     minimum MAPQ [>=#] or second(score)/best(score)<%.2f\n", opt->pairqc_opt.minFractionalMAPQ);
		fprintf(stderr, "       -a INT     max. number of secondary alternative allowed [<=%d]\n", opt->pairqc_opt.maxX1);
		fprintf(stderr, "\n");
		fprintf(stderr, "       -i         show span bin (for debugging)\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		span_opt_terminate(opt);
		free(opt);
		return 1;
	}
	
	// open file handlers
	if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open \"%s\" for reading.\n", __func__, argv[optind]);
		
		free(fn_list);
		span_opt_terminate(opt);
		free(opt);
		
		return 1;
	}
	if (in->header == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to read the header from \"%s\".\n", __func__, argv[optind]);
		
		free(fn_list);
		// close files, free and return
		samclose(in);
		
		span_opt_terminate(opt);
		free(opt);
		
		return 1;
	}
	
	if (opt->toGroup) {
		// construct the out filename from the in filename
		{
			if (NULL==opt->outputPrefix.s || 0==opt->outputPrefix.l) {
				ks_set(argv[optind], &(opt->outputPrefix));
				char *pch = strrchr(opt->outputPrefix.s, '.');
				if (NULL==pch) {
					// there is no extension, use as-is
				} else {
					if (0==strcmp(pch, ".bam") || 0==strcmp(pch, ".sam")) {
						*pch = '\0';
						opt->outputPrefix.l -= 4;
					} else if (0==strcmp(pch, ".gz")) {
						*pch = '\0';
						opt->outputPrefix.l -= 3;
						char *pch1 = strrchr(opt->outputPrefix.s, '.');
						if (NULL!=pch1) {
							// remove one more extension
							*pch1 = '\0';
							opt->outputPrefix.l = pch1 - opt->outputPrefix.s;
						}
					} else {
						*pch = '\0';
						opt->outputPrefix.l = pch - opt->outputPrefix.s;
					}
				}
			}
		}
		if (init_CPSpan_Outputs (opt->outputPrefix.s, out_mode, in->header, outfds, opt->selfLigation)) {
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open for writing.\n", __func__);
			
			free(fn_list);
			// close files, free and return
			samclose(in);
			
			span_opt_terminate(opt);
			free(opt);
			
			return 1;
		}
		
		{
			if (opt->n_threads > 1) samthreads(outfds[CPU_SPAN_OUTPUT_SELF_LIGATION], opt->n_threads, 256);
			if (opt->n_threads > 1) samthreads(outfds[CPU_SPAN_OUTPUT_INTRA_LIGATION], opt->n_threads, 256);
			if (opt->n_threads > 1) samthreads(outfds[CPU_SPAN_OUTPUT_INTER_LIGATION], opt->n_threads, 256);
		}
	}
	
	// PROCESSING
	span_histogram_t spans; init_Span_Histogram (opt->lower, opt->upper, opt->n_bins, &spans);
	if (argc == optind + 1) { // convert/print the entire file
		//int nNumRequested = opt->chunk_size * opt->n_threads;
		int nNumRequested = opt->chunk_size;
		t_timepointProcessing = realtime();
		while ((bamRecs = readCPUBam(&n, nNumRequested, in, &bamsBuffer)) != 0) {
			t_timepointIO = realtime();
			t_diffIO = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d records..\n", __func__, n-bamsBuffer.unprocessed);
			
			// processing the pairing
			cpu_readid_t *readids = calloc(n, sizeof(cpu_readid_t));
			tag_group_t *tagGroups = calloc(n, sizeof(tag_group_t));
			
			span_histogram_t *tspans = calloc(opt->n_threads, sizeof(span_histogram_t));
			for(i=0; i<opt->n_threads; ++i) { init_Span_Histogram (opt->lower, opt->upper, opt->n_bins, &(tspans[i])); }
			
			int numTagGroups = 0;
			span_process_alns(opt, n_processed, n, bamRecs, readids, &numTagGroups, tagGroups, &bamsBuffer, tspans, &spans, pairQCCheck);
			t_timepointProcessing = realtime();
			t_diffProcessing = t_timepointProcessing - t_timepointIO;
			
			if (opt->toGroup) {
				// TODO: write the bam results
				for(i=0; i<numTagGroups; ++i) {
					if (CPU_SPAN_OUTPUT_NONE!=tagGroups[i].outputClass) {
						samfile_t *outfd = outfds[tagGroups[i].outputClass];
						for(j=0, k=tagGroups[i].readIndex; j<tagGroups[i].n_reads; ++j, ++k) {
							samwrite(outfd, &bamRecs[k]);
						}
					}
				}
				t_timepointIO = realtime();
				t_diffIO += (t_timepointIO - t_timepointProcessing);
			}
			
			// clean up
			for(i=0; i<n; ++i) free(bamRecs[i].data);
			free(bamRecs);
			free(readids);
			for(i=0; i<opt->n_threads; ++i) { terminate_Span_Histogram (&(tspans[i])); }
			free(tspans);
			free(tagGroups);
			
			n_processed += (n-bamsBuffer.unprocessed);
			n_pairProcessed += numTagGroups;
			if (bwa_verbose >= 3) {
				t_diff = realtime() - G_t_real;
				fprintf(stderr, "[M::%s] %lld tags %lld pairs processed, %.0f tags/sec %.0f pairs/sec, %.2f min, i/o %.2f sec, processing %.2f sec..\n", __func__, n_processed, n_pairProcessed, 1.0*n_processed/t_diff, 1.0*n_pairProcessed/t_diff, t_diff/60.0, t_diffIO, t_diffProcessing);
			}
			
			t_timepointProcessing = realtime();
		}
	}
	// END - PROCESSING
	
	free(fn_list);
	
	// close files, free and return
	samclose(in);
	
	terminate_CPUBamBuffer(&bamsBuffer);
	
	// TODO: report the span statistics
	int64_t withLinker = 0;
	int64_t withoutLinker = 0;
	for(i=0;i<spans.nBins;i+=NUM_OF_FIELDS_PER_SPAN_BIN) {
		if (spans.bins[i+NUM_PAIR_ORIENTATIONS]>0) withLinker+=spans.bins[i+NUM_PAIR_ORIENTATIONS];
		if (spans.bins[i+NUM_PAIR_ORIENTATIONS+1]>0) withLinker+=spans.bins[i+NUM_PAIR_ORIENTATIONS+1];
		if (spans.bins[i+NUM_PAIR_ORIENTATIONS+2]>0) withLinker+=spans.bins[i+NUM_PAIR_ORIENTATIONS+2];
		if (spans.bins[i+NUM_PAIR_ORIENTATIONS+3]>0) withLinker+=spans.bins[i+NUM_PAIR_ORIENTATIONS+3];
		if (spans.bins[i+NUM_PAIR_ORIENTATIONS+4]>0) withoutLinker+=spans.bins[i+NUM_PAIR_ORIENTATIONS+4];
	}
	int linkerTypeSeen = 0;
	if (withLinker) linkerTypeSeen |= 1;
	if (withoutLinker) linkerTypeSeen |= 2;
	if (3==linkerTypeSeen) {
		if (bwa_verbose >= 2)
			fprintf(stderr, "[W::%s] Mixture of read-pairs with and without linkers detected (%lld vs %lld)\n", __func__, withLinker, withoutLinker);
	}
	int64_t totalInBins = 0;
	for(i=0;i<spans.nBins;i+=NUM_OF_FIELDS_PER_SPAN_BIN) {
		// 4 pair-orientation
		totalInBins += spans.bins[i];
		totalInBins += spans.bins[i+1];
		totalInBins += spans.bins[i+2];
		totalInBins += spans.bins[i+3];
	}
	if (0==totalInBins) {
		if (bwa_verbose >= 2)
		fprintf(stderr, "[W::%s] %lld in bins of %lld pairs (%lld tags). Has \"%s\" been paired?\n", __func__, totalInBins, n_pairProcessed, n_processed, argv[optind]);
	}
	fprintf(stdout, "##CPU\t%s\n", CPU_version);
	{
		fprintf(stdout, "##COMMAND\t");
		for (i = 0; i < argc; ++i)
		fprintf(stdout, " %s", argv[i]);
		fprintf(stdout, "\n");
	}
	
	long nTotal = n_pairProcessed;
	if (0==nTotal) nTotal = 1; // prevent division by zero!
	fprintf(stdout, ">>Library information\n");
	fprintf(stdout, "#Measure\tValue\n");
	fprintf(stdout, "Filename\t%s\n", argv[optind]);
	fprintf(stdout, ">>END\n");
	
	fprintf(stdout, ">>General Statistics\n");
	fprintf(stdout, "#Category\tCount\tPercent\n");
	fprintf(stdout, "Total tags\t%lld\n", n_processed);
	fprintf(stdout, "Total pairs\t%lld\n", n_pairProcessed);
	
	// general
	if (-1==opt->pairqc_opt.minFractionalMAPQ) {
		fprintf(stdout, "mapq>=%d\tCount\tPercent\tPassed\t%%Passed\tFailed\t%%Failed\n", opt->pairqc_opt.minMAPQ);
	} else {
		fprintf(stdout, "second/best<%.2f\tCount\tPercent\tPassed\t%%Passed\tFailed\t%%Failed\n", opt->pairqc_opt.minFractionalMAPQ);
	}
	long lValue = spans.n_differentChrom[PF_STATES_PASS]+spans.n_differentChrom[PF_STATES_FAIL];
	fprintf(stdout, "Different chromosomes\t%ld\t%.2f%%\t%ld\t%.2f%%\t%ld\t%.2f%%\n", lValue, lValue*100.0f/nTotal,
			spans.n_differentChrom[PF_STATES_PASS],spans.n_differentChrom[PF_STATES_PASS]*100.0f/nTotal,
			spans.n_differentChrom[PF_STATES_FAIL],spans.n_differentChrom[PF_STATES_FAIL]*100.0f/nTotal);
	lValue = spans.n_sameChrom[PF_STATES_PASS]+spans.n_sameChrom[PF_STATES_FAIL];
	fprintf(stdout, "Same chromosome\t%ld\t%.2f%%\t%ld\t%.2f%%\t%ld\t%.2f%%\n", lValue, lValue*100.0f/nTotal,
			spans.n_sameChrom[PF_STATES_PASS],spans.n_sameChrom[PF_STATES_PASS]*100.0f/nTotal,
			spans.n_sameChrom[PF_STATES_FAIL],spans.n_sameChrom[PF_STATES_FAIL]*100.0f/nTotal);
	{
		// adjust for same chromosome breakdown
		long lTotal = spans.n_sameChrom[PF_STATES_PASS];
		if (0==lTotal) lTotal = 1; // prevent division by zero!
		kstring_t condition = {0,0,0};
		if (0==(opt->selfLigation%1000)) {
			ksprintf(&condition, "%d Kb", opt->selfLigation/1000);
		} else if (0==(opt->selfLigation%100)) {
			ksprintf(&condition, "%.1f Kb", opt->selfLigation/100.0);
		} else {
			ksprintf(&condition, "%d b", opt->selfLigation);
		}
		lValue = spans.n_sameChromLow[PF_STATES_PASS]+spans.n_sameChromLow[PF_STATES_FAIL];
		fprintf(stdout, "Same chromosome (<%s)\t%ld\t%.2f%%\t%ld\t%.2f%%\t%ld\t%.2f%%\n", condition.s, lValue, lValue*100.0f/lTotal,
				spans.n_sameChromLow[PF_STATES_PASS],spans.n_sameChromLow[PF_STATES_PASS]*100.0f/nTotal,
				spans.n_sameChromLow[PF_STATES_FAIL],spans.n_sameChromLow[PF_STATES_FAIL]*100.0f/nTotal);
		lValue = spans.n_sameChromHigh[PF_STATES_PASS]+spans.n_sameChromHigh[PF_STATES_FAIL];
		fprintf(stdout, "Same chromosome (>=%s)\t%ld\t%.2f%%\t%ld\t%.2f%%\t%ld\t%.2f%%\n", condition.s, lValue, lValue*100.0f/lTotal,
				spans.n_sameChromHigh[PF_STATES_PASS],spans.n_sameChromHigh[PF_STATES_PASS]*100.0f/nTotal,
				spans.n_sameChromHigh[PF_STATES_FAIL],spans.n_sameChromHigh[PF_STATES_FAIL]*100.0f/nTotal);
		free(condition.s);
	}
	// general histogram statistics
	fprintf(stdout, "Histogram range\t[%ld,%ld]\n", (long)spans.lower, (long)spans.upper);
	fprintf(stdout, "Number of histogram bin\t%d\n", spans.nBins);
	fprintf(stdout, "Histogram interval\t%ld\n", (long)spans.interval);
	fprintf(stdout, "Total in bins\t%lld\t%.2g%%\n", totalInBins, totalInBins*100.0/nTotal);
	int lower, upper;
	binInterval (&spans, -1, &lower, &upper);
	// ++0,+-1,-+2,--3 --> +-1,-+2,++0,--3
	//fprintf(stdout, "\tcount\tpercent\tcount[++]\tpercent[++]\tcount[+-]\tpercent[+-]\tcount[-+]\tpercent[-+]\tcount[--]\tpercent[--]\n");
	fprintf(stdout, "\tcount\tpercent\tcount[+-]\tpercent[+-]\tcount[-+]\tpercent[-+]\tcount[++]\tpercent[++]\tcount[--]\tpercent[--]\n");
	if (-1==upper) {
		fprintf(stdout, "Total in underflow bin\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\n");
	} else {
		lValue = 0;
		for(i=0; i<NUM_PAIR_ORIENTATIONS; i++) { lValue += spans.bin_under[i]; }
		fprintf(stdout, "Total in underflow bin [0,%d]\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\n", upper,
				lValue, lValue*100.0/nTotal,
				spans.bin_under[1], spans.bin_under[1]*100.0/nTotal,
				spans.bin_under[2], spans.bin_under[2]*100.0/nTotal,
				spans.bin_under[0], spans.bin_under[0]*100.0/nTotal,
				spans.bin_under[3], spans.bin_under[3]*100.0/nTotal);
	}
	binInterval (&spans, spans.nBins, &lower, &upper);
	lValue = 0;
	for(i=0; i<NUM_PAIR_ORIENTATIONS; i++) { lValue += spans.bin_over[i]; }
	fprintf(stdout, "Total in overflow bin [%d,inf)\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\n", lower,
			lValue, lValue*100.0/nTotal,
			spans.bin_over[1], spans.bin_over[1]*100.0/nTotal,
			spans.bin_over[2], spans.bin_over[2]*100.0/nTotal,
			spans.bin_over[0], spans.bin_over[0]*100.0/nTotal,
			spans.bin_over[3], spans.bin_over[3]*100.0/nTotal);
	
	// report the linker class
	fprintf(stdout, "\tcount\tpercent");
	if (1==(1&linkerTypeSeen)) fprintf(stdout, "\tcount[NC]\tpercent[NC]\tcount[CM]\tpercent[CM]\tcount[AA]\tpercent[AA]\tcount[BB]\tpercent[BB]\tcount[AB]\tpercent[AB]\tcount[BA]\tpercent[BA]");
	if (2==(2&linkerTypeSeen)) fprintf(stdout, "\tcount[no]\tpercent[no]");
	fprintf(stdout, "\n");
	if (-1==upper) {
		fprintf(stdout, "Total in underflow bin\tn.a.\tn.a.");
		if (1==(1&linkerTypeSeen)) fprintf(stdout, "\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.\tn.a.");
		if (2==(2&linkerTypeSeen)) fprintf(stdout, "\tn.a.\tn.a.");
		fprintf(stdout, "\n");
	} else {
		lValue = 0;
		for(i=0; i<NUM_LINKER_COMBINATIONS; i++) { lValue += spans.bin_under[i]; }
		fprintf(stdout, "Total in underflow bin [0,%d]\t%ld\t%.2g%%", upper, lValue, lValue*100.0/nTotal);
		if (1==(1&linkerTypeSeen)) fprintf(stdout, "\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%",
				spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AA)]+spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BB)], (spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AA)]+spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BB)])*100.0/nTotal,
				spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AB)]+spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BA)], (spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AB)]+spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BA)])*100.0/nTotal,
				spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AA)], spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AA)]*100.0/nTotal,
				spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BB)], spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BB)]*100.0/nTotal,
				spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AB)], spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AB)]*100.0/nTotal,
				spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BA)], spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BA)]*100.0/nTotal);
		if (2==(2&linkerTypeSeen)) fprintf(stdout, "\t%ld\t%.2g%%",
				spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_NONE)], spans.bin_under[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_NONE)]*100.0/nTotal);
		fprintf(stdout, "\n");
	}
	binInterval (&spans, spans.nBins, &lower, &upper);
	lValue = 0;
	for(i=0; i<NUM_LINKER_COMBINATIONS; i++) { lValue += spans.bin_over[i]; }
	fprintf(stdout, "Total in overflow bin [%d,inf)\t%ld\t%.2g%%", lower, lValue, lValue*100.0/nTotal);
	if (1==(1&linkerTypeSeen)) fprintf(stdout, "\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%",
			spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AA)]+spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BB)], (spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AA)]+spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BB)])*100.0/nTotal,
			spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AB)]+spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BA)], (spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AB)]+spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BA)])*100.0/nTotal,
			spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AA)], spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AA)]*100.0/nTotal,
			spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BB)], spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BB)]*100.0/nTotal,
			spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AB)], spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_AB)]*100.0/nTotal,
			spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BA)], spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_BA)]*100.0/nTotal);
	if (2==(2&linkerTypeSeen)) fprintf(stdout, "\t%ld\t%.2g%%",
			spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_NONE)], spans.bin_over[LINKER_OFFSET_IN_BIN(0,CPU_SPAN_LINKER_NONE)]*100.0/nTotal);
	fprintf(stdout, "\n");
	// end of linker class
	fprintf(stdout, ">>END\n");
	
	if (totalInBins>0) {
		fprintf(stdout, ">>Span-frequency pair orientation distribution\n");
		// ++0,+-1,-+2,--3 --> +-1,-+2,++0,--3
		//fprintf(stdout, "#bin\tlower\tupper\trepSpan\tcount\tpercent\tcount[++]\tpercent[++]\tcount[+-]\tpercent[+-]\tcount[-+]\tpercent[-+]\tcount[--]\tpercent[--]\n");
		fprintf(stdout, "#bin\tlower\tupper\trepSpan\tcount\tpercent\tcount[+-]\tpercent[+-]\tcount[-+]\tpercent[-+]\tcount[++]\tpercent[++]\tcount[--]\tpercent[--]\tRatio\n");
		// we do not use #sameChrom so as to get the general sense of usable %data
		for(i=0; i<spans.nBins;++i) {
			binInterval (&spans, i, &lower, &upper);
			if (lower<0) lower = 0;
			lValue = 0; k=i*NUM_OF_FIELDS_PER_SPAN_BIN; for(j=k; j<k+NUM_PAIR_ORIENTATIONS; j++) { lValue += spans.bins[j]; }
			// TODO: consider pairing
			float ratio = ((spans.bins[k+2]+spans.bins[k]+spans.bins[k+3])*1.0f) / (spans.bins[k+1] * 1.0f);
			fprintf(stdout, "%d\t%d\t%d\t%d\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%.1f\n", i, lower, upper, (int)((lower+upper)/2),
					lValue, lValue*100.0/nTotal,
					spans.bins[k+1], spans.bins[k+1]*100.0/totalInBins,
					spans.bins[k+2], spans.bins[k+2]*100.0/totalInBins,
					spans.bins[k], spans.bins[k]*100.0/totalInBins,
					spans.bins[k+3], spans.bins[k+3]*100.0/totalInBins,
					ratio);
		}
		fprintf(stdout, ">>END\n");
		
		
		
		fprintf(stdout, ">>Span-frequency linker class distribution\n");
		fprintf(stdout, "#bin\tlower\tupper\trepSpan\tcount\tpercent");
		if (1==(1&linkerTypeSeen)) fprintf(stdout, "\tcount[NC]\tpercent[NC]\tcount[CM]\tpercent[CM]\tcount[AA]\tpercent[AA]\tcount[BB]\tpercent[BB]\tcount[AB]\tpercent[AB]\tcount[BA]\tpercent[BA]\tRatio");
		if (2==(2&linkerTypeSeen)) fprintf(stdout, "\tcount[no]\tpercent[no]");
		fprintf(stdout, "\n");
		// we do not use #sameChrom so as to get the general sense of usable %data
		for(i=0; i<spans.nBins;++i) {
			binInterval (&spans, i, &lower, &upper);
			if (lower<0) lower = 0;
			lValue = 0; k=i*NUM_OF_FIELDS_PER_SPAN_BIN; for(j=k; j<k+NUM_LINKER_COMBINATIONS; j++) { lValue += spans.bins[j]; }
			// TODO: consider pairing
			fprintf(stdout, "%d\t%d\t%d\t%d\t%ld\t%.2g%%", i, lower, upper, (int)((lower+upper)/2), lValue, lValue*100.0/nTotal);
			if (1==(1&linkerTypeSeen)) {
				float ratio = ((spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AA)]+spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BB)])*1.0f) / ((spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AB)]+spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BA)]) * 1.0f);
				fprintf(stdout, "\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%ld\t%.2g%%\t%.1f",
					(spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AA)]+spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BB)]), (spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AA)]+spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BB)])*100.0/totalInBins,
					(spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AB)]+spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BA)]), (spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AB)]+spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BA)])*100.0/totalInBins,
					spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AA)], spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AA)]*100.0/totalInBins,
					spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BB)], spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BB)]*100.0/totalInBins,
					spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AB)], spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_AB)]*100.0/totalInBins,
					spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BA)], spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_BA)]*100.0/totalInBins,
					ratio);
			}
			if (2==(2&linkerTypeSeen)) fprintf(stdout, "\t%ld\t%.2g%%",
				spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_NONE)], spans.bins[LINKER_OFFSET_IN_BIN(i,CPU_SPAN_LINKER_NONE)]*100.0/totalInBins);
			fprintf(stdout, "\n");
		}
		
		fprintf(stdout, ">>END\n");
	}
	
	// END - report the span statistics
	terminate_Span_Histogram (&spans);
	
	if (opt->toGroup) {
		terminate_CPSpan_Outputs(outfds);
	}
	
	span_opt_terminate(opt);
	free(opt);
	
	return ret;
}


