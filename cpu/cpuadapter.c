#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h> // for multi-threading
#include "utils.h"
#include "kseq.h"
#include "kstring.h"
#include "bwa.h"
KSEQ_DECLARE(gzFile)

#include "bwamem.h" // for MEM_F_PE
//WCH
//extern unsigned char nst_nt4_table[256];

#include "ssw.h"

#include "cpuadapter.h"
#include "cputag.h"

extern double G_t_real;

// available in >=v1.2.4
//#define USE_LARGER_GZIP_BUFFER
#define G_GZIP_BUFFER_SIZE (128*1024)

#define CPU_ADAPTER_OUTPUT_PREFIX     "cpadapters"
#define CPU_ADAPTER_OUTPUT_PAIRED          0x0000
#define CPU_ADAPTER_OUTPUT_UNPAIRED        0x0001
#define CPU_ADAPTER_OUTPUT_PAIRED_SHORT    0x0002
#define CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT  0x0003
#define CPU_ADAPTER_OUTPUT_COUNT           0x0004

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

typedef struct {
	int32_t adapter_type;
	s_align *paa;
} adapter_align_t;

typedef struct {
	int32_t adapter_type;
	int32_t adapter_bufsize;
	int8_t *adapter_num;
} adapter_t;

typedef struct {
	int classification;
	kstring_t fq;
} cpu_adapter_fq_t;

typedef struct {
	kstring_t adapter5;
	kstring_t adapter3;
	kstring_t rc_adapter5;
	kstring_t rc_adapter3;
	
	int8_t  aln_flag;
	int32_t match;
	int32_t mismatch;
	int32_t gap_open;
	int32_t gap_extension;
	
	int32_t nNT;
	int8_t* mata;
	
	int32_t path;
	int32_t filter;
	
	adapter_t align_adapter5;
	adapter_t align_adapter3;
	adapter_t align_rcAdapter5;
	adapter_t align_rcAdapter3;
	
	int32_t filterAdapter5;
	int32_t filterAdapter3;
#ifdef POSITIVE_SCORE_MATRIX
#else
	int32_t minAD5Overlap;
	int32_t filterAD5Overlap;
	int32_t minAD3Overlap;
	int32_t filterAD3Overlap;
#endif
	
	int minQuality;
	int minQualityLen;
	int minReadLen;
	
	int chunk_size;
	int n_threads;
	int T;
	int flag;               // see MEM_F_* macros
	int cpuflag;               // see MEM_F_* macros
	
	kstring_t outputPrefix;
} adapter_opt_t;

void initAdapter(adapter_opt_t *o)
{
	ks_set(o->adapter5.s, &(o->rc_adapter5));
	ks_set(o->adapter3.s, &(o->rc_adapter3));
	reverse_comple(o->adapter5.s, o->rc_adapter5.s, o->adapter5.l);
	reverse_comple(o->adapter3.s, o->rc_adapter3.s, o->adapter3.l);
}

void init_AdapterScoreMatrix (int8_t *mat, int32_t match, int32_t mismatch)
{
	int32_t l, m, k;
	// initialize scoring matrix
	for (l = k = 0; LIKELY(l < 4); ++l) {
		for (m = 0; LIKELY(m < 4); ++m) mat[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
		mat[k++] = 0; // ambiguous base
	}
	for (m = 0; LIKELY(m < 5); ++m) mat[k++] = 0;
}

void init_ssw_Adapter(adapter_t *assw, const kstring_t *ant, int32_t adapterType)
{
	assw->adapter_type = adapterType;
	assw->adapter_bufsize = (int32_t) ant->l;
	assw->adapter_num = (int8_t*)realloc(assw->adapter_num, assw->adapter_bufsize);
	int32_t m;
	for (m = 0; m < assw->adapter_bufsize; ++m) assw->adapter_num[m] = nt_table[(int)ant->s[m]];
}

void init_AlignmentAdapter(adapter_opt_t *o)
{
	initAdapter(o);
	
	init_ssw_Adapter(&(o->align_adapter5), &(o->adapter5), ADAPTER_5P);
	init_ssw_Adapter(&(o->align_adapter3), &(o->adapter3), ADAPTER_3P);
	init_ssw_Adapter(&(o->align_rcAdapter5), &(o->rc_adapter5), ADAPTER_5Prc);
	init_ssw_Adapter(&(o->align_rcAdapter3), &(o->rc_adapter3), ADAPTER_3Prc);
}

void initAdapterAlignmentParameters(adapter_opt_t *o)
{
	// re-set up the scoring matrix (user override from command line)
	init_AdapterScoreMatrix(o->mata, o->match, o->mismatch);
	
	init_AlignmentAdapter(o);
}

adapter_opt_t *adapter_opt_init()
{
	adapter_opt_t *o;
	o = calloc(1, sizeof(adapter_opt_t));
	
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->T = 21;
	o->flag = 0;
	o->cpuflag = CPU_FASTQ;
	o->cpuflag &= ~CPU_COMPRESS;
	
	memset(&(o->adapter5), 0, sizeof(kstring_t));
	memset(&(o->adapter3), 0, sizeof(kstring_t));
	memset(&(o->rc_adapter5), 0, sizeof(kstring_t));
	memset(&(o->rc_adapter3), 0, sizeof(kstring_t));
	
	ks_set(G_NEXTERA_5PAdapter, &(o->adapter5));
	ks_set(G_NEXTERA_3PAdapter, &(o->adapter3));
	
#ifdef POSITIVE_SCORE_MATRIX
	o->match = 1;
	o->mismatch = 1;
	o->gap_open = 6;
	o->gap_extension = 1;
#else
	o->match = 1;
	o->mismatch = 1;
	o->gap_open = 1;
	o->gap_extension = 1;
#endif
	
	o->nNT = 4+1; // A, C, G, T, N
	o->mata = (int8_t*)calloc(o->nNT*o->nNT, sizeof(int8_t));
	init_AdapterScoreMatrix(o->mata, o->match, o->mismatch);
	
	o->aln_flag = 1;
	o->path = 0; // TODO
	o->filter = 0; // TODO
	
	memset(&(o->align_adapter5), 0, sizeof(adapter_t));
	memset(&(o->align_adapter3), 0, sizeof(adapter_t));
	memset(&(o->align_rcAdapter5), 0, sizeof(adapter_t));
	memset(&(o->align_rcAdapter3), 0, sizeof(adapter_t));
	
#ifdef POSITIVE_SCORE_MATRIX
	o->filterAdapter5 = 21; // TODO: to determine based on the scroing matrix
	o->filterAdapter3 = 21; // TODO: to determine based on the scroing matrix
#else
	o->filterAdapter5 = 27; // TODO: to determine based on the scroing matrix
	o->filterAdapter3 = 28; // TODO: to determine based on the scroing matrix
#if 1
	o->minAD5Overlap = 12;
	o->filterAD5Overlap = (o->minAD5Overlap) * o->match - o->mismatch; // TODO: need to re-calculate
	o->minAD3Overlap = 12;
	o->filterAD3Overlap = (o->minAD3Overlap) * o->match - o->mismatch; // TODO: need to re-calculate
#else
	o->minAD5Overlap = 22;
	o->filterAD5Overlap = (o->minAD5Overlap) * o->match - o->mismatch; // TODO: need to re-calculate
	o->minAD3Overlap = 22;
	o->filterAD3Overlap = (o->minAD3Overlap) * o->match - o->mismatch; // TODO: need to re-calculate
#endif
#endif
	
	initAdapter(o);
	
	o->minQuality = 0;
	o->minQualityLen = 0;
	o->minReadLen = 18;
	
	ks_set(CPU_ADAPTER_OUTPUT_PREFIX, &(o->outputPrefix));
	
	return o;
}

void adapter_opt_terminate(adapter_opt_t *o)
{
	free(o->adapter5.s);
	free(o->adapter3.s);
	free(o->rc_adapter5.s);
	free(o->rc_adapter3.s);

	free(o->mata);
	
	free(o->align_adapter5.adapter_num);
	free(o->align_adapter3.adapter_num);
	free(o->align_rcAdapter5.adapter_num);
	free(o->align_rcAdapter3.adapter_num);
	
	free(o->outputPrefix.s);
}

// TODO: handle single-end and paired-end reads
static inline int init_CPAdapters_Ouputs (int PEreads, int flag, const char *outPrefix, gzFile *outfds) {
	const char *mode = (flag & CPU_COMPRESS) ? "w6" : "wT"; //"w0";
	const char *outSuffix = (flag & CPU_COMPRESS) ? ".gz" : "";
	
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.%s.fastq%s", outPrefix, (0!=PEreads ? "pe" : "se"), outSuffix);
		outfds[CPU_ADAPTER_OUTPUT_PAIRED] = gzopen (filename.s, mode);
		if (outfds[CPU_ADAPTER_OUTPUT_PAIRED] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_ADAPTER_OUTPUT_PAIRED+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_ADAPTER_OUTPUT_PAIRED], G_GZIP_BUFFER_SIZE);
	}
	
	if (0!=PEreads) {
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.pe.unpaired.fastq%s", outPrefix, outSuffix);
		outfds[CPU_ADAPTER_OUTPUT_UNPAIRED] = gzopen (filename.s, mode);
		if (outfds[CPU_ADAPTER_OUTPUT_UNPAIRED] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_ADAPTER_OUTPUT_UNPAIRED+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_ADAPTER_OUTPUT_UNPAIRED], G_GZIP_BUFFER_SIZE);
	}
	
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.%s.tooshort.fastq%s", outPrefix, (0!=PEreads ? "pe" : "se"), outSuffix);
		outfds[CPU_ADAPTER_OUTPUT_PAIRED_SHORT] = gzopen (filename.s, mode);
		if (outfds[CPU_ADAPTER_OUTPUT_PAIRED_SHORT] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_ADAPTER_OUTPUT_PAIRED_SHORT+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_ADAPTER_OUTPUT_PAIRED_SHORT], G_GZIP_BUFFER_SIZE);
	}
	
	if (0!=PEreads) {
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.pe.unpaired.tooshort.fastq%s", outPrefix, outSuffix);
		outfds[CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT] = gzopen (filename.s, mode);
		if (outfds[CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT], G_GZIP_BUFFER_SIZE);
	}
	
	return 0;
}

static inline int terminate_CPAdapters_Ouputs (gzFile *outfds) {
	int nFailed = 0;
	if (outfds[CPU_ADAPTER_OUTPUT_PAIRED]) {
		if  (gzclose(outfds[CPU_ADAPTER_OUTPUT_PAIRED])) {
			fprintf(stderr, "[E::%s] fail to close CPAdapters-1st gzip stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_ADAPTER_OUTPUT_UNPAIRED]) {
		if  (gzclose(outfds[CPU_ADAPTER_OUTPUT_UNPAIRED])) {
			fprintf(stderr, "[E::%s] fail to close CPAdapters-2nd gzip stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_ADAPTER_OUTPUT_PAIRED_SHORT]) {
		if  (gzclose(outfds[CPU_ADAPTER_OUTPUT_PAIRED_SHORT])) {
			fprintf(stderr, "[E::%s] fail to close CPAdapters-3rd gzip stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT]) {
		if  (gzclose(outfds[CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT])) {
			fprintf(stderr, "[E::%s] fail to close CPAdapters-4th gzip stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	
	return nFailed;
}

typedef struct {
	const adapter_opt_t *opt;
	bseq1_t *seqs;
	adapter_align_t *aas;
	int *read3primes;
	cpu_adapter_fq_t *fqs;
	int64_t n_processed;
} worker_t;

static inline int nexteraRCAdapter3Overlap3prime(const s_align *pa, const adapter_opt_t *o, const bseq1_t *seqs)
{
	if (!pa) return 0;
	
	if (pa->score1<o->filterAD3Overlap) return 0; // did not reach the minimum score
	
	if (0!=pa->ref_begin1) return 0; // the adapter's 5' end was not the start of match!
	
	// TODO: we cannot test end of read as we have not quality trimmed yet
	// TODO: so, for now, we used quality-checking as valid sign that the read is ending
	if (pa->read_end1!=seqs->l_seq-1) return 0; // did not end at 3'-end of read
	
	// TODO: need an array for the score thresholding!!! based on 10% error rate!
	int refSpan = pa->ref_end1 - pa->ref_begin1 + 1;
	if (refSpan>=28) {
		if (refSpan-pa->score1<=6) {
			return 1;
		}
#if 1
	} else if (refSpan>=18) {
		if (refSpan-pa->score1<=4) {
			return 1;
		}
	} else if (refSpan>=12) {
		if (refSpan-pa->score1<=2) {
			return 1;
		}
#else
	} else if (refSpan>=24) {
		if (refSpan-pa->score1<=4) {
			return 1;
		}
	} else if (refSpan>=23) {
		if (refSpan-pa->score1<=2) {
			return 1;
		}
	} else if (refSpan>=22) {
		if (refSpan-pa->score1==0) {
			return 1;
		}
#endif
	} else {
		return 0;
	}
	return 0;
}

static inline int nexteraRCAdapter5Overlap3prime(const s_align *pa, const adapter_opt_t *o, const bseq1_t *seqs)
{
	if (!pa) return 0;
	
	if (pa->score1<o->filterAD5Overlap) return 0; // did not reach the minimum score
	
	if (0!=pa->ref_begin1) return 0; // the adapter's 5' end was not the start of match!
	
	// TODO: we cannot test end of read as we have not quality trimmed yet
	// TODO: so, for now, we used quality-checking as valid sign that the read is ending
	if (pa->read_end1!=seqs->l_seq-1) return 0; // did not end at 3'-end of read
	
	// TODO: need an array for the score thresholding!!! based on 10% error rate!
	int refSpan = pa->ref_end1 - pa->ref_begin1 + 1;
	if (refSpan>=28) {
		if (refSpan-pa->score1<=6) {
			return 1;
		}
#if 1
	} else if (refSpan>=18) {
		if (refSpan-pa->score1<=4) {
			return 1;
		}
	} else if (refSpan>=12) {
		if (refSpan-pa->score1<=2) {
			return 1;
		}
#else
	} else if (refSpan>=24) {
		if (refSpan-pa->score1<=4) {
			return 1;
		}
	} else if (refSpan>=23) {
		if (refSpan-pa->score1<=2) {
			return 1;
		}
	} else if (refSpan>=22) {
		if (refSpan-pa->score1==0) {
			return 1;
		}
#endif
	} else {
		return 0;
	}
	return 0;
}

static void adapter_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		//mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a);
		//mem_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
		//free(w->regs[i].a);
	} else {
		// Process Read/1
		adapter_align_t *pAARec = &(w->aas[i<<1|0]);
		pAARec->adapter_type = ADAPTER_NONE; pAARec->paa = 0;
		
		s_profile *p1 = 0;
		int32_t readLen = w->seqs[i<<1|0].l_seq;
		int32_t maskLen = readLen / 2;
		//int32_t maskLen = w->opt->halfLinkerA.l / 2;
		int8_t *num = (int8_t*)malloc(readLen+1);
		int32_t m;
		for (m = 0; m < readLen; ++m) num[m] = nt_table[(int)w->seqs[i<<1|0].seq[m]];
		p1 = ssw_init(num, readLen, w->opt->mata, w->opt->nNT, 2);
		
		// TODO: this should be saved and processed centrally by the caller
		// test linkers alignment
		s_align *resNTrc3p = 0; const adapter_t *adapterNTrc3p = &(w->opt->align_rcAdapter3);
		resNTrc3p = ssw_align (p1, adapterNTrc3p->adapter_num, adapterNTrc3p->adapter_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
		pAARec->paa = resNTrc3p;
		if (resNTrc3p && resNTrc3p->score1>=w->opt->filterAdapter3) {
			pAARec->adapter_type = ADAPTER_3Prc;
		} else {
			if (0!=nexteraRCAdapter3Overlap3prime(pAARec->paa, w->opt, &w->seqs[i<<1|0])) {
				// okie, we have a prefix suffix case
				pAARec->adapter_type = ADAPTER_3Prc;
			}
		}
		init_destroy(p1);
		free(num);
		
		w->read3primes[i<<1|0] = trim_read_3prime(w->opt->minQuality, w->opt->minQualityLen, &(w->seqs[i<<1|0]));
		
		// Process Read/2
		pAARec = &(w->aas[i<<1|1]);
		pAARec->adapter_type = ADAPTER_NONE; pAARec->paa = 0;
		
		s_profile *p2 = 0;
		readLen = w->seqs[i<<1|1].l_seq;
		maskLen = readLen / 2;
		//int32_t maskLen = w->opt->halfLinkerA.l / 2;
		num = (int8_t*)malloc(readLen+1);
		for (m = 0; m < readLen; ++m) num[m] = nt_table[(int)w->seqs[i<<1|1].seq[m]];
		p2 = ssw_init(num, readLen, w->opt->mata, w->opt->nNT, 2);
		
		// TODO: this should be saved and processed centrally by the caller
		// test linkers alignment
		s_align *resNTrc5p = 0; const adapter_t *adapterNTrc5p = &(w->opt->align_rcAdapter5);
		resNTrc5p = ssw_align (p2, adapterNTrc5p->adapter_num, adapterNTrc5p->adapter_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
		pAARec->paa = resNTrc5p;
		if (resNTrc5p && resNTrc5p->score1>=w->opt->filterAdapter5) {
			pAARec->adapter_type = ADAPTER_5Prc;
		} else {
			if (0!=nexteraRCAdapter5Overlap3prime(pAARec->paa, w->opt, &w->seqs[i<<1|0])) {
				// okie, we have a prefix suffix case
				pAARec->adapter_type = ADAPTER_5Prc;
			}
		}
		init_destroy(p2);
		free(num);
		
		w->read3primes[i<<1|1] = trim_read_3prime(w->opt->minQuality, w->opt->minQualityLen, &(w->seqs[i<<1|1]));
		
		{
			adapter_align_t *pAARec;
			pAARec = &(w->aas[i<<1|0]);
			int len1 = (ADAPTER_3Prc==pAARec->adapter_type) ? pAARec->paa->read_begin1 : w->seqs[i<<1|0].l_seq;
			if (w->read3primes[i<<1|0]<len1) len1 = w->read3primes[i<<1|0];
			
			pAARec = &(w->aas[i<<1|1]);
			int len2 = (ADAPTER_5Prc==pAARec->adapter_type) ? pAARec->paa->read_begin1 : w->seqs[i<<1|1].l_seq;
			if (w->read3primes[i<<1|1]<len2) len2 = w->read3primes[i<<1|1];
			
			kstring_t *pfq = &(w->fqs[i<<1|0].fq);
			kputc('@', pfq); kputs(w->seqs[i<<1|0].name, pfq); kputs("/1", pfq); kputc('\n', pfq);
			kputsn(w->seqs[i<<1|0].seq, len1, pfq); kputc('\n', pfq);
			kputs("+\n", pfq);
			kputsn(w->seqs[i<<1|0].qual, len1, pfq); kputc('\n', pfq);
		
			pfq = &(w->fqs[i<<1|1].fq);
			kputc('@', pfq); kputs(w->seqs[i<<1|1].name, pfq); kputs("/2", pfq); kputc('\n', pfq);
			kputsn(w->seqs[i<<1|1].seq, len1, pfq); kputc('\n', pfq);
			kputs("+\n", pfq);
			kputsn(w->seqs[i<<1|1].qual, len1, pfq); kputc('\n', pfq);

			//classification
			if (len1>w->opt->minReadLen) {
				if (len2>w->opt->minReadLen) {
					w->fqs[i<<1|0].classification = w->fqs[i<<1|1].classification = CPU_ADAPTER_OUTPUT_PAIRED;
				} else {
					w->fqs[i<<1|0].classification = CPU_ADAPTER_OUTPUT_UNPAIRED;
					w->fqs[i<<1|1].classification = CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT;
				}
			} else {
				if (len2>w->opt->minReadLen) {
					w->fqs[i<<1|0].classification = CPU_ADAPTER_OUTPUT_UNPAIRED_SHORT;
					w->fqs[i<<1|1].classification = CPU_ADAPTER_OUTPUT_UNPAIRED;
				} else {
					w->fqs[i<<1|0].classification = CPU_ADAPTER_OUTPUT_PAIRED_SHORT;
					w->fqs[i<<1|1].classification = CPU_ADAPTER_OUTPUT_PAIRED_SHORT;
				}
			}
		}
	}
}

void adapter_process_seqs(const adapter_opt_t *opt, bseq1_t *seqs, adapter_align_t *aas, int *read3primes, cpu_adapter_fq_t *fqs, int64_t n_processed, int n)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	w.opt = opt;
	w.seqs = seqs;
	w.aas = aas;
	w.read3primes = read3primes;
	w.fqs = fqs;
	//w.n_processed = n_processed;
	w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	kt_for(opt->n_threads, adapter_worker, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
}

int main_adapter(int argc, char *argv[])
{
	adapter_opt_t *opt;
	char *p;
	int fd, fd2, i, c, n;
	int copy_comment = 0;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	void *ko = 0, *ko2 = 0;
	int64_t n_processed = 0;
	
	double t_diff;
	
	gzFile outfds[CPU_ADAPTER_OUTPUT_COUNT] = {0};

	opt = adapter_opt_init();
	while ((c = getopt(argc, argv, "NCSQp5:3:t:T:q:r:O:a")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'T') opt->T = atoi(optarg), opt->T = opt->T > 0 ? opt->T : 21;
		else if (c == '5') { ks_set(optarg, &(opt->adapter5)); }
		else if (c == '4') { ks_set(optarg, &(opt->adapter3)); }
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'N') opt->cpuflag |= CPU_NAME;
		else if (c == 'S') opt->cpuflag |= CPU_SEQ;
		else if (c == 'Q') opt->cpuflag |= CPU_QUALITY;
		else if (c == 'C') opt->cpuflag |= CPU_COMPRESS;
		else if (c == 'F') opt->cpuflag &= ~CPU_FASTQ;
		else if (c == 'q') {
			opt->minQuality = (int32_t) strtol(optarg, &p, 10);
			if (0!=*p && ispunct(*p) && isdigit(p[1])) opt->minQualityLen = (int32_t) strtol(p+1, &p, 10);
			opt->minQuality = opt->minQuality > 0 ? opt->minQuality : 0;
			opt->minQualityLen = opt->minQualityLen > 0 ? opt->minQualityLen : 0;
		}
		else if (c == 'r') opt->minReadLen = atoi(optarg), opt->minReadLen = opt->minReadLen > 0 ? opt->minReadLen : 0;
		else if (c == 'O') { ks_set(optarg, &(opt->outputPrefix)); }
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else {
			adapter_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu adapter [options] <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -5 STR        5' adapter [%s]\n", opt->adapter5.s);
		fprintf(stderr, "       -3 STR        3' adapter [%s]\n", opt->adapter3.s);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -O STR        output prefix [%s]\n", opt->outputPrefix.s);
		fprintf(stderr, "       -C            compress .fastq\n");
		fprintf(stderr, "       -q INT[,INT]  quality threshold for read [%d], trimming down to [%dbp]\n", opt->minQuality, opt->minQualityLen);
		fprintf(stderr, "       -r INT        min. read len to output [%dbp]\n", opt->minReadLen);
		fprintf(stderr, "       -a            out both paired-end and single end reads\n");
		fprintf(stderr, "       -F            skip .fastq writing\n");
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		//fprintf(stderr, "       -N            write read name\n");
		fprintf(stderr, "       -S            write sequence\n");
		fprintf(stderr, "       -Q            write quality\n");
		fprintf(stderr, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		adapter_opt_terminate(opt);
		free(opt);
		return 1;
	}

	ko = kopen(argv[optind], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind]);
		
		adapter_opt_terminate(opt);
		free(opt);
		return 1;
	}
	fp = gzdopen(fd, "r");
#ifdef USE_LARGER_GZIP_BUFFER
	gzbuffer(fp, G_GZIP_BUFFER_SIZE);
#endif
	ks = kseq_init(fp);
	if (optind + 1 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file will be ignored.\n", __func__);
		} else {
			ko2 = kopen(argv[optind + 1], &fd2);
			if (ko2 == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
				
				adapter_opt_terminate(opt);
				free(opt);
				return 1;
			}
			fp2 = gzdopen(fd2, "r");
#ifdef USE_LARGER_GZIP_BUFFER
			gzbuffer(fp2, G_GZIP_BUFFER_SIZE);
#endif
			ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}
	
	initAdapter(opt);
	initAdapterAlignmentParameters(opt);

	int PEreads = (opt->flag&MEM_F_PE || 0!=ko2);
	if (opt->cpuflag & CPU_FASTQ) {
		if (init_CPAdapters_Ouputs(PEreads, opt->cpuflag, opt->outputPrefix.s, outfds)) {
			terminate_CPAdapters_Ouputs(outfds);
			
			adapter_opt_terminate(opt);
			free(opt);
			kseq_destroy(ks);
			err_gzclose(fp); kclose(ko);
			if (ks2) {
				kseq_destroy(ks2);
				err_gzclose(fp2); kclose(ko2);
			}
			
			return 1;
		}
	}
	
	// TODO: should initialize shared resources across threads
	
	while ((seqs = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
		int64_t size = 0;
		if ((opt->flag & MEM_F_PE) && (n&1) == 1) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] odd number of reads in the PE mode; last read dropped\n", __func__);
			n = n>>1<<1;
		}
		if (!copy_comment)
			for (i = 0; i < n; ++i) {
				free(seqs[i].comment); seqs[i].comment = 0;
			}
		for (i = 0; i < n; ++i) size += seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp)..\n", __func__, n, (long)size);
		
		// might need to pass statistics holder
		adapter_align_t *aas = calloc(n, sizeof(adapter_align_t));
		int *read3primes = calloc(n, sizeof(int));
		cpu_adapter_fq_t *fqs = calloc(n, sizeof(cpu_adapter_fq_t));
		adapter_process_seqs(opt, seqs, aas, read3primes, fqs, n_processed, n);
		
		// OUTPUT:
		if (opt->cpuflag & CPU_FASTQ) {
			for (i = 0; i < n; ++i) {
				gzwrite(outfds[fqs[i].classification], fqs[i].fq.s, fqs[i].fq.l);
			}
		}

		if (0!=(opt->cpuflag & CPU_DEBUG_MASK)) {
			// OLD code for text based output (debugging purposes/just counting)
			for (i = 0; i < n; i+=2) {
				fprintf(stdout, "%s", seqs[i].name);
				adapter_align_t *pAARec = &(aas[i]);
				if (ADAPTER_NONE==pAARec->adapter_type) {
					fprintf(stdout, "\t-1\t-1\t-1\t-1\t-1\t.");
				} else {
					fprintf(stdout, "\t%d\t%d\t%d\t%d\t%d\t%s", pAARec->paa->ref_begin1, pAARec->paa->ref_end1, pAARec->paa->read_begin1, pAARec->paa->read_end1, pAARec->paa->score1, adapterTypeToStr(pAARec->adapter_type));
					align_destroy(pAARec->paa);
				}
				pAARec++;
				if (ADAPTER_NONE==pAARec->adapter_type) {
					fprintf(stdout, "\t-1\t-1\t-1\t-1\t-1\t.");
				} else {
					fprintf(stdout, "\t%d\t%d\t%d\t%d\t%d\t%s", pAARec->paa->ref_begin1, pAARec->paa->ref_end1, pAARec->paa->read_begin1, pAARec->paa->read_end1, pAARec->paa->score1, adapterTypeToStr(pAARec->adapter_type));
					align_destroy(pAARec->paa);
				}
				
				// TODO: use a switch to decide when to output the sequence
				if (opt->cpuflag & CPU_SEQ) {
					fprintf(stdout, "\t%s\t%s", seqs[i].seq, seqs[i+1].seq);
				}
				if (opt->cpuflag & CPU_QUALITY) {
					fprintf(stdout, "\t%s\t%s", seqs[i].qual, seqs[i+1].qual);
				}
				fprintf(stdout, "\n");
			}
		}
		// END - OUTPUT
		
		free(aas);
		free(read3primes);
		
		for (i = 0; i < n; ++i) {
			//err_fputs(seqs[i].sam, stdout);
			free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual);
			//free(seqs[i].sam);
			free(fqs[i].fq.s);
		}
		free(seqs);
		free(fqs);
		
		n_processed += n;
		if (bwa_verbose >= 3) {
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s] %lld sequences processed, %.0f seq/sec, %.2f min..\n", __func__, n_processed, 1.0*n_processed/t_diff, t_diff/60.0);
		}
	}
	
	terminate_CPAdapters_Ouputs(outfds);
	
	adapter_opt_terminate(opt);
	free(opt);
	kseq_destroy(ks);
	err_gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	
	return 0;
}


