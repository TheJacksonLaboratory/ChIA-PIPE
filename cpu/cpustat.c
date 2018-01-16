// TODO:
// dynamics allocation of read length

#include "zlib.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h> // for multi-threading
#include <errno.h>
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
#include "cpulinker.h"
#include "cputag.h"
#include "cpustag.h"

extern double G_t_real;
extern const char *CPU_version;

#define G_GZIP_LINE_BUFFER (8*1024)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

// TODO: clean up tags at the end
uint8_t catOrders[CPU_OUTPUT_COUNT]={
CPU_OUTPUT_FL_NC_2TAGS,
CPU_OUTPUT_FL_NC_1TAG,
CPU_OUTPUT_FL_C_2TAGS,
CPU_OUTPUT_FL_C_1TAG,
CPU_OUTPUT_NONE,
CPU_OUTPUT_CONFLICT,
CPU_OUTPUT_HL,
CPU_OUTPUT_TIED};

#define REPORT_LINKER_COUNT 8
uint8_t linkerOrders[REPORT_LINKER_COUNT]={
	LINKER_A_PRIV | LINKER_RC_A_PRIV,
	LINKER_B_PRIV | LINKER_RC_B_PRIV,
	LINKER_A_PRIV | LINKER_RC_B_PRIV,
	LINKER_B_PRIV | LINKER_RC_A_PRIV,
	LINKER_A_PRIV,
	LINKER_RC_A_PRIV,
	LINKER_B_PRIV,
	LINKER_RC_B_PRIV
};
char *linkerLabels[REPORT_LINKER_COUNT]={
	"Aa",
	"Bb",
	"Ab",
	"Ba",
	"A",
	"a",
	"B",
	"b"
};

uint8_t slCatOrders[CPU_OUTPUT_SL_COUNT]={
CPU_OUTPUT_SL_2TAGS,
CPU_OUTPUT_SL_1TAG,
CPU_OUTPUT_NONE,
CPU_OUTPUT_CONFLICT,
CPU_OUTPUT_TIED};

#define REPORT_SLINKER_COUNT 2
uint8_t slinkerOrders[REPORT_SLINKER_COUNT]={
	LINKER_SINGLE_PRIV,
	LINKER_RC_SINGLE_PRIV
};
char *slinkerLabels[REPORT_SLINKER_COUNT]={
	"SL",
	"ls"
};

static char outputClassUserBuffer[11];
static const char *outputClassToUserStr(int16_t t) {
	if (CPU_OUTPUT_NONE==t)             return "No Linker";
	else if (CPU_OUTPUT_TIED==t)        return ">1 Linker in read";
	else if (CPU_OUTPUT_CONFLICT==t)    return "Conflict";
	else if (CPU_OUTPUT_HL==t)          return "Half Linker (R1/R2)";
	else if (CPU_OUTPUT_FL_C_2TAGS==t)  return "Full Linker 2 tags (AB/BA)";
	else if (CPU_OUTPUT_FL_C_1TAG==t)   return "Full Linker 1 tag (AB/BA)";
	else if (CPU_OUTPUT_FL_NC_2TAGS==t) return "Full Linker 2 tags (AA/BB)";
	else if (CPU_OUTPUT_FL_NC_1TAG==t)  return "Full Linker 1 tag (AA/BB)";
	else {
		//return "?";
		sprintf(outputClassUserBuffer, "%#x", t);
		return outputClassUserBuffer;
	}
}

static const char *outputSLClassToUserStr(int16_t t) {
	if (CPU_OUTPUT_NONE==t)          return "No Linker";
	else if (CPU_OUTPUT_TIED==t)     return ">1 Linker in read";
	else if (CPU_OUTPUT_CONFLICT==t) return "Conflict";
	else if (CPU_OUTPUT_SL_2TAGS==t) return "Single Linker 2 tags (SL/ls)";
	else if (CPU_OUTPUT_SL_1TAG==t)  return "Single Linker 1 tag (SL/ls)";
	else {
		//return "?";
		sprintf(outputClassUserBuffer, "%#x", t);
		return outputClassUserBuffer;
	}
}

typedef struct {
	int chunk_size;
	int n_threads;
	int T;
	int flag;               // see MEM_F_* macros
	
	int cpuflag;               // see MEM_F_* macros
	int minTagLen;
	int linkerFlag;
	int nTagFamily;
	int classType;
} stat_opt_t;

// TODO: handle larger capacity?
// TODO: linker count for DiLinker and SingleLinker is different
#define LINKER_COUNT 16
#if 0
#define MAX_READ_LEN (151+1)
typedef struct {
	uint32_t classification[CPU_OUTPUT_COUNT];
	uint32_t classTags[CPU_OUTPUT_COUNT];
	//uint32_t readLen;
	// TODO: to get the proper length of the possible read position distribution
	uint32_t endAdapter[2][MAX_READ_LEN];
	uint32_t linker[2][LINKER_COUNT][MAX_READ_LEN];
	uint32_t taglen[2][2][MAX_READ_LEN];
} stat_classification_t;
#else
#define MAX_READ_LEN (151)
typedef struct {
	uint32_t classification[CPU_OUTPUT_COUNT];
	uint32_t classTags[CPU_OUTPUT_COUNT];
	//uint32_t readLen;
	// TODO: to get the proper length of the possible read position distribution
	uint32_t *endAdapter[2];
	uint32_t *linker[2][LINKER_COUNT];
	uint32_t *taglen[2][2];
} stat_classification_t;
#endif

// NOTE: not with re-entrant; memory loss
void stat_classification_init(stat_classification_t *s, int maxReadLength)
{
	int i;
	int nSize;
	
	memset(s, 0, sizeof(stat_classification_t));
	
	nSize = maxReadLength + 1;
	
	s->endAdapter[0] = calloc(nSize, sizeof(uint32_t));
	s->endAdapter[1] = calloc(nSize, sizeof(uint32_t));
	
	for(i=0; i<LINKER_COUNT; ++i) {
		s->linker[0][i] = calloc(nSize, sizeof(uint32_t));
		s->linker[1][i] = calloc(nSize, sizeof(uint32_t));
	}
	
	s->taglen[0][0] = calloc(nSize, sizeof(uint32_t));
	s->taglen[0][1] = calloc(nSize, sizeof(uint32_t));
	s->taglen[1][0] = calloc(nSize, sizeof(uint32_t));
	s->taglen[1][1] = calloc(nSize, sizeof(uint32_t));
}

void stat_classification_terminate(stat_classification_t *s)
{
	int i;
	
	free(s->endAdapter[0]);
	free(s->endAdapter[1]);
	
	for(i=0; i<LINKER_COUNT; ++i) {
		free(s->linker[0][i]);
		free(s->linker[1][i]);
	}
	
	free(s->taglen[0][0]);
	free(s->taglen[0][1]);
	free(s->taglen[1][0]);
	free(s->taglen[1][1]);
}

stat_opt_t *stat_opt_init()
{
	stat_opt_t *o;
	o = calloc(1, sizeof(stat_opt_t));
	
	o->chunk_size = 1000000;
	o->n_threads = 1;
	o->T = 21;
	o->flag = 0;
	
	// TODO: minimal tag length for fastq generation!
	o->minTagLen = 20;
	
	o->linkerFlag = 2;
	o->nTagFamily = CPU_TAGTYPE_DILINKER;
	o->classType = -1;
	
	return o;
}

void stat_opt_terminate(stat_opt_t *o)
{
}

typedef struct {
	const stat_opt_t *opt;
	int64_t n_processed;
	
	cpu_record_t *cpus;
	stat_classification_t *stat_classes;
	int *maxReadLengths;
} worker_t;

static void stat_worker_maxReadLength (void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (w->cpus[i].readLen>w->maxReadLengths[tid]) w->maxReadLengths[tid] = w->cpus[i].readLen;
	} else {
		// read length determination
		if (w->cpus[i<<1|0].readLen>w->maxReadLengths[tid]) w->maxReadLengths[tid] = w->cpus[i<<1|0].readLen;
		if (w->cpus[i<<1|1].readLen>w->maxReadLengths[tid]) w->maxReadLengths[tid] = w->cpus[i<<1|1].readLen;
	}
}

void stat_process_maxReadLen(const stat_opt_t *opt, cpu_record_t *cpus, int *maxReadLengths, int64_t n_processed, int n, void (*stat_worker)(void*,int,int), int *maxReadLength)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	int i;
	worker_t w;
	
	w.opt = opt;
	w.cpus = cpus;
	w.stat_classes = NULL;
	w.maxReadLengths = maxReadLengths;
	//w.n_processed = n_processed;
	w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	kt_for(opt->n_threads, stat_worker, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
	
	// accummulate the summary
	for (i=0; i<opt->n_threads; ++i)
		if (maxReadLengths[i]>*maxReadLength) *maxReadLength = maxReadLengths[i];
}

// TODO: clean up the code as the computation end is the same
static void stat_worker_sl(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Not implemented yet.\n", __func__);
	} else {
		processPairedSTag (i + w->n_processed + 1,
						   &(w->cpus[i<<1|0]), NULL, NULL,
						   &(w->cpus[i<<1|1]), NULL, NULL,
						   w->opt->minTagLen, 0, w->opt->nTagFamily);
		
		// classification is for a read-pair, so we do not process [i<<1|1]
		w->stat_classes[tid].classification[w->cpus[i<<1|0].classification]++;
		w->stat_classes[tid].classTags[w->cpus[i<<1|0].classification] += w->cpus[i<<1|0].tags_n;
		w->stat_classes[tid].classTags[w->cpus[i<<1|0].classification] += w->cpus[i<<1|1].tags_n;
		
		if (-1==w->opt->classType || w->opt->classType==w->cpus[i<<1|0].classification) {
			// work on the end-adapter position
			if (ADAPTER_NONE!=w->cpus[i<<1|0].adapter_type)
				w->stat_classes[tid].endAdapter[0][w->cpus[i<<1|0].adapter_read_begin1]++;
			if (ADAPTER_NONE!=w->cpus[i<<1|1].adapter_type)
				w->stat_classes[tid].endAdapter[1][w->cpus[i<<1|1].adapter_read_begin1]++;
			
			// work on the linker position
			if (LINKER_NONE!=w->cpus[i<<1|0].linker_type)
				w->stat_classes[tid].linker[0][w->cpus[i<<1|0].linker_type & LINKER_MASK][w->cpus[i<<1|0].linker_read_begin1]++;
			if (LINKER_NONE!=w->cpus[i<<1|1].linker_type)
				w->stat_classes[tid].linker[1][w->cpus[i<<1|1].linker_type & LINKER_MASK][w->cpus[i<<1|1].linker_read_begin1]++;
			
			// work on the taglen
			if (LINKER_NONE==w->cpus[i<<1|0].linker_type) {
				w->stat_classes[tid].taglen[0][0][w->cpus[i<<1|0].len[0]]++;
			} else {
				// either 1 (L/R) tag or 2 tags (both L&R)
				w->stat_classes[tid].taglen[0][0][w->cpus[i<<1|0].len[0]]++;
				w->stat_classes[tid].taglen[0][1][w->cpus[i<<1|0].len[1]]++;
			}
			if (LINKER_NONE==w->cpus[i<<1|1].linker_type) {
				w->stat_classes[tid].taglen[1][0][w->cpus[i<<1|1].len[0]]++;
			} else {
				// either 1 (L/R) tag or 2 tags (both L&R)
				w->stat_classes[tid].taglen[1][0][w->cpus[i<<1|1].len[0]]++;
				w->stat_classes[tid].taglen[1][1][w->cpus[i<<1|1].len[1]]++;
			}
		}
	}
}

static void stat_worker_dl(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		processSingleTag(i + w->n_processed + 1,
						 &(w->cpus[i]), NULL, NULL,
						 w->opt->minTagLen, 0, w->opt->nTagFamily);
		
		// classification is for a read-pair, so we do not process [i<<1|1]
		w->stat_classes[tid].classification[w->cpus[i].classification]++;
		w->stat_classes[tid].classTags[w->cpus[i].classification] += w->cpus[i].tags_n;
		w->stat_classes[tid].classTags[w->cpus[i].classification] += w->cpus[i].tags_n;
		
		if (-1==w->opt->classType || w->opt->classType==w->cpus[i].classification) {
			// work on the end-adapter position
			if (ADAPTER_NONE!=w->cpus[i].adapter_type)
				w->stat_classes[tid].endAdapter[0][w->cpus[i].adapter_read_begin1]++;
			
			// work on the linker position
			if (LINKER_NONE!=w->cpus[i].linker_type)
				w->stat_classes[tid].linker[0][w->cpus[i].linker_type & LINKER_MASK][w->cpus[i].linker_read_begin1]++;
			
			// work on the taglen
			if (LINKER_NONE==w->cpus[i].linker_type) {
				w->stat_classes[tid].taglen[0][0][w->cpus[i].len[0]]++;
			} else {
				// either 1 (L/R) tag or 2 tags (both L&R)
				w->stat_classes[tid].taglen[0][0][w->cpus[i].len[0]]++;
				w->stat_classes[tid].taglen[0][1][w->cpus[i].len[1]]++;
			}
		}
		
	} else {
		// TODO: not implemented yet
		//       as the parameters can be different from the final .fastq produced
		//       it is much better to re-process the classification
		//       This makes it possible to explore the data characteristics w/o re-running the tag detection and generation

		/*
		 // this is for dual-linker system
		 // this is for lmp; a kind of dual-linker system, just that we have only Aa, NO Bb
		 */
		
		processPairedTag (i + w->n_processed + 1,
						  &(w->cpus[i<<1|0]), NULL, NULL,
						  &(w->cpus[i<<1|1]), NULL, NULL,
						  w->opt->minTagLen, 0, w->opt->nTagFamily);
		
		// classification is for a read-pair, so we do not process [i<<1|1]
		w->stat_classes[tid].classification[w->cpus[i<<1|0].classification]++;
		w->stat_classes[tid].classTags[w->cpus[i<<1|0].classification] += w->cpus[i<<1|0].tags_n;
		w->stat_classes[tid].classTags[w->cpus[i<<1|0].classification] += w->cpus[i<<1|1].tags_n;

		
		if (-1==w->opt->classType || w->opt->classType==w->cpus[i<<1|0].classification) {
			// work on the end-adapter position
			if (ADAPTER_NONE!=w->cpus[i<<1|0].adapter_type)
				w->stat_classes[tid].endAdapter[0][w->cpus[i<<1|0].adapter_read_begin1]++;
			if (ADAPTER_NONE!=w->cpus[i<<1|1].adapter_type)
				w->stat_classes[tid].endAdapter[1][w->cpus[i<<1|1].adapter_read_begin1]++;
			
			// work on the linker position
			if (LINKER_NONE!=w->cpus[i<<1|0].linker_type)
				w->stat_classes[tid].linker[0][w->cpus[i<<1|0].linker_type & LINKER_MASK][w->cpus[i<<1|0].linker_read_begin1]++;
			if (LINKER_NONE!=w->cpus[i<<1|1].linker_type)
				w->stat_classes[tid].linker[1][w->cpus[i<<1|1].linker_type & LINKER_MASK][w->cpus[i<<1|1].linker_read_begin1]++;
			
			// work on the taglen
			if (LINKER_NONE==w->cpus[i<<1|0].linker_type) {
				w->stat_classes[tid].taglen[0][0][w->cpus[i<<1|0].len[0]]++;
			} else {
				// either 1 (L/R) tag or 2 tags (both L&R)
				w->stat_classes[tid].taglen[0][0][w->cpus[i<<1|0].len[0]]++;
				w->stat_classes[tid].taglen[0][1][w->cpus[i<<1|0].len[1]]++;
			}
			if (LINKER_NONE==w->cpus[i<<1|1].linker_type) {
				w->stat_classes[tid].taglen[1][0][w->cpus[i<<1|1].len[0]]++;
			} else {
				// either 1 (L/R) tag or 2 tags (both L&R)
				w->stat_classes[tid].taglen[1][0][w->cpus[i<<1|1].len[0]]++;
				w->stat_classes[tid].taglen[1][1][w->cpus[i<<1|1].len[1]]++;
			}
		}
	}
}

void stat_process_seqs(const stat_opt_t *opt, cpu_record_t *cpus, stat_classification_t *stat_classes, int64_t n_processed, int n, int maxReadLength, void (*stat_worker)(void*,int,int), stat_classification_t *summary)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	int i, j, k;
	worker_t w;
	
	w.opt = opt;
	w.cpus = cpus;
	w.stat_classes = stat_classes;
	w.maxReadLengths = NULL;
	//w.n_processed = n_processed;
	w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	kt_for(opt->n_threads, stat_worker, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions

	// accummulate the summary
	//for (i=0; i<CPU_OUTPUT_COUNT; ++i) {
	for (i=0; i<opt->n_threads; ++i) {
		//if (stat_classes[i].readLen>summary->readLen) summary->readLen = stat_classes[i].readLen;
		for (j=0; j<CPU_OUTPUT_COUNT; ++j) {
			summary->classification[j] += stat_classes[i].classification[j];
			summary->classTags[j] += stat_classes[i].classTags[j];
		}
		for (j=0; j<=maxReadLength; ++j) {
			summary->endAdapter[0][j] += stat_classes[i].endAdapter[0][j];
			summary->endAdapter[1][j] += stat_classes[i].endAdapter[1][j];
		}
		for (j=0; j<=maxReadLength; ++j) {
			for (k=0; k<LINKER_COUNT; ++k) {
				summary->linker[0][k][j] += stat_classes[i].linker[0][k][j];
				summary->linker[1][k][j] += stat_classes[i].linker[1][k][j];
			}
		}
		for (j=0; j<=maxReadLength; ++j) {
			summary->taglen[0][0][j] += stat_classes[i].taglen[0][0][j];
			summary->taglen[0][1][j] += stat_classes[i].taglen[0][1][j];
			summary->taglen[1][0][j] += stat_classes[i].taglen[1][0][j];
			summary->taglen[1][1][j] += stat_classes[i].taglen[1][1][j];
		}
	}
}

cpu_record_t *stat_readcpu(int *nread, int n_, gzFile ci1_)
{
	cpu_record_t *cpus = 0;
	cpus = realloc(cpus, n_ * sizeof(cpu_record_t));
	int i = gzread(ci1_, cpus, n_ * sizeof(cpu_record_t));
	if (-1==i) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to read %d cpu record(s). -1 returned.\n", __func__, n_);
		free(cpus);
		*nread = 0;
		return NULL;
	}
	i /= sizeof(cpu_record_t);
	if (i != n_) {
		//if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to read %d cpu record(s). %d read instead.\n", __func__, n_, i);
	}
	*nread = i;
	if (0==i) {
		free(cpus); cpus=0;
	}
	return cpus;
}

int main_stat(int argc, char *argv[])
{
	int c, n, i, j;
	stat_opt_t *opt;
	cpu_record_t *cpus;
	int fdci;
	gzFile fpci = 0;
	void *ci = 0;
	int64_t n_processed = 0;
	
	double t_diff;

	opt = stat_opt_init();
	while ((c = getopt(argc, argv, "sdpt:T:c:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'T') opt->minTagLen = atoi(optarg), opt->minTagLen = opt->minTagLen >= 0 ? opt->minTagLen : 20;
		else if (c == 'c') opt->classType = atoi(optarg), opt->classType = -1==opt->classType ? opt->n_threads : (opt->classType>CPU_OUTPUT_SL_COUNT ? -1 : opt->classType);
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'd') opt->linkerFlag = 2;
		else if (c == 's') opt->linkerFlag = 1;
		else {
			stat_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu stat [options] <in.cpu> <in.fq>\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p         first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -T INT     min. tag len to be written to .fastq\n");
		fprintf(stderr, "       -d         di-linker system metrics\n");
		fprintf(stderr, "       -s         single-linker system metrics\n");
		fprintf(stderr, "       -c INT     classification to summarize\n");
		fprintf(stderr, "                  None (0), Tied (1), Conflict (2)\n");
		fprintf(stderr, "                  SL: 2 tags (3), 1 tag (4)\n");
		fprintf(stderr, "                  DL: HL (3), FLC 2 tags (4), FLC 1 tag (5), FLNC 2 tags (5), FLNC 1 tag (7)\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		stat_opt_terminate(opt);
		free(opt);
		return 1;
	}

	ci = kopen(argv[optind], &fdci);
	if (ci == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open .cpu file `%s'.\n", __func__, argv[optind]);
		stat_opt_terminate(opt);
		free(opt);
		return 1;
	}
	fpci = gzdopen(fdci, "r");
	gzbuffer(fpci, G_GZIP_BUFFER_SIZE);
	char gzBuffer[G_GZIP_LINE_BUFFER]; gzBuffer[0] = 0;
	gzgets(fpci, gzBuffer, G_GZIP_LINE_BUFFER);
	size_t nCPUlen = strlen(gzBuffer);
	if ('\n'==gzBuffer[nCPUlen-1] || '\r'==gzBuffer[nCPUlen-1]) {
		gzBuffer[nCPUlen-1] = '\0';
	}
	//#define CPU_RECORD_VERSION "0.1a"
	//gzprintf(outfd, "%s\t%d\t%s\n", CPU_RECORD_VERSION, sizeof(cpu_record_t), (opt->cpuflag & CPU_DEBUG_MASK) ? "Text" : "Binary");
	if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] '%s' says '%s'..\n", __func__, argv[optind], gzBuffer);
	int nRecordSize = atoi(gzBuffer);
	if (nRecordSize > sizeof(cpu_record_t)) {
		// larger buffer, it is likely from new version
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Newer version data, %d vs expected %lu\n", __func__, nRecordSize, sizeof(cpu_record_t));
		stat_opt_terminate(opt);
		free(opt);
		return 1;
	} else if (nRecordSize == sizeof(cpu_record_t)) {
		// matching size
	} else {
		// buffer is smaller than current version
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Earlier version data, %d vs expected %lu\n", __func__, nRecordSize, sizeof(cpu_record_t));
		stat_opt_terminate(opt);
		free(opt);
		return 1;
	}
	
	// parse for other information
	kstring_t swVersion; memset(&swVersion, 0, sizeof(kstring_t));
	kstring_t swModule; memset(&swModule, 0, sizeof(kstring_t));
	{
		int *fields, n;
		kstring_t s; memset(&s, 0, sizeof(kstring_t));
		kputs(gzBuffer, &s);
		fields = ksplit(&s, 0, &n);
		if (n>=3) kputs(s.s+fields[1], &swVersion);
		if (n>=4) kputs(s.s+fields[2], &swModule);
		free(s.s); free(fields);
	}
	
	// let's try to determine if the appropriate parameters has been used for the module
	if (swModule.s) {
		if (0==strcmp("CPDL", swModule.s)) {
			opt->nTagFamily = CPU_TAGTYPE_DILINKER;
			if (2!=opt->linkerFlag && bwa_verbose >= 2) fprintf(stderr, "[W::%s] module `%s', using linkerFlag=%d. Shouldn't it be 2?\n", __func__, swModule.s, opt->linkerFlag);
		} else if (0==strcmp("CPSL", swModule.s)) {
			opt->nTagFamily = CPU_TAGTYPE_SILINKER;
			if (1!=opt->linkerFlag && bwa_verbose >= 2) fprintf(stderr, "[W::%s] module `%s', using linkerFlag=%d. Shouldn't it be 1?\n", __func__, swModule.s, opt->linkerFlag);
		} else if (0==strcmp("CPLMP", swModule.s)) {
			opt->nTagFamily = CPU_TAGTYPE_LMP;
			if (2!=opt->linkerFlag && bwa_verbose >= 2) fprintf(stderr, "[W::%s] module `%s', using linkerFlag=%d. Shouldn't it be 2?\n", __func__, swModule.s, opt->linkerFlag);
		} else {
			if (bwa_verbose >= 2) fprintf(stderr, "[W::%s] Unknown module `%s', using linkerFlag=%d.\n", __func__, swModule.s, opt->linkerFlag);
		}
	}

	
	//read off the command parameters
	char szCommand[G_GZIP_LINE_BUFFER]; szCommand[0] = 0;
	{
		gzBuffer[0] = 0;
		gzgets(fpci, gzBuffer, G_GZIP_LINE_BUFFER);
		size_t nCPUlen = strlen(gzBuffer);
		if ('\n'==gzBuffer[nCPUlen-1] || '\r'==gzBuffer[nCPUlen-1]) {
			gzBuffer[nCPUlen-1] = '\0';
		}
		char *szValue = gzBuffer;
		while ('\0' != *szValue && ' ' != *szValue) szValue++; // locate space
		while (' ' == *szValue) szValue++; // skip space
		strcpy(szCommand, szValue);
	}

	// TODO: let's determine the maximum read length and initialize accordingly
	int maxReadLength = MAX_READ_LEN;
	{
		//gzseek can be extremely slow; so we are going to open the file twice!
		//z_off_t startOfRecords = gztell(fpci);
		
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] Determine max read length..\n", __func__);
		
		while ((cpus = stat_readcpu(&n, opt->chunk_size * opt->n_threads, fpci)) != 0) {
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] read %d records..\n", __func__, n);
			
			// might need to pass statistics holder
			int *maxReadLengths = calloc(opt->n_threads, sizeof(int));
			stat_process_maxReadLen(opt, cpus, maxReadLengths, n_processed, n, stat_worker_maxReadLength, &maxReadLength);
			
			// TODO: accummulate the statistics metrics?
			
			free(cpus);
			free(maxReadLengths);
			
			n_processed += n;
			if (bwa_verbose >= 3) {
				t_diff = realtime() - G_t_real;
				fprintf(stderr, "[M::%s] %lld sequences processed, %.0f seq/sec, %.2f min..\n", __func__, n_processed, 1.0*n_processed/t_diff, t_diff/60.0);
			}
		}
		n_processed = 0;
		
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] Max read length = %d.\n", __func__, maxReadLength);
		
		/*
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] Rewinding to start of stream..\n", __func__);
		if (-1==gzseek(fpci, startOfRecords, SEEK_SET)) {
			// ERROR: cannot seek to be start of records
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Fail to rewind to start of records\n", __func__);
			stat_opt_terminate(opt);
			free(opt);
			
			free(swVersion.s); free(swModule.s);
			
			int ret = gzclose(fpci);
			if (Z_OK != ret)
			{
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to close fpci `%s'.\n", __func__, Z_ERRNO == ret ? strerror(errno) : zError(ret));
			} else {
				kclose(ci);
			}
		}
		*/
		
		int ret = gzclose(fpci);
		if (Z_OK != ret)
		{
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to close fpci `%s'.\n", __func__, Z_ERRNO == ret ? strerror(errno) : zError(ret));
		} else {
			kclose(ci);
		}
		
		// let's re-open the file
		ci = kopen(argv[optind], &fdci);
		if (ci == 0) {
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open .cpu file `%s'.\n", __func__, argv[optind]);
			stat_opt_terminate(opt);
			free(opt);
			return 1;
		}
		fpci = gzdopen(fdci, "r");
		gzbuffer(fpci, G_GZIP_BUFFER_SIZE);
		char gzBuffer[G_GZIP_LINE_BUFFER]; gzBuffer[0] = 0;
		gzgets(fpci, gzBuffer, G_GZIP_LINE_BUFFER);
		size_t nCPUlen = strlen(gzBuffer);
		if ('\n'==gzBuffer[nCPUlen-1] || '\r'==gzBuffer[nCPUlen-1]) {
			gzBuffer[nCPUlen-1] = '\0';
		}
		//#define CPU_RECORD_VERSION "0.1a"
		//gzprintf(outfd, "%s\t%d\t%s\n", CPU_RECORD_VERSION, sizeof(cpu_record_t), (opt->cpuflag & CPU_DEBUG_MASK) ? "Text" : "Binary");
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] '%s' says '%s'..\n", __func__, argv[optind], gzBuffer);
		int nRecordSize = atoi(gzBuffer);
		if (nRecordSize > sizeof(cpu_record_t)) {
			// larger buffer, it is likely from new version
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Newer version data, %d vs expected %lu\n", __func__, nRecordSize, sizeof(cpu_record_t));
			stat_opt_terminate(opt);
			free(opt);
			return 1;
		} else if (nRecordSize == sizeof(cpu_record_t)) {
			// matching size
		} else {
			// buffer is smaller than current version
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Earlier version data, %d vs expected %lu\n", __func__, nRecordSize, sizeof(cpu_record_t));
			stat_opt_terminate(opt);
			free(opt);
			return 1;
		}
		
		// parse for other information
		kstring_t swVersion; memset(&swVersion, 0, sizeof(kstring_t));
		kstring_t swModule; memset(&swModule, 0, sizeof(kstring_t));
		{
			int *fields, n;
			kstring_t s; memset(&s, 0, sizeof(kstring_t));
			kputs(gzBuffer, &s);
			fields = ksplit(&s, 0, &n);
			if (n>=3) kputs(s.s+fields[1], &swVersion);
			if (n>=4) kputs(s.s+fields[2], &swModule);
			free(s.s); free(fields);
		}
		
		// let's try to determine if the appropriate parameters has been used for the module
		if (swModule.s) {
			if (0==strcmp("CPDL", swModule.s)) {
				opt->nTagFamily = CPU_TAGTYPE_DILINKER;
				if (2!=opt->linkerFlag && bwa_verbose >= 2) fprintf(stderr, "[W::%s] module `%s', using linkerFlag=%d. Shouldn't it be 2?\n", __func__, swModule.s, opt->linkerFlag);
			} else if (0==strcmp("CPSL", swModule.s)) {
				opt->nTagFamily = CPU_TAGTYPE_SILINKER;
				if (1!=opt->linkerFlag && bwa_verbose >= 2) fprintf(stderr, "[W::%s] module `%s', using linkerFlag=%d. Shouldn't it be 1?\n", __func__, swModule.s, opt->linkerFlag);
			} else if (0==strcmp("CPLMP", swModule.s)) {
				opt->nTagFamily = CPU_TAGTYPE_LMP;
				if (2!=opt->linkerFlag && bwa_verbose >= 2) fprintf(stderr, "[W::%s] module `%s', using linkerFlag=%d. Shouldn't it be 2?\n", __func__, swModule.s, opt->linkerFlag);
			} else {
				if (bwa_verbose >= 2) fprintf(stderr, "[W::%s] Unknown module `%s', using linkerFlag=%d.\n", __func__, swModule.s, opt->linkerFlag);
			}
		}
		
		
		//read off the command parameters
		char szCommand[G_GZIP_LINE_BUFFER]; szCommand[0] = 0;
		{
			gzBuffer[0] = 0;
			gzgets(fpci, gzBuffer, G_GZIP_LINE_BUFFER);
			size_t nCPUlen = strlen(gzBuffer);
			if ('\n'==gzBuffer[nCPUlen-1] || '\r'==gzBuffer[nCPUlen-1]) {
				gzBuffer[nCPUlen-1] = '\0';
			}
			char *szValue = gzBuffer;
			while ('\0' != *szValue && ' ' != *szValue) szValue++; // locate space
			while (' ' == *szValue) szValue++; // skip space
			strcpy(szCommand, szValue);
		}
		// TODO: clean up the above code
	}
	
	void (*stat_worker)(void*,int,int) = (2==opt->linkerFlag) ? stat_worker_dl : stat_worker_sl;
	//stat_classification_t stat_class; memset(&stat_class, 0, sizeof(stat_classification_t));
	stat_classification_t stat_class; stat_classification_init(&stat_class, maxReadLength);
	
	while ((cpus = stat_readcpu(&n, opt->chunk_size * opt->n_threads, fpci)) != 0) {
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d records..\n", __func__, n);
		
		// might need to pass statistics holder
		stat_classification_t *stat_classes = calloc(opt->n_threads, sizeof(stat_classification_t));
		for(i=0; i<opt->n_threads; ++i) stat_classification_init(&(stat_classes[i]), maxReadLength);
		stat_process_seqs(opt, cpus, stat_classes, n_processed, n, maxReadLength, stat_worker, &stat_class);
		
		// TODO: accummulate the statistics metrics?
			
		free(cpus);
		for(i=0; i<opt->n_threads; ++i) stat_classification_terminate(&(stat_classes[i]));
		free(stat_classes);
		
		n_processed += n;
		if (bwa_verbose >= 3) {
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s] %lld sequences processed, %.0f seq/sec, %.2f min..\n", __func__, n_processed, 1.0*n_processed/t_diff, t_diff/60.0);
		}
	}
	
	// OUTPUT: final statistics summary
	int64_t ulPairs = n_processed;
	if (opt->flag&MEM_F_PE) ulPairs = ulPairs >> 1;
	
	fprintf(stdout, "##CPU\t%s\n", swVersion.s ? swVersion.s : CPU_version);
	fprintf(stdout, "##MODULE\t%s\n", swModule.s);
	fprintf(stdout, "##COMMAND\t%s\n", szCommand);
	fprintf(stdout, ">>Library information\n");
	fprintf(stdout, "#Measure\tValue\n");
	fprintf(stdout, "Filename\t%s\n", argv[optind]);
	fprintf(stdout, "Total reads\t%lld\n", n_processed);
	if (opt->flag&MEM_F_PE) fprintf(stdout, "Total pairs\t%lld\n", ulPairs);
	//fprintf(stdout, "Read length\t%d\n", stat_class.readLen);
	fprintf(stdout, "Read length\t%d\n", maxReadLength);
	if (-1==opt->classType) {
		fprintf(stdout, "Statistics filter\tn.a.\n");
	} else {
		fprintf(stdout, "Statistics filter\t%d, %s\n", opt->classType, ((1==opt->linkerFlag) ? outputSLClassToUserStr(opt->classType) : outputClassToUserStr(opt->classType)));
	}
	fprintf(stdout, ">>END\n");
	
	fprintf(stdout, ">>Linker category\n");
	if (opt->flag&MEM_F_PE) fprintf(stdout, "#Category\t#Pair\t%%Pair\t#Tags\n");
	else fprintf(stdout, "#Category\t#Read\t%%Read\t#Tags\n");
	if (1==opt->linkerFlag) {
		int64_t ulReported = 0;
		// we report total number of pairs with linker for easier comparison
		int64_t ulWithLinkerPairs = 0;
		int64_t ulWithLinkerTags = 0;
		for (i=0; i<CPU_OUTPUT_COUNT; ++i) {
			if (i!=CPU_OUTPUT_NONE) {
				ulWithLinkerPairs+=stat_class.classification[i];
				ulWithLinkerTags+=stat_class.classTags[i];
			}
		}
		fprintf(stdout, "%s\t%lld\t%.2f%%\t%lld\n", "Linker detected", ulWithLinkerPairs, ulWithLinkerPairs*100.0f/ulPairs, ulWithLinkerTags);
		for (i=0; i<CPU_OUTPUT_SL_COUNT; ++i) {
			int j=slCatOrders[i];
			fprintf(stdout, "%s\t%d\t%.2f%%\t%d\n", outputSLClassToUserStr(j), stat_class.classification[j], stat_class.classification[j]*100.0f/ulPairs, stat_class.classTags[j]);
			ulReported += stat_class.classification[j];
		}
		if (ulReported != ulPairs) {
			ulReported = ulPairs - ulReported;
			fprintf(stdout, "%s\t%lld\t%.2f%%\tn.a.\n", "Un-accounted For ?!", ulReported, ulReported*100.0f/ulPairs);
		}
	} else {
		int64_t ulReported = 0;
		// we report total number of pairs with linker for easier comparison
		int64_t ulWithLinkerPairs = 0;
		int64_t ulWithLinkerTags = 0;
		for (i=0; i<CPU_OUTPUT_COUNT; ++i) {
			if (i!=CPU_OUTPUT_NONE) {
				ulWithLinkerPairs+=stat_class.classification[i];
				ulWithLinkerTags+=stat_class.classTags[i];
			}
		}
		fprintf(stdout, "%s\t%lld\t%.2f%%\t%lld\n", "Linker detected", ulWithLinkerPairs, ulWithLinkerPairs*100.0f/ulPairs, ulWithLinkerTags);
		for (i=0; i<CPU_OUTPUT_COUNT; ++i) {
			int j=catOrders[i];
			fprintf(stdout, "%s\t%d\t%.2f%%\t%d\n", outputClassToUserStr(j), stat_class.classification[j], stat_class.classification[j]*100.0f/ulPairs, stat_class.classTags[j]);
			ulReported += stat_class.classification[j];
		}
		if (ulReported != ulPairs) {
			ulReported = ulPairs - ulReported;
			fprintf(stdout, "%s\t%lld\t%.2f%%\tn.a.\n", "Un-accounted For ?!", ulReported, ulReported*100.0f/ulPairs);
		}
	}
	fprintf(stdout, ">>END\n");
	
	// end-adapter summary
	fprintf(stdout, ">>End-adapter position\n");
	if (opt->flag&MEM_F_PE) {
		fprintf(stdout, "#BasePosition\t#R/1\t#R/2\t%%R/1\t%%R/2\n");
		for (i=0; i<=maxReadLength; ++i) {
			fprintf(stdout, "%d\t%d\t%d\t%.2f%%\t%.2f%%\n", i, stat_class.endAdapter[0][i], stat_class.endAdapter[1][i], stat_class.endAdapter[0][i]*100.0f/ulPairs, stat_class.endAdapter[1][i]*100.0f/ulPairs);
		}
	} else {
		fprintf(stdout, "#BasePosition\t#R/1\t%%R/1\n");
		for (i=0; i<=maxReadLength; ++i) {
			fprintf(stdout, "%d\t%d\t%.2f%%\n", i, stat_class.endAdapter[0][i], stat_class.endAdapter[0][i]*100.0f/ulPairs);
		}
	}
	fprintf(stdout, ">>END\n");

	// linker summary
	fprintf(stdout, ">>Linker position\n");
	fprintf(stdout, "#BasePosition");
	if (1==opt->linkerFlag) {
		for (i=0; i<REPORT_SLINKER_COUNT; ++i) { fprintf(stdout, "\t#%s-R/1", slinkerLabels[i]); }
		if (opt->flag&MEM_F_PE) for (i=0; i<REPORT_SLINKER_COUNT; ++i) { fprintf(stdout, "\t#%s-R/2", slinkerLabels[i]); }
		for (i=0; i<REPORT_SLINKER_COUNT; ++i) { fprintf(stdout, "\t%%%s-R/1", slinkerLabels[i]); }
		if (opt->flag&MEM_F_PE) for (i=0; i<REPORT_SLINKER_COUNT; ++i) { fprintf(stdout, "\t%%%s-R/2", slinkerLabels[i]); }
	} else {
		for (i=0; i<REPORT_LINKER_COUNT; ++i) { fprintf(stdout, "\t#%s-R/1", linkerLabels[i]); }
		if (opt->flag&MEM_F_PE) for (i=0; i<REPORT_LINKER_COUNT; ++i) { fprintf(stdout, "\t#%s-R/2", linkerLabels[i]); }
		for (i=0; i<REPORT_LINKER_COUNT; ++i) { fprintf(stdout, "\t%%%s-R/1", linkerLabels[i]); }
		if (opt->flag&MEM_F_PE) for (i=0; i<REPORT_LINKER_COUNT; ++i) { fprintf(stdout, "\t%%%s-R/2", linkerLabels[i]); }
	}
	fprintf(stdout, "\n");
	if (1==opt->linkerFlag) {
		for (i=0; i<=maxReadLength; ++i) {
			fprintf(stdout, "%d", i);
			// numbers
			for (j=0; j<REPORT_SLINKER_COUNT; ++j) { fprintf(stdout, "\t%d", stat_class.linker[0][slinkerOrders[j]][i]); }
			if (opt->flag&MEM_F_PE) for (j=0; j<REPORT_SLINKER_COUNT; ++j) { fprintf(stdout, "\t%d", stat_class.linker[1][slinkerOrders[j]][i]); }
			// percentage
			for (j=0; j<REPORT_SLINKER_COUNT; ++j) { fprintf(stdout, "\t%.2f%%", stat_class.linker[0][slinkerOrders[j]][i]*100.0f/ulPairs); }
			if (opt->flag&MEM_F_PE) for (j=0; j<REPORT_SLINKER_COUNT; ++j) { fprintf(stdout, "\t%.2f%%", stat_class.linker[1][slinkerOrders[j]][i]*100.0f/ulPairs); }
			fprintf(stdout, "\n");
		}
	} else {
		for (i=0; i<=maxReadLength; ++i) {
			fprintf(stdout, "%d", i);
			// numbers
			for (j=0; j<REPORT_LINKER_COUNT; ++j) { fprintf(stdout, "\t%d", stat_class.linker[0][linkerOrders[j]][i]); }
			if (opt->flag&MEM_F_PE) for (j=0; j<REPORT_LINKER_COUNT; ++j) { fprintf(stdout, "\t%d", stat_class.linker[1][linkerOrders[j]][i]); }
			// percentage
			for (j=0; j<REPORT_LINKER_COUNT; ++j) { fprintf(stdout, "\t%.2f%%", stat_class.linker[0][linkerOrders[j]][i]*100.0f/ulPairs); }
			if (opt->flag&MEM_F_PE) for (j=0; j<REPORT_LINKER_COUNT; ++j) { fprintf(stdout, "\t%.2f%%", stat_class.linker[1][linkerOrders[j]][i]*100.0f/ulPairs); }
			fprintf(stdout, "\n");
		}
	}
	fprintf(stdout, ">>END\n");

	// taglen summary
	fprintf(stdout, ">>Tag length distribution\n");
	if (opt->flag&MEM_F_PE) {
		fprintf(stdout, "#BasePosition\t#TagLeft-R/1\t#TagRight-R/1\t#TagLeft-R/2\t#TagRight-R/2\t%%TagLeft-R/1\t%%TagRight-R/1\t%%TagLeft-R/2\t%%TagRight-R/2\n");
		for (i=0; i<=maxReadLength; ++i) {
			fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n", i,
					stat_class.taglen[0][0][i], stat_class.taglen[0][1][i], stat_class.taglen[1][0][i], stat_class.taglen[1][1][i],
					stat_class.taglen[0][0][i]*100.0f/ulPairs, stat_class.taglen[0][1][i]*100.0f/ulPairs,
					stat_class.taglen[1][0][i]*100.0f/ulPairs, stat_class.taglen[1][1][i]*100.0f/ulPairs);
		}
	} else {
		fprintf(stdout, "#BasePosition\t#TagLeft-R/1\t#TagRight-R/1\t%%TagLeft-R/1\t%%TagRight-R/1\n");
		for (i=0; i<=maxReadLength; ++i) {
			fprintf(stdout, "%d\t%d\t%d\t%.2f%%\t%.2f%%\n", i,
					stat_class.taglen[0][0][i], stat_class.taglen[0][1][i],
					stat_class.taglen[0][0][i]*100.0f/ulPairs, stat_class.taglen[0][1][i]*100.0f/ulPairs);
		}
	}
	fprintf(stdout, ">>END\n");
	// END - OUTPUT
	
	stat_opt_terminate(opt);
	free(opt);
	
	stat_classification_terminate(&stat_class);
	free(swVersion.s); free(swModule.s);

	//err_gzclose(fpci); kclose(ci);
	int ret = gzclose(fpci);
	if (Z_OK != ret)
	{
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to close fpci `%s'.\n", __func__, Z_ERRNO == ret ? strerror(errno) : zError(ret));
	} else {
		kclose(ci);
	}
	
	return 0;
}


