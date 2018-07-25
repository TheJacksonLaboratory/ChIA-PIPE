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
#include "cpufq.h"

extern double G_t_real;

#define G_GZIP_LINE_BUFFER (8*1024)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

// TODO: clean up tags at the end

typedef struct {
	int chunk_size;
	int n_threads;
	int T;
	int flag;               // see MEM_F_* macros
	
	kstring_t outputPrefix;
	
	int cpuflag;               // see MEM_F_* macros
	int minTagLen;
	int nTagFamily;
} fq_opt_t;

fq_opt_t *fq_opt_init()
{
	fq_opt_t *o;
	o = calloc(1, sizeof(fq_opt_t));
	
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->T = 21;
	o->flag = 0;
	
	ks_set(CPU_OUTPUT_PREFIX, &(o->outputPrefix));
	
	o->cpuflag = CPU_FASTQ;
	o->cpuflag &= ~CPU_COMPRESS;
	
	// TODO: minimal tag length for fastq generation!
	o->minTagLen = 18;
	
	o->nTagFamily = CPU_TAGTYPE_DILINKER;
	
	return o;
}

void fq_opt_terminate(fq_opt_t *o)
{
	free(o->outputPrefix.s);
}

typedef struct {
	const fq_opt_t *opt;

	bseq1_t *seqs;
	
	int64_t n_processed;
	
    cpu_record_t *cpus;
	cpu_fq_t *cpu_fqs;
} worker_t;


static void fq_worker_sl(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
	} else {
		processPairedSTag (i + w->n_processed + 1,
						  &(w->cpus[i<<1|0]), &(w->seqs[i<<1|0]), &(w->cpu_fqs[i<<1|0]),
						  &(w->cpus[i<<1|1]), &(w->seqs[i<<1|1]), &(w->cpu_fqs[i<<1|1]),
						  w->opt->minTagLen, w->opt->cpuflag & (CPU_FASTQ|CPU_RAW), w->opt->nTagFamily);
	}
}

static void fq_worker_dl(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
	} else {
		// TODO: output type consideration
		// 1) full fastq based on their classification
		//    GOOD for relooking into the classification
		// 2) fastq based on classification but the read is split into their respective tags
		//    GOOD for performing mapping
		// CPU_OUTPUT_NONE     : just write R/1 follow by R/2
		// CPU_OUTPUT_TIED     : just write R/1 follow by R/2 for now
		//                       in future, to trouble shoot the cases, we need them to be written out
		//                       in to proper pairing information
		// CPU_OUTPUT_CONFLICT : just write R/1 follow by R/2 for now
		//                       for troubleshooting, we need them to be written out with proper pairing information
		// CPU_OUTPUT_HL       : just write R/1 follow by R/2 for now
		//                       for troubleshooting, we need them to be written out with proper pairing information
		//                       to get more value from the data, we need this to be paired (for ohter analysis)
		// CPU_OUTPUT_FL_C     : write out the proper pairing
		//                       for consistency checking, we can check consecutive records
		// CPU_OUTPUT_FL_NC    : write out the proper pairing
		//                     : for consistency checking, we can check consecutive records
		//
		// CONSIDERATIONS:
		// 1. for scaling, it will have been best to have 1 read-pair so that there is not record to record dependency
		// 2. also, pairing is necessary if we wish to recover |----<FL>---<gap>-------| cases
		// 3. so, in the best case, we only have 2 records to consider when looking at the results
		// 4. but more often, we will have multiple records to consider which means simpler ranked partition will not work
		// 5. Thus, it is best to keep all of them separated, i.e. single-end read treatment
		// 6. We can then have another routine to sieve thru' these results
		// 7. @<pair-id>:<read1/2>:<#Tag>:<TagId>:<linkerType>:<tiedLinkertype>:<finalLinkerType>:<classification>
		//

		processPairedTag (i + w->n_processed + 1,
							  &(w->cpus[i<<1|0]), &(w->seqs[i<<1|0]), &(w->cpu_fqs[i<<1|0]),
							  &(w->cpus[i<<1|1]), &(w->seqs[i<<1|1]), &(w->cpu_fqs[i<<1|1]),
							  w->opt->minTagLen, w->opt->cpuflag & (CPU_FASTQ|CPU_RAW|CPU_RC_READ), w->opt->nTagFamily);
	}
}

void fq_process_seqs(const fq_opt_t *opt, bseq1_t *seqs, cpu_record_t *cpus, cpu_fq_t *cpu_fqs, int64_t n_processed, int n, void (*fq_worker)(void*,int,int))
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	w.opt = opt;
	w.seqs = seqs;
	w.cpus = cpus;
	w.cpu_fqs = cpu_fqs;
	//w.n_processed = n_processed;
	w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	kt_for(opt->n_threads, fq_worker, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
}

cpu_record_t *fq_readcpu(int *nread, int n_, gzFile ci1_)
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
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to read %d cpu record(s). %d read instead.\n", __func__, n_, i);
	}
	*nread = i;
	return cpus;
}

int main_fq(int argc, char *argv[])
{
	int i, c, n;
	fq_opt_t *opt;
	gzFile fp = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	cpu_record_t *cpus;
	int fd;
	void *ko = 0;
	int fdci;
	gzFile fpci = 0;
	void *ci = 0;
	int64_t n_processed = 0;
	
	gzFile outfds[CPU_OUTPUT_COUNT] = {0};
	
	double t_diff;

	opt = fq_opt_init();
	while ((c = getopt(argc, argv, "WCpt:O:T:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'T') opt->minTagLen = atoi(optarg), opt->minTagLen = opt->minTagLen >= 0 ? opt->minTagLen : 20;
		else if (c == 'O') { ks_set(optarg, &(opt->outputPrefix)); }
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'C') opt->cpuflag |= CPU_COMPRESS;
		else if (c == 'W') opt->cpuflag |= CPU_RC_READ;
		else {
			fq_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu fq [options] <in.cpu> <in.fq>\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p         first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -O STR     output prefix [%s]\n", opt->outputPrefix.s);
		fprintf(stderr, "       -T INT     min. tag len to be written to .fastq [%d]\n", opt->minTagLen);
		fprintf(stderr, "       -C         compress .fastq\n");
		fprintf(stderr, "       -W         reverse complement the read for .fastq\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		fq_opt_terminate(opt);
		free(opt);
		return 1;
	}

	ci = kopen(argv[optind], &fdci);
	if (ci == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open .cpu file `%s'.\n", __func__, argv[optind]);
		fq_opt_terminate(opt);
		free(opt);
		return 1;
	}
#if 0
	fpci = gzdopen(fdci, "r");
	gzbuffer(fpci, G_GZIP_BUFFER_SIZE);
	//ks = kseq_init(fpci);
	char gzBuffer[G_GZIP_LINE_BUFFER]; gzBuffer[0] = 0;
	gzgets(fpci, gzBuffer, G_GZIP_LINE_BUFFER);
	size_t nCPUlen = strlen(gzBuffer);
	if ('\n'==gzBuffer[nCPUlen-1] || '\r'==gzBuffer[nCPUlen-1]) {
		gzBuffer[nCPUlen-1] = '\0';
	}
	// TODO: check version?!
	//#define CPU_RECORD_VERSION "0.1a"
	//gzprintf(outfd, "%s\t%d\t%s\n", CPU_RECORD_VERSION, sizeof(cpu_record_t), (opt->cpuflag & CPU_DEBUG_MASK) ? "Text" : "Binary");
	if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] .cpu file '%s' says '%s'..\n", __func__, argv[optind], gzBuffer);
	int nRecordSize = atoi(gzBuffer);
	if (nRecordSize > sizeof(cpu_record_t)) {
		// larger buffer, it is likely from new version
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Newer version data, %d vs expected %lu\n", __func__, nRecordSize, sizeof(cpu_record_t));
		fq_opt_terminate(opt);
		free(opt);
		return 1;
	} else if (nRecordSize == sizeof(cpu_record_t)) {
		// matching size
	} else {
		// buffer is smaller than current version
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Earlier version data, %d vs expected %lu\n", __func__, nRecordSize, sizeof(cpu_record_t));
		fq_opt_terminate(opt);
		free(opt);
		return 1;
	}
#else
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
		fq_opt_terminate(opt);
		free(opt);
		return 1;
	} else if (nRecordSize == sizeof(cpu_record_t)) {
		// matching size
	} else {
		// buffer is smaller than current version
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] Earlier version data, %d vs expected %lu\n", __func__, nRecordSize, sizeof(cpu_record_t));
		fq_opt_terminate(opt);
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
			if (bwa_verbose >= 3) fprintf(stderr, "[M::%s] module `%s', using TagFamily=%d.\n", __func__, swModule.s, opt->nTagFamily);
		} else if (0==strcmp("CPSL", swModule.s)) {
			opt->nTagFamily = CPU_TAGTYPE_SILINKER;
			if (bwa_verbose >= 3) fprintf(stderr, "[M::%s] module `%s', using TagFamily=%d.\n", __func__, swModule.s, opt->nTagFamily);
		} else if (0==strcmp("CPLMP", swModule.s)) {
			opt->nTagFamily = CPU_TAGTYPE_LMP;
			if (bwa_verbose >= 3) fprintf(stderr, "[M::%s] module `%s', using TagFamily=%d.\n", __func__, swModule.s, opt->nTagFamily);
		} else {
			if (bwa_verbose >= 2) fprintf(stderr, "[W::%s] Unknown module `%s', using TagFamily=%d.\n", __func__, swModule.s, opt->nTagFamily);
		}
	} else {
		if (bwa_verbose >= 2) fprintf(stderr, "[W::%s] No module, default to TagFamily=%d.\n", __func__, opt->nTagFamily);
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
	
	void (*fq_worker)(void*,int,int) = (CPU_TAGTYPE_DILINKER==opt->nTagFamily || CPU_TAGTYPE_LMP==opt->nTagFamily) ? fq_worker_dl : fq_worker_sl;
#endif

	ko = kopen(argv[optind+1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind+1]);
		fq_opt_terminate(opt);
		free(opt);
		return 1;
	}
	fp = gzdopen(fd, "r");
	gzbuffer(fp, G_GZIP_BUFFER_SIZE);
	ks = kseq_init(fp);
	
	// TODO: should initialize shared resources across threads
	if (init_CPTags_Ouputs(opt->cpuflag, opt->outputPrefix.s, outfds)) {
		terminate_CPTags_Ouputs(outfds);
		
		fq_opt_terminate(opt);
		free(opt);
		kseq_destroy(ks);
		err_gzclose(fp); kclose(ko);
		err_gzclose(fpci); kclose(ci);
		
		return 1;
	}
	
	while ((seqs = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
		int64_t size = 0;
		if ((opt->flag & MEM_F_PE) && (n&1) == 1) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] odd number of reads in the PE mode; last read dropped\n", __func__);
			n = n>>1<<1;
		}
		for (i = 0; i < n; ++i) size += seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp)..\n", __func__, n, (long)size);
		
		// might need to pass statistics holder
		int nRead = 0;
		cpus = fq_readcpu(&nRead, n, fpci);
		if (nRead != n) {
			for (i = 0; i < n; ++i) {
				free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual);
			}
			free(seqs);
			break;
		}
		
		cpu_fq_t *cpu_fqs = calloc(n>0?n:1, sizeof(cpu_fq_t));
		fq_process_seqs(opt, seqs, cpus, cpu_fqs, n_processed, n, fq_worker);
		
		// OUTPUT:
		for (i = 0; i < n; ++i) {
			if (cpu_fqs[i].fq.l>=opt->minTagLen) gzwrite(outfds[cpus[i].classification], cpu_fqs[i].fq.s, cpu_fqs[i].fq.l);
		}
		// END - OUTPUT
		
		for (i = 0; i < n; ++i) {
			free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual);
			free(cpu_fqs[i].fq.s);
		}
		free(seqs);
		
		free(cpus);
		free(cpu_fqs);
		
		n_processed += n;
		if (bwa_verbose >= 3) {
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s] %lld sequences processed, %.0f seq/sec, %.2f min..\n", __func__, n_processed, 1.0*n_processed/t_diff, t_diff/60.0);
		}
	}
	
	terminate_CPTags_Ouputs(outfds);
	
	fq_opt_terminate(opt);
	free(opt);
	kseq_destroy(ks);

	//err_gzclose(fp); kclose(ko);
	int ret = gzclose(fp);
	if (Z_OK != ret)
	{
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to close fp `%s'.\n", __func__, Z_ERRNO == ret ? strerror(errno) : zError(ret));
	} else {
		kclose(ko);
	}
	//err_gzclose(fpci); kclose(ci);
	ret = gzclose(fpci);
	if (Z_OK != ret)
	{
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to close fpci `%s'.\n", __func__, Z_ERRNO == ret ? strerror(errno) : zError(ret));
	} else {
		kclose(ci);
	}
	
	return 0;
}


