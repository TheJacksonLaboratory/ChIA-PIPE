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

extern double G_t_real;

// available in >=v1.2.4
//#define USE_LARGER_GZIP_BUFFER
#define G_GZIP_BUFFER_SIZE (128*1024)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

typedef struct {
	int chunk_size;
	int n_threads;
	int T;
	int flag;               // see MEM_F_* macros
} read_opt_t;

read_opt_t *read_opt_init()
{
	read_opt_t *o;
	o = calloc(1, sizeof(read_opt_t));
	
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->T = 21;
	o->flag = 0;
	return o;
}

typedef struct {
	const read_opt_t *opt;
	bseq1_t *seqs;
	
	int64_t n_processed;
} worker_t;

static void read_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;

	if (!(w->opt->flag&MEM_F_PE)) {
	} else {
	}
}

void read_process_seqs(const read_opt_t *opt, bseq1_t *seqs, int64_t n_processed, int n)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	w.opt = opt;
	w.seqs = seqs;
	//w.n_processed = n_processed;
	w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	kt_for(opt->n_threads, read_worker, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
}

int main_read(int argc, char *argv[])
{
	read_opt_t *opt;
	int fd, fd2, i, c, n;
	int copy_comment = 0;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	void *ko = 0, *ko2 = 0;
	int64_t n_processed = 0;
	
	double t_diff;

	opt = read_opt_init();
	while ((c = getopt(argc, argv, "p5:3:t:T:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else {
			free(opt);
			return 1;
		}
	}
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu read [options] <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p         first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		free(opt);
		return 1;
	}

	ko = kopen(argv[optind], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind]);
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
		read_process_seqs(opt, seqs, n_processed, n);
		
		// OUTPUT:
		for (i = 0; i < n; i+=2) {
			fprintf(stdout, "%s", seqs[i].name);
			fprintf(stdout, "\n");
		}
		// END - OUTPUT
		
		for (i = 0; i < n; ++i) {
			//err_fputs(seqs[i].sam, stdout);
			free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual);
			//free(seqs[i].sam);
		}
		free(seqs);
		
		n_processed += n;
		if (bwa_verbose >= 3) {
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s] %lld sequences processed, %.0f seq/sec, %.2f min..\n", __func__, n_processed, 1.0*n_processed/t_diff, t_diff/60.0);
		}
	}
	
	free(opt);
	kseq_destroy(ks);
	err_gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	
	return 0;
}


