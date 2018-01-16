#include <stdio.h>
#include <string.h>
#include "utils.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1a-r1"
#endif

#include "kstring.h"

int main_adapter(int argc, char *argv[]);
int main_linker(int argc, char *argv[]);
int main_stag(int argc, char *argv[]);
int main_tag(int argc, char *argv[]);
int main_lmp(int argc, char *argv[]);
int main_read(int argc, char *argv[]);
int main_fq(int argc, char *argv[]);
int main_stat(int argc, char *argv[]);
int main_memaln(int argc, char *argv[]);
int main_pair(int argc, char *argv[]);
int main_span(int argc, char *argv[]);
int main_dedup(int argc, char *argv[]);
int main_mtdedup(int argc, char *argv[]);
int main_cluster(int argc, char *argv[]);

char *bwa_pg;
const char *CPU_version=PACKAGE_VERSION;

double G_t_real;

/* This table is used to transform nucleotide letters into numbers. */
int8_t nt_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: cpu (chIA-pET utilities)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Chee-Hong WONG <cheehongsg@gmail.com>\n\n");
	fprintf(stderr, "Usage:   cpu <command> [options]\n\n");
	fprintf(stderr, "Command: adapter       locate adapter sequence\n");
	fprintf(stderr, "         linker        locate linker sequence\n");
	fprintf(stderr, "         stag          locate single-tag sequence\n");
	fprintf(stderr, "         tag           locate di-tag sequence\n");
	fprintf(stderr, "         lmp           locate LMP sequence\n");
	fprintf(stderr, "         read          time read loop\n");
	fprintf(stderr, "         fq            generate fastq\n");
	fprintf(stderr, "         stat          generate statistics\n");
	fprintf(stderr, "         memaln        perform mapping\n");
	fprintf(stderr, "         pair          conclude pairing for mapped reads\n");
	fprintf(stderr, "         span          compute intra-chromosomal span distribution\n");
	fprintf(stderr, "         dedup         PETs de-duplication\n");
	fprintf(stderr, "         mtdedup       PETs de-duplication (multi-threaded)\n");
	fprintf(stderr, "         cluster       cluster PETs\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Note: To use CPU, <todo...>\n");
	fprintf(stderr, "      <todo...>\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int i, ret;
	kstring_t pg = {0,0,0};
	G_t_real = realtime();
	ksprintf(&pg, "@PG\tID:cpu\tPN:cpu\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "adapter") == 0) ret = main_adapter(argc-1, argv+1);
	else if (strcmp(argv[1], "linker") == 0) ret = main_linker(argc-1, argv+1);
	else if (strcmp(argv[1], "stag") == 0) ret = main_stag(argc-1, argv+1);
	else if (strcmp(argv[1], "tag") == 0) ret = main_tag(argc-1, argv+1);
	else if (strcmp(argv[1], "lmp") == 0) ret = main_lmp(argc-1, argv+1);
	else if (strcmp(argv[1], "read") == 0) ret = main_read(argc-1, argv+1);
	else if (strcmp(argv[1], "fq") == 0) ret = main_fq(argc-1, argv+1);
	else if (strcmp(argv[1], "stat") == 0) ret = main_stat(argc-1, argv+1);
	else if (strcmp(argv[1], "memaln") == 0) ret = main_memaln(argc-1, argv+1);
	else if (strcmp(argv[1], "pair") == 0) ret = main_pair(argc-1, argv+1);
	else if (strcmp(argv[1], "span") == 0) ret = main_span(argc-1, argv+1);
	else if (strcmp(argv[1], "dedup") == 0) ret = main_dedup(argc-1, argv+1);
	else if (strcmp(argv[1], "mtdedup") == 0) ret = main_mtdedup(argc-1, argv+1);
	else if (strcmp(argv[1], "cluster") == 0) ret = main_cluster(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	err_fflush(stdout);
	err_fclose(stdout);
	if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - G_t_real, cputime());
	}
	return ret;
}
