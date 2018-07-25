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

#include "cpulinker.h"
#include "cputag.h"

extern double G_t_real;

// available in >=v1.2.4
//#define USE_LARGER_GZIP_BUFFER
#define G_GZIP_BUFFER_SIZE (128*1024)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

typedef struct {
	int32_t linker_type;
	s_align *pla;
} linker_align_t;

typedef struct {
	int32_t linker_type;
	int32_t linker_bufsize;
	int8_t *linker_num;
} linker_t;

typedef struct {
	kstring_t halfLinkerA;
	kstring_t halfLinkerB;
	kstring_t rc_halfLinkerA;
	kstring_t rc_halfLinkerB;
	kstring_t fullLinkerAA;
	kstring_t fullLinkerBB;
	kstring_t fullLinkerAB;
	kstring_t fullLinkerBA;
	
	int8_t  aln_flag;
	int32_t match;
	int32_t mismatch;
	int32_t gap_open;
	int32_t gap_extension;
	
	int32_t nNT;
	int8_t* mata;
	
	int32_t path;
	int32_t filter;

	linker_t align_hlA;
	linker_t align_hlB;
	linker_t align_hlrcA;
	linker_t align_hlrcB;
	linker_t align_flAA;
	linker_t align_flBB;
	linker_t align_flAB;
	linker_t align_flBA;
	
	int32_t filterHL;
	int32_t filterFL;
	
	int chunk_size;
	int n_threads;
	int T;
	int flag;               // see MEM_F_* macros
	int cpuflag;               // see MEM_F_* macros
} linker_opt_t;

void initLinker(linker_opt_t *o)
{
	ks_set(o->halfLinkerA.s, &(o->rc_halfLinkerA));
	ks_set(o->halfLinkerB.s, &(o->rc_halfLinkerB));
	reverse_comple(o->halfLinkerA.s, o->rc_halfLinkerA.s, o->halfLinkerA.l);
	reverse_comple(o->halfLinkerB.s, o->rc_halfLinkerB.s, o->halfLinkerB.l);
	
	// generate full linker
	kstring_t tAA; memset(&tAA, 0, sizeof(kstring_t));
	ksprintf(&tAA, "%s%s", o->halfLinkerA.s, o->rc_halfLinkerA.s); ks_set(tAA.s, &(o->fullLinkerAA)); free(tAA.s);
	kstring_t tBB; memset(&tBB, 0, sizeof(kstring_t));
	ksprintf(&tBB, "%s%s", o->halfLinkerB.s, o->rc_halfLinkerB.s); ks_set(tBB.s, &(o->fullLinkerBB)); free(tBB.s);
	kstring_t tAB; memset(&tAB, 0, sizeof(kstring_t));
	ksprintf(&tAB, "%s%s", o->halfLinkerA.s, o->rc_halfLinkerB.s); ks_set(tAB.s, &(o->fullLinkerAB)); free(tAB.s);
	kstring_t tBA; memset(&tBA, 0, sizeof(kstring_t));
	ksprintf(&tBA, "%s%s", o->halfLinkerB.s, o->rc_halfLinkerA.s); ks_set(tBA.s, &(o->fullLinkerBA)); free(tBA.s);
}

void init_LinkerScoreMatrix (int8_t *mat, int32_t match, int32_t mismatch)
{
	int32_t l, m, k;
	// initialize scoring matrix
	for (l = k = 0; LIKELY(l < 4); ++l) {
		for (m = 0; LIKELY(m < 4); ++m) mat[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
		mat[k++] = 0; // ambiguous base
	}
	for (m = 0; LIKELY(m < 5); ++m) mat[k++] = 0;
}

void init_ssw_Linker(linker_t *lssw, const kstring_t *lnt, int32_t linkerType)
{
	lssw->linker_type = linkerType;
	lssw->linker_bufsize = (int32_t) lnt->l;
	lssw->linker_num = (int8_t*)realloc(lssw->linker_num, lssw->linker_bufsize);
	int32_t m;
	for (m = 0; m < lssw->linker_bufsize; ++m) lssw->linker_num[m] = nt_table[(int)lnt->s[m]];
}

void init_AlignmentLinker(linker_opt_t *o)
{
	initLinker(o);
	
	init_ssw_Linker(&(o->align_hlA), &(o->halfLinkerA), LINKER_A);
	init_ssw_Linker(&(o->align_hlB), &(o->halfLinkerB), LINKER_B);
	init_ssw_Linker(&(o->align_hlrcA), &(o->rc_halfLinkerA), LINKER_RC_A);
	init_ssw_Linker(&(o->align_hlrcB), &(o->rc_halfLinkerB), LINKER_RC_B);
	
	init_ssw_Linker(&(o->align_flAA), &(o->fullLinkerAA), LINKER_AA);
	init_ssw_Linker(&(o->align_flBB), &(o->fullLinkerBB), LINKER_BB);
	init_ssw_Linker(&(o->align_flAB), &(o->fullLinkerAB), LINKER_AB);
	init_ssw_Linker(&(o->align_flBA), &(o->fullLinkerBA), LINKER_BA);
}

void initLinkerAlignmentParameters(linker_opt_t *o)
{
	// re-set up the scoring matrix (user override from command line)
	init_LinkerScoreMatrix(o->mata, o->match, o->mismatch);
	
	init_AlignmentLinker(o);
}

linker_opt_t *linker_opt_init()
{
	linker_opt_t *o;
	o = calloc(1, sizeof(linker_opt_t));
	
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->T = 21;
	o->flag = 0;
	o->cpuflag = 0;
	memset(&(o->halfLinkerA), 0, sizeof(kstring_t));
	memset(&(o->halfLinkerB), 0, sizeof(kstring_t));
	memset(&(o->rc_halfLinkerA), 0, sizeof(kstring_t));
	memset(&(o->rc_halfLinkerB), 0, sizeof(kstring_t));
	
	ks_set(G_CP_HalfLinkerA, &(o->halfLinkerA));
	ks_set(G_CP_HalfLinkerB, &(o->halfLinkerB));

	o->match = 1;
	o->mismatch = 1;
	o->gap_open = 6;
	o->gap_extension = 1;
	
	o->nNT = 4+1; // A, C, G, T, N
	o->mata = (int8_t*)calloc(o->nNT*o->nNT, sizeof(int8_t));
	init_LinkerScoreMatrix(o->mata, o->match, o->mismatch);
	
	o->aln_flag = 1;
	o->path = 0; // TODO
	o->filter = 0; // TODO
	
	memset(&(o->align_hlA), 0, sizeof(linker_t));
	memset(&(o->align_hlB), 0, sizeof(linker_t));
	memset(&(o->align_hlrcA), 0, sizeof(linker_t));
	memset(&(o->align_hlrcB), 0, sizeof(linker_t));

	memset(&(o->align_flAA), 0, sizeof(linker_t));
	memset(&(o->align_flBB), 0, sizeof(linker_t));
	memset(&(o->align_flAB), 0, sizeof(linker_t));
	memset(&(o->align_flBA), 0, sizeof(linker_t));
	
	o->filterHL = 15; // TODO: to determine based on the scroing matrix
	o->filterFL = 26; // TODO: to determine based on the scroing matrix
	
	initLinker(o);
	
	return o;
}

void linker_opt_terminate(linker_opt_t *o)
{
	free(o->halfLinkerA.s);
	free(o->halfLinkerB.s);
	free(o->rc_halfLinkerA.s);
	free(o->rc_halfLinkerB.s);
	
	free(o->fullLinkerAA.s);
	free(o->fullLinkerBB.s);
	free(o->fullLinkerAB.s);
	free(o->fullLinkerBA.s);

	free(o->mata);
	
	free(o->align_hlA.linker_num);
	free(o->align_hlB.linker_num);
	free(o->align_hlrcA.linker_num);
	free(o->align_hlrcB.linker_num);
	
	free(o->align_flAA.linker_num);
	free(o->align_flBB.linker_num);
	free(o->align_flAB.linker_num);
	free(o->align_flBA.linker_num);
	
}

typedef struct {
	const linker_opt_t *opt;
	
	bseq1_t *seqs;

	linker_align_t *las;
	
	int64_t n_processed;
} worker_t;

static void linker_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		//mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a);
		//mem_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
		//free(w->regs[i].a);
	} else {
		// TODO: performance, move all the one time initialization out
		// WCH: simplest implementation will require O(4xFL) or O(8xHL)
		//      there are several implementation possible to reduce the number of alignment needed
		//      1) fatest: O(2 5/7 xHL) but will require stitching the alignment information
		//      2) faster: O(3 4/7 xHL) but will not require stitching the alignmet information
		
		// Process Read/1
		linker_align_t *pLARec = &(w->las[i<<1|0]);
		pLARec->linker_type = LINKER_NONE; pLARec->pla = 0;
		
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
		s_align *reshlA1 = 0; const linker_t *hlA1 = &(w->opt->align_hlA);
		reshlA1 = ssw_align (p1, hlA1->linker_num, hlA1->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
		if (reshlA1 && reshlA1->score1>=w->opt->filterHL) {
			// can be AA / AB
			s_align *resflAA1 = 0; const linker_t *flAA1 = &(w->opt->align_flAA);
			resflAA1 = ssw_align (p1, flAA1->linker_num, flAA1->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
			if (resflAA1 && resflAA1->score1>=w->opt->filterFL) {
				// CONCLUSION: linkerAA
				// TODO: what if we only have the left partial, such that we cannot tell AA vs AB?
				pLARec->linker_type = LINKER_AA; pLARec->pla = resflAA1;
			} else {
				s_align *resflAB1 = 0; const linker_t *flAB1 = &(w->opt->align_flAB);
				resflAB1 = ssw_align (p1, flAB1->linker_num, flAB1->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
				if (resflAB1 && resflAB1->score1>=w->opt->filterFL) {
					// CONCLUSION: linkerAB
					// TODO: what if we only have the left partial, such that we cannot tell AB vs AA?
					pLARec->linker_type = LINKER_AB; pLARec->pla = resflAB1;
				} else {
					// CONCLUSION: left partial of AA, i.e. half linker A
					pLARec->linker_type = LINKER_A; pLARec->pla = reshlA1;
				}
				if (pLARec->pla != resflAB1) align_destroy(resflAB1);
			}
			if (pLARec->pla != resflAA1) align_destroy(resflAA1);
		} else {
			// TODO: do we need to reperform 'ssw_init()' ?
			
			s_align *reshlB1 = 0; const linker_t *hlB1 = &(w->opt->align_hlB);
			reshlB1 = ssw_align (p1, hlB1->linker_num, hlB1->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
			if (reshlB1 && reshlB1->score1>=w->opt->filterHL) {
				// can be BB / BA
				s_align *resflBB1 = 0; const linker_t *flBB1 = &(w->opt->align_flBB);
				resflBB1 = ssw_align (p1, flBB1->linker_num, flBB1->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
				if (resflBB1 && resflBB1->score1>=w->opt->filterFL) {
					// CONCLUSION: linkerBB
					// TODO: what if we only have the left partial, such that we cannot tell BB vs BA?
					pLARec->linker_type = LINKER_BB; pLARec->pla = resflBB1;
				} else {
					s_align *resflBA1 = 0; const linker_t *flBA1 = &(w->opt->align_flBA);
					resflBA1 = ssw_align (p1, flBA1->linker_num, flBA1->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
					if (resflBA1 && resflBA1->score1>=w->opt->filterFL) {
						// CONCLUSION: linkerBA
						// TODO: what if we only have the left partial, such that we cannot tell BA vs BB?
						pLARec->linker_type = LINKER_BA; pLARec->pla = resflBA1;
					} else {
						// CONCLUSION: left partial of BB, i.e. half linker B
						pLARec->linker_type = LINKER_B; pLARec->pla = reshlB1;
					}
					if (pLARec->pla != resflBA1) align_destroy(resflBA1);
				}
				if (pLARec->pla != resflBB1) align_destroy(resflBB1);
			} else {
				// not a full linker or partial left of full linker
				s_align *reshlA1rc = 0; const linker_t *hlA1rc = &(w->opt->align_hlrcA);
				reshlA1rc = ssw_align (p1, hlA1rc->linker_num, hlA1rc->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
				if (reshlA1rc && reshlA1rc->score1>=w->opt->filterHL) {
					// right most of the AA, i.e. right partial of AA
					pLARec->linker_type = LINKER_RC_A; pLARec->pla = reshlA1rc;
				} else {
					s_align *reshlB1rc = 0; const linker_t *hlB1rc = &(w->opt->align_hlrcB);
					reshlB1rc = ssw_align (p1, hlB1rc->linker_num, hlB1rc->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
					if (reshlB1rc && reshlB1rc->score1>=w->opt->filterHL) {
						// right most of the BB, i.e. right partial of BB
						pLARec->linker_type = LINKER_RC_B; pLARec->pla = reshlB1rc;
					} else {
						// no linker of any form
					}
					if (pLARec->pla != reshlB1rc) align_destroy(reshlB1rc);
				}
				if (pLARec->pla != reshlA1rc) align_destroy(reshlA1rc);
			}
			if (pLARec->pla != reshlB1) align_destroy(reshlB1);
		}
		if (pLARec->pla != reshlA1) align_destroy(reshlA1);
		
		// release profile for Read/1
		init_destroy(p1);
		free(num);
		
		// Process Read/2
		pLARec = &(w->las[i<<1|1]);
		pLARec->linker_type = LINKER_NONE; pLARec->pla = 0;
		
		s_profile *p2 = 0;
		readLen = w->seqs[i<<1|1].l_seq;
		maskLen = readLen / 2;
		//int32_t maskLen = w->opt->halfLinkerA.l / 2;
		num = (int8_t*)malloc(readLen+1);
		for (m = 0; m < readLen; ++m) num[m] = nt_table[(int)w->seqs[i<<1|1].seq[m]];
		p2 = ssw_init(num, readLen, w->opt->mata, w->opt->nNT, 2);
		
		// TODO: this should be saved and processed centrally by the caller
		// test linkers alignment
		s_align *reshlA2 = 0; const linker_t *hlA2 = &(w->opt->align_hlA);
		reshlA2 = ssw_align (p2, hlA2->linker_num, hlA2->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
		if (reshlA2 && reshlA2->score1>=w->opt->filterHL) {
			// can be AA / AB
			s_align *resflAA2 = 0; const linker_t *flAA2 = &(w->opt->align_flAA);
			resflAA2 = ssw_align (p2, flAA2->linker_num, flAA2->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
			if (resflAA2 && resflAA2->score1>=w->opt->filterFL) {
				// CONCLUSION: linkerAA
				// TODO: what if we only have the left partial, such that we cannot tell AA vs AB?
				pLARec->linker_type = LINKER_AA; pLARec->pla = resflAA2;
			} else {
				s_align *resflAB2 = 0; const linker_t *flAB2 = &(w->opt->align_flAB);
				resflAB2 = ssw_align (p2, flAB2->linker_num, flAB2->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
				if (resflAB2 && resflAB2->score1>=w->opt->filterFL) {
					// CONCLUSION: linkerAB
					// TODO: what if we only have the left partial, such that we cannot tell AB vs AA?
					pLARec->linker_type = LINKER_AB; pLARec->pla = resflAB2;
				} else {
					// CONCLUSION: left partial of AA, i.e. half linker A
					pLARec->linker_type = LINKER_A; pLARec->pla = reshlA2;
				}
				if (pLARec->pla != resflAB2) align_destroy(resflAB2);
			}
			if (pLARec->pla != resflAA2) align_destroy(resflAA2);
		} else {
			// TODO: do we need to reperform 'ssw_init()' ?
			
			s_align *reshlB2 = 0; const linker_t *hlB2 = &(w->opt->align_hlB);
			reshlB2 = ssw_align (p2, hlB2->linker_num, hlB2->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
			if (reshlB2 && reshlB2->score1>=w->opt->filterHL) {
				// can be BB / BA
				s_align *resflBB2 = 0; const linker_t *flBB2 = &(w->opt->align_flBB);
				resflBB2 = ssw_align (p2, flBB2->linker_num, flBB2->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
				if (resflBB2 && resflBB2->score1>=w->opt->filterFL) {
					// CONCLUSION: linkerBB
					// TODO: what if we only have the left partial, such that we cannot tell BB vs BA?
					pLARec->linker_type = LINKER_BB; pLARec->pla = resflBB2;
				} else {
					s_align *resflBA2 = 0; const linker_t *flBA2 = &(w->opt->align_flBA);
					resflBA2 = ssw_align (p2, flBA2->linker_num, flBA2->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
					if (resflBA2 && resflBA2->score1>=w->opt->filterFL) {
						// CONCLUSION: linkerBA
						// TODO: what if we only have the left partial, such that we cannot tell BA vs BB?
						pLARec->linker_type = LINKER_BA; pLARec->pla = resflBA2;
					} else {
						// CONCLUSION: left partial of BB, i.e. half linker B
						pLARec->linker_type = LINKER_B; pLARec->pla = reshlB2;
					}
					if (pLARec->pla != resflBA2) align_destroy(resflBA2);
				}
				if (pLARec->pla != resflBB2) align_destroy(resflBB2);
			} else {
				// not a full linker or partial left of full linker
				s_align *reshlA2rc = 0; const linker_t *hlA2rc = &(w->opt->align_hlrcA);
				reshlA2rc = ssw_align (p2, hlA2rc->linker_num, hlA2rc->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
				if (reshlA2rc && reshlA2rc->score1>=w->opt->filterHL) {
					// right most of the AA, i.e. right partial of AA
					pLARec->linker_type = LINKER_RC_A; pLARec->pla = reshlA2rc;
				} else {
					s_align *reshlB2rc = 0; const linker_t *hlB2rc = &(w->opt->align_hlrcB);
					reshlB2rc = ssw_align (p2, hlB2rc->linker_num, hlB2rc->linker_bufsize, w->opt->gap_open, w->opt->gap_extension, w->opt->aln_flag, w->opt->filter, 0, maskLen);
					if (reshlB2rc && reshlB2rc->score1>=w->opt->filterHL) {
						// right most of the BB, i.e. right partial of BB
						pLARec->linker_type = LINKER_RC_B; pLARec->pla = reshlB2rc;
					} else {
						// no linker of any form
					}
					if (pLARec->pla != reshlB2rc) align_destroy(reshlB2rc);
				}
				if (pLARec->pla != reshlA2rc) align_destroy(reshlA2rc);
			}
			if (pLARec->pla != reshlB2) align_destroy(reshlB2);
		}
		if (pLARec->pla != reshlA2) align_destroy(reshlA2);
		
		// release profile for Read/2
		init_destroy(p2);
		free(num);
	}
}

void linker_process_seqs(const linker_opt_t *opt, bseq1_t *seqs, linker_align_t *las, int64_t n_processed, int n)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	w.opt = opt;
	w.seqs = seqs;
	w.las = las;
	//w.n_processed = n_processed;
	w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	kt_for(opt->n_threads, linker_worker, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
}

int main_linker(int argc, char *argv[])
{
	linker_opt_t *opt;
	int fd, fd2, i, c, n;
	int copy_comment = 0;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	void *ko = 0, *ko2 = 0;
	int64_t n_processed = 0;

	double t_diff;
	
	opt = linker_opt_init();
	while ((c = getopt(argc, argv, "SpA:B:t:T:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'T') opt->T = atoi(optarg), opt->T = opt->T > 0 ? opt->T : 21;
		else if (c == 'A') { ks_set(optarg, &(opt->halfLinkerA)); }
		else if (c == 'B') { ks_set(optarg, &(opt->halfLinkerB)); }
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'S') opt->cpuflag |= CPU_SEQ;
		else {
			linker_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu linker [options] <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -A STR     half linker A [%s]\n", opt->halfLinkerA.s);
		fprintf(stderr, "       -B STR     half linker B [%s]\n", opt->halfLinkerB.s);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p         first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT     minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		linker_opt_terminate(opt);
		free(opt);
		return 1;
	}

	ko = kopen(argv[optind], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind]);
		
		linker_opt_terminate(opt);
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
				
				linker_opt_terminate(opt);
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
	
	initLinker(opt); // TODO: is this necessary?
	initLinkerAlignmentParameters(opt);
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
		linker_align_t *las = calloc(n, sizeof(linker_align_t));
		linker_process_seqs(opt, seqs, las, n_processed, n);
		
		// OUTPUT:
		for (i = 0; i < n; i+=2) {
			// TODO: this should be saved and processed centrally by the caller
			// rcAdapt3<start,end>,read/1<start,end>,<score>
			// rcAdapt5<start,end>,read/2<start,end>,<score>
			// fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", result1->ref_begin1, result1->ref_end1, result1->read_begin1, result1->read_end1, result1->score1, result2->ref_begin1, result2->ref_end1, result2->read_begin1, result2->read_end1, result2->score1);
			
			fprintf(stdout, "%s", seqs[i].name);
			
			linker_align_t *pLARec1 = &(las[i]);
			int32_t r1linker = pLARec1->linker_type;
			if (LINKER_NONE==pLARec1->linker_type) {
				fprintf(stdout, "\t-1\t-1\t-1\t-1\t-1\t.");
			} else {
				fprintf(stdout, "\t%d\t%d\t%d\t%d\t%d\t%s", pLARec1->pla->ref_begin1, pLARec1->pla->ref_end1, pLARec1->pla->read_begin1, pLARec1->pla->read_end1, pLARec1->pla->score1, linkerTypeToStr(pLARec1->linker_type));
				align_destroy(pLARec1->pla);
			}
			linker_align_t *pLARec2 = &(las[i+1]);
			int32_t r2linker = pLARec2->linker_type;
			if (LINKER_NONE==pLARec2->linker_type) {
				fprintf(stdout, "\t-1\t-1\t-1\t-1\t-1\t.");
			} else {
				fprintf(stdout, "\t%d\t%d\t%d\t%d\t%d\t%s", pLARec2->pla->ref_begin1, pLARec2->pla->ref_end1, pLARec2->pla->read_begin1, pLARec2->pla->read_end1, pLARec2->pla->score1, linkerTypeToStr(pLARec2->linker_type));
				align_destroy(pLARec2->pla);
			}
			if (r1linker == r2linker) {
				fprintf(stdout, "\t=");
			} else {
				if (LINKER_AB==r1linker && LINKER_BA==r2linker) {
					fprintf(stdout, "\t=");
				} else if (LINKER_BA==r1linker && LINKER_AB==r2linker) {
					fprintf(stdout, "\t=");
				} else {
					fprintf(stdout, "\t?");
				}
			}
			
			// TODO: use a switch to decide when to output the sequence
			fprintf(stdout, "\t%s\t%s", seqs[i].seq, seqs[i+1].seq);
			
			fprintf(stdout, "\n");
		}
		// END - OUTPUT

		free(las);
		
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
	
	linker_opt_terminate(opt);
	free(opt);
	kseq_destroy(ks);
	err_gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	
	return 0;
}


