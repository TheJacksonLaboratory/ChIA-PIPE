#include <zlib.h>
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

#include "ssw.h"

#include "cpuadapter.h"
#include "cpulinker.h"
#include "cputag.h"

// |G_CP_LMP_HalfLinkerA| = 19
#define G_CP_LMP_HalfLinkerA "CTGTCTCTTATACACATCT"

extern double G_t_real;

// available in >=v1.2.4
#define G_GZIP_BUFFER_SIZE (128*1024)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

typedef struct {
	int16_t adapter_type;
	s_align *paa;
	int16_t linker_type;
	s_align *pla;
	int16_t tied_linker_type;
} tag_align_t;

typedef struct {
	int16_t tag_type;
	int16_t tag_bufsize;
	int8_t *tag_num;
} tag_t;

typedef struct {
	kstring_t adapter5;
	kstring_t adapter3;
	kstring_t rc_adapter5;
	kstring_t rc_adapter3;
	
	kstring_t halfLinkerA;
	kstring_t rc_halfLinkerA;
	kstring_t fullLinkerAA;
	
	int8_t  linker_aln_flag;
	int32_t linker_match;
	int32_t linker_mismatch;
	int32_t linker_gap_open;
	int32_t linker_gap_extension;
	
	int32_t linker_path;
	int32_t linker_filter;

	int8_t  adapter_aln_flag;
	int32_t adapter_match;
	int32_t adapter_mismatch;
	int32_t adapter_gap_open;
	int32_t adapter_gap_extension;
	
	int32_t adapter_path;
	int32_t adapter_filter;
	
	int32_t nNT;
	int8_t* adapter_mata;
	int8_t* linker_mata;
	
	tag_t align_hlA;
	tag_t align_hlrcA;
	tag_t align_flAA;
	
	int32_t filterHL;
	int32_t filterFL;
#ifdef POSITIVE_SCORE_MATRIX
#else
	int32_t minHLOverlap;
	int32_t filterHLOverlap;
	int32_t minFLOverlap;
	int32_t filterFLOverlap;
#endif
	
	tag_t align_adapter5;
	tag_t align_adapter3;
	tag_t align_rcAdapter5;
	tag_t align_rcAdapter3;
	
	int32_t filterAdapter5;
	int32_t filterAdapter3;
#ifdef POSITIVE_SCORE_MATRIX
#else
	int32_t minAD5Overlap;
	int32_t filterAD5Overlap;
	int32_t minAD3Overlap;
	int32_t filterAD3Overlap;
#endif
	
	//int aas_n; // starting amount of reads alignment that we allocated for
	//int las_n; // starting amount of reads alignment that we allocated for
	int cpus_n; // starting amount of reads alignment that we allocated for
	
	int chunk_size;
	int n_threads;
	int T;
	int flag;               // see MEM_F_* macros
	int cpuflag;               // see MEM_F_* macros
	
	int minTagLen;
	int cpu_decisions_n;
	kstring_t outputPrefix;
} lmp_opt_t;

void initLmp(lmp_opt_t *o)
{
	ks_set(o->adapter5.s, &(o->rc_adapter5));
	ks_set(o->adapter3.s, &(o->rc_adapter3));
	reverse_comple(o->adapter5.s, o->rc_adapter5.s, o->adapter5.l);
	reverse_comple(o->adapter3.s, o->rc_adapter3.s, o->adapter3.l);
	
	ks_set(o->halfLinkerA.s, &(o->rc_halfLinkerA));
	reverse_comple(o->halfLinkerA.s, o->rc_halfLinkerA.s, o->halfLinkerA.l);
	
	// generate full linker
	kstring_t tAA; memset(&tAA, 0, sizeof(kstring_t));
	ksprintf(&tAA, "%s%s", o->halfLinkerA.s, o->rc_halfLinkerA.s); ks_set(tAA.s, &(o->fullLinkerAA)); free(tAA.s);
}

void init_LmpScoreMatrix (int8_t *mat, int32_t match, int32_t mismatch)
{
	int32_t l, m, k;
	// initialize scoring matrix
	for (l = k = 0; LIKELY(l < 4); ++l) {
		for (m = 0; LIKELY(m < 4); ++m) mat[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
		mat[k++] = 0; // ambiguous base
	}
	for (m = 0; LIKELY(m < 5); ++m) mat[k++] = 0;
}

void init_ssw_Lmp(tag_t *tssw, const kstring_t *tnt, int32_t tagType)
{
	tssw->tag_type = tagType;
	tssw->tag_bufsize = (int32_t) tnt->l;
	tssw->tag_num = (int8_t*)realloc(tssw->tag_num, tssw->tag_bufsize);
	int32_t m;
	for (m = 0; m < tssw->tag_bufsize; ++m) tssw->tag_num[m] = nt_table[(int)tnt->s[m]];
}

void init_AlignmentLmp(lmp_opt_t *o)
{
	initLmp(o);
	
	init_ssw_Lmp(&(o->align_adapter5), &(o->adapter5), ADAPTER_5P);
	init_ssw_Lmp(&(o->align_adapter3), &(o->adapter3), ADAPTER_3P);
	init_ssw_Lmp(&(o->align_rcAdapter5), &(o->rc_adapter5), ADAPTER_5Prc);
	init_ssw_Lmp(&(o->align_rcAdapter3), &(o->rc_adapter3), ADAPTER_3Prc);
	
	init_ssw_Lmp(&(o->align_hlA), &(o->halfLinkerA), LINKER_A);
	init_ssw_Lmp(&(o->align_hlrcA), &(o->rc_halfLinkerA), LINKER_RC_A);
	
	init_ssw_Lmp(&(o->align_flAA), &(o->fullLinkerAA), LINKER_AA);
}

void initLmpAlignmentParameters(lmp_opt_t *o)
{
	// re-set up the scoring matrix (user override from command line)
	init_LmpScoreMatrix(o->adapter_mata, o->adapter_match, o->adapter_mismatch);
	init_LmpScoreMatrix(o->linker_mata, o->linker_match, o->linker_mismatch);
	
	init_AlignmentLmp(o);
}

lmp_opt_t *lmp_opt_init()
{
	lmp_opt_t *o;
	o = calloc(1, sizeof(lmp_opt_t));
	
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
	
	ksprintf(&(o->adapter5), "%s", G_NEXTERA_5PAdapter);
	ksprintf(&(o->adapter3), "%s", G_NEXTERA_3PAdapter);
	
	memset(&(o->halfLinkerA), 0, sizeof(kstring_t));
	memset(&(o->rc_halfLinkerA), 0, sizeof(kstring_t));
	
	ksprintf(&(o->halfLinkerA), "%s", G_CP_LMP_HalfLinkerA);
	
#ifdef POSITIVE_SCORE_MATRIX
	o->adapter_match = 1;
	o->adapter_mismatch = 1;
	o->adapter_gap_open = 6;
	o->adapter_gap_extension = 1;
#else
	o->adapter_match = 1;
	o->adapter_mismatch = 1;
	o->adapter_gap_open = 1;
	o->adapter_gap_extension = 1;
#endif
	
#ifdef POSITIVE_SCORE_MATRIX
	o->linker_match = 1;
	o->linker_mismatch = 1;
	o->linker_gap_open = 6;
	o->linker_gap_extension = 1;
#else
	o->linker_match = 1;
	o->linker_mismatch = 1;
	o->linker_gap_open = 1;
	o->linker_gap_extension = 1;
#endif
	
	o->nNT = 4+1; // A, C, G, T, N
	o->adapter_mata = (int8_t*)calloc(o->nNT*o->nNT, sizeof(int8_t));
	init_LmpScoreMatrix(o->adapter_mata, o->adapter_match, o->adapter_mismatch);
	o->linker_mata = (int8_t*)calloc(o->nNT*o->nNT, sizeof(int8_t));
	init_LmpScoreMatrix(o->linker_mata, o->linker_match, o->linker_mismatch);
	
	o->adapter_aln_flag = 8;
	o->adapter_path = 0; // TODO
	o->adapter_filter = 0; // TODO
	
	o->linker_aln_flag = 8;
	o->linker_path = 0; // TODO
	o->linker_filter = 0; // TODO
	
	memset(&(o->align_adapter5), 0, sizeof(tag_t));
	memset(&(o->align_adapter3), 0, sizeof(tag_t));
	memset(&(o->align_rcAdapter5), 0, sizeof(tag_t));
	memset(&(o->align_rcAdapter3), 0, sizeof(tag_t));
	
	memset(&(o->align_hlA), 0, sizeof(tag_t));
	memset(&(o->align_hlrcA), 0, sizeof(tag_t));

	memset(&(o->align_flAA), 0, sizeof(tag_t));
	
#ifdef POSITIVE_SCORE_MATRIX
	o->filterAdapter5 = 21; // TODO: to determine based on the scroing matrix
	o->filterAdapter3 = 21; // TODO: to determine based on the scroing matrix
#else
	o->filterAdapter5 = 27; // TODO: to determine based on the scroing matrix
	o->filterAdapter3 = 28; // TODO: to determine based on the scroing matrix
#if 1
	//o->minAD5Overlap = 12; // NOTE: the 5'-end of rc(AD5) is actually the single junction adapter (19bp)
	o->minAD5Overlap = 23;   // NOTE: we require 5 additional bases match from 19bp w/o mismatches
	// NOTE: thus, 23(0x, Score 23), 24(1x, Score=22), 25(2x, Score=21), 26(2x,Score=22), 27(2x,Score=23), 28(3x, Score=22), ...
	o->filterAD5Overlap = (o->minAD5Overlap) * o->adapter_match - 2*o->adapter_mismatch; // TODO: need to re-calculate
	o->minAD3Overlap = 12;
	o->filterAD3Overlap = (o->minAD3Overlap) * o->adapter_match - o->adapter_mismatch; // TODO: need to re-calculate
#else
	o->minAD5Overlap = 22;
	o->filterAD5Overlap = (o->minAD5Overlap) * o->adapter_match - o->adapter_mismatch; // TODO: need to re-calculate
	o->minAD3Overlap = 22;
	o->filterAD3Overlap = (o->minAD3Overlap) * o->adapter_match - o->adapter_mismatch; // TODO: need to re-calculate
#endif
#endif

	// TODO: for LMP, we do not keep just 30 or 4 mismatches for 38 bases
	//       will have to explore the data to tell really!!!
#ifdef POSITIVE_SCORE_MATRIX
	o->filterHL = 15; // TODO: to determine based on the scroing matrix
	o->filterFL = 26; // TODO: to determine based on the scroing matrix
#else
	o->filterHL = 14; // TODO: to determine based on the scroing matrix
	o->filterFL = 30; // TODO: to determine based on the scroing matrix
	o->minHLOverlap = 12;
	o->filterHLOverlap = o->minHLOverlap * o->linker_match - o->linker_mismatch;
	o->minFLOverlap = 32;
	o->filterFLOverlap = (o->minFLOverlap - 3) * o->linker_match - 3 * o->linker_mismatch;
#endif
	
	
	o->cpus_n = 1000; // TODO: we should allocate more for testing
	
	o->cpu_decisions_n = 1000;
	ks_set(CPU_OUTPUT_PREFIX, &(o->outputPrefix));
	o->minTagLen = 20;
	
	initLmp(o);
	
	return o;
}

void lmp_opt_terminate(lmp_opt_t *o)
{
	free(o->adapter5.s);
	free(o->adapter3.s);
	free(o->rc_adapter5.s);
	free(o->rc_adapter3.s);
	
	free(o->adapter_mata);
	
	free(o->align_adapter5.tag_num);
	free(o->align_adapter3.tag_num);
	free(o->align_rcAdapter5.tag_num);
	free(o->align_rcAdapter3.tag_num);

	
	free(o->halfLinkerA.s);
	free(o->rc_halfLinkerA.s);
	
	free(o->fullLinkerAA.s);

	free(o->linker_mata);
	
	free(o->align_hlA.tag_num);
	free(o->align_hlrcA.tag_num);
	
	free(o->align_flAA.tag_num);
	
	free(o->outputPrefix.s);
}

typedef struct {
	const lmp_opt_t *opt;
	
	bseq1_t *seqs;
	
	int64_t n_processed;
	
	cpu_record_t *cpus;
	cpu_fq_t *cpu_fqs;
} worker_t;

static inline int nexteraRCAdapter3Overlap3prime(const s_align *pa, const lmp_opt_t *o, const bseq1_t *seqs)
{
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

static inline int nexteraRCAdapter5Overlap3prime(const s_align *pa, const lmp_opt_t *o, const bseq1_t *seqs)
{
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
#if 1
		// tailored for LMP
	} else if (refSpan>=25) {
		if (refSpan-pa->score1<=4) {
			return 1;
		}
	} else if (refSpan==24) {
		if (refSpan-pa->score1<=2) {
			return 1;
		}
	} else if (refSpan==23) {
		if (refSpan-pa->score1==0) {
			return 1;
		}
#else
	} else if (refSpan>=18) {
		if (refSpan-pa->score1<=4) {
			return 1;
		}
	} else if (refSpan>=12) {
		if (refSpan-pa->score1<=2) {
			return 1;
		}
#endif
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

static inline int flOverlap5prime(const s_align *pa, const lmp_opt_t *o, const bseq1_t *seqs)
{
	if (pa->score1<o->filterFLOverlap) return 0; // did not reach the minimum score
	
	if (0!=pa->read_begin1) return 0; // did not start from 5'-end of read
	
	if (o->fullLinkerAA.l-pa->ref_end1>1) return 0; // the adapter's 3' end was not the end of match!
	
	// TODO: need an array for the score thresholding!!! based on 10% error rate!
	int refSpan = pa->ref_end1 - pa->ref_begin1 + 1;
	if (refSpan>=36) {
		if (refSpan-pa->score1<=8) {
			return 1;
		}
	} else if (refSpan>=32) {
		if (refSpan-pa->score1<=6) {
			return 1;
		}
	} else {
		return 0;
	}
	return 0;
}

static inline int flOverlap3prime(const s_align *pa, const lmp_opt_t *o, const bseq1_t *seqs)
{
	if (pa->score1<o->filterFLOverlap) return 0; // did not reach the minimum score
	
	if (0!=pa->ref_begin1) return 0; // the adapter's 5' end was not the start of match!
	
	// TODO: we cannot test end of read as we have not quality trimmed yet
	// TODO: so, for now, we used quality-checking as valid sign that the read is ending
	if (pa->read_end1!=seqs->l_seq-1) return 0; // did not end at 3'-end of read
	
	// TODO: need an array for the score thresholding!!! based on 10% error rate!
	int refSpan = pa->ref_end1 - pa->ref_begin1 + 1;
	if (refSpan>=36) {
		if (refSpan-pa->score1<=8) {
			return 1;
		}
	} else if (refSpan>=32) {
		if (refSpan-pa->score1<=6) {
			return 1;
		}
	} else {
		return 0;
	}
	return 0;
}

static inline int hlOverlap5prime(const s_align *pa, const lmp_opt_t *o, const bseq1_t *seqs)
{
	if (pa->score1<o->filterHLOverlap) return 0; // did not reach the minimum score
	
	if (0!=pa->read_begin1) return 0; // did not start from 5'-end of read
	
	if (o->halfLinkerA.l-pa->ref_end1>1) return 0; // the adapter's 3' end was not the end of match!
	
	// TODO: need an array for the score thresholding!!! based on 10% error rate!
	int refSpan = pa->ref_end1 - pa->ref_begin1 + 1;
	if (refSpan>=17) {
		if (refSpan-pa->score1<=4) {
			return 1;
		}
	} else if (refSpan>=14) {
		if (refSpan-pa->score1<=2) {
			return 1;
		}
	} else if (refSpan==13) {
		if (refSpan-pa->score1==2) {
			return 1;
		}
	} else {
		return 0;
	}
	return 0;
}

static inline int hlOverlap3prime(const s_align *pa, const lmp_opt_t *o, const bseq1_t *seqs)
{
	if (pa->score1<o->filterHLOverlap) return 0; // did not reach the minimum score
	
	if (0!=pa->ref_begin1) return 0; // the adapter's 5' end was not the start of match!
	
	// TODO: we cannot test end of read as we have not quality trimmed yet
	// TODO: so, for now, we used quality-checking as valid sign that the read is ending
	if (pa->read_end1!=seqs->l_seq-1) return 0; // did not end at 3'-end of read
	
	// TODO: need an array for the score thresholding!!! based on 10% error rate!
	int refSpan = pa->ref_end1 - pa->ref_begin1 + 1;
	if (refSpan>=17) {
		if (refSpan-pa->score1<=4) {
			return 1;
		}
	} else if (refSpan>=14) {
		if (refSpan-pa->score1<=2) {
			return 1;
		}
	} else if (refSpan==13) {
		if (refSpan-pa->score1==2) {
			return 1;
		}
	} else {
		return 0;
	}
	return 0;
}

static void lmp_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		//mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a);
		//mem_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
		//free(w->regs[i].a);
	} else {
		// TODO: performance, move all the one time initialization out
		tag_align_t aRec;
		// WCH: simplest implementation will require O(4xFL) or O(8xHL)
		//      there are several implementation possible to reduce the number of alignment needed
		//      1) fatest: O(2 5/7 xHL) but will require stitching the alignment information
		//      2) faster: O(3 4/7 xHL) but will not require stitching the alignmet information
		
		// Process Read/1
		memset(&aRec, 0, sizeof(tag_align_t));
		//aRec.adapter_type = ADAPTER_NONE; aRec.paa = 0;
		//aRec.linker_type = LINKER_NONE; aRec.pla = 0; aRec.tied_linker_type = LINKER_NONE;
		
		s_profile *p1 = 0;
		int32_t readLen = w->seqs[i<<1|0].l_seq;
		int32_t maskLen = readLen / 2;
		//int32_t maskLen = w->opt->halfLinkerA.l / 2;
		int8_t *num = (int8_t*)malloc(readLen+1);
		int32_t m;
		for (m = 0; m < readLen; ++m) num[m] = nt_table[(int)w->seqs[i<<1|0].seq[m]];
		
		/* Adapter */
		p1 = ssw_init(num, readLen, w->opt->linker_mata, w->opt->nNT, 2);
		
		// TODO: this should be saved and processed centrally by the caller
		// test linkers alignment
		s_align *resNTrc3p = 0; const tag_t *adapterNTrc3p = &(w->opt->align_rcAdapter3);
		resNTrc3p = ssw_align (p1, adapterNTrc3p->tag_num, adapterNTrc3p->tag_bufsize, w->opt->adapter_gap_open, w->opt->adapter_gap_extension, w->opt->adapter_aln_flag, w->opt->adapter_filter, 0, maskLen);
		aRec.paa = resNTrc3p;
		if (resNTrc3p->score1>=w->opt->filterAdapter3) {
			aRec.adapter_type = ADAPTER_3Prc;
			if (0!=nexteraRCAdapter3Overlap3prime(aRec.paa, w->opt, &w->seqs[i<<1|0])) {
				// okie, we have a prefix suffix case
				aRec.adapter_type = ADAPTER_3Prc;
			}
		}
		init_destroy(p1);
		/* END - Adapter */
		
		/* Linker */
		// TODO: might have to use a half width to cater for possible boundary
		p1 = ssw_init(num, (ADAPTER_NONE==aRec.adapter_type) ? readLen : aRec.paa->read_begin1, w->opt->linker_mata, w->opt->nNT, 2);
		
		s_align *resflAA1 = 0; const tag_t *flAA1 = &(w->opt->align_flAA);
		resflAA1 = ssw_align (p1, flAA1->tag_num, flAA1->tag_bufsize, w->opt->linker_gap_open, w->opt->linker_gap_extension, w->opt->linker_aln_flag, w->opt->linker_filter, 0, maskLen);
		aRec.linker_type = LINKER_AA; aRec.pla = resflAA1;
		
		if (aRec.pla->score1>=w->opt->filterFL) {
			// okie, we have a winner for full linker
		} else {
			// full linker is NOT a good solution
			
			if (0!=flOverlap5prime(aRec.pla, w->opt, &w->seqs[i<<1|0])) {
				// okie, we have a prefix fix case
			} else if (0!=flOverlap3prime(aRec.pla, w->opt, &w->seqs[i<<1|0])) {
				// okie, we have a suffix fix case
			} else {
				
				align_destroy(aRec.pla);
				aRec.linker_type = LINKER_NONE; aRec.pla = 0; aRec.tied_linker_type = LINKER_NONE;
				
				// test the half linker!
				s_align *reshlA1 = 0; const tag_t *hlA1 = &(w->opt->align_hlA);
				reshlA1 = ssw_align (p1, hlA1->tag_num, hlA1->tag_bufsize, w->opt->linker_gap_open, w->opt->linker_gap_extension, w->opt->linker_aln_flag, w->opt->linker_filter, 0, maskLen);
				
				s_align *reshlA1rc = 0; const tag_t *hlA1rc = &(w->opt->align_hlrcA);
				reshlA1rc = ssw_align (p1, hlA1rc->tag_num, hlA1rc->tag_bufsize, w->opt->linker_gap_open, w->opt->linker_gap_extension, w->opt->linker_aln_flag, w->opt->linker_filter, 0, maskLen);
				
				if (reshlA1->score1 > reshlA1rc->score1) {
					aRec.linker_type = LINKER_A; aRec.pla = reshlA1;
					align_destroy(reshlA1rc);
				} else if (reshlA1->score1 < reshlA1rc->score1) {
					aRec.linker_type = LINKER_RC_A; aRec.pla = reshlA1rc;
					align_destroy(reshlA1);
				} else {
					// a tie, does it matter?
					aRec.linker_type = LINKER_A; aRec.pla = reshlA1; aRec.tied_linker_type = LINKER_RC_A;
					align_destroy(reshlA1rc);
				}
				
				if (aRec.pla->score1>=w->opt->filterHL) {
					// okie, we have a winner for half linker
				} else {
					if (0!=hlOverlap5prime(aRec.pla, w->opt, &w->seqs[i<<1|0])) {
						// okie, we have a prefix fix case
					} else if (0!=hlOverlap3prime(aRec.pla, w->opt, &w->seqs[i<<1|0])) {
						// okie, we have a suffix fix case
					} else {
						// half linker is NOT a good solution
						//align_destroy(aRec.pla);
						aRec.linker_type = LINKER_NONE;
						//aRec.pla = 0;
						aRec.tied_linker_type = LINKER_NONE;
					}
				}
			}
		}
		
		// release profile for Read/1
		init_destroy(p1);
		free(num);
		/* END - Linker */
		
		/* prepare the record for binary compression writing */
		cpu_record_t *rec = &(w->cpus[i<<1|0]);
		// adapter 1
		rec->adapter_type = aRec.adapter_type;
		s_align *pa = aRec.paa;
		rec->adapter_score1 = pa->score1;
		rec->adapter_ref_begin1 = pa->ref_begin1;
		rec->adapter_ref_end1 = pa->ref_end1;
		rec->adapter_read_begin1 = pa->read_begin1;
		rec->adapter_read_end1 = pa->read_end1;
		// linker 1
		rec->linker_type = aRec.linker_type;
		pa = aRec.pla;
		rec->linker_score1 = pa->score1;
		rec->linker_ref_begin1 = pa->ref_begin1;
		rec->linker_ref_end1 = pa->ref_end1;
		rec->linker_read_begin1 = pa->read_begin1;
		rec->linker_read_end1 = pa->read_end1;
		rec->tied_linker_type = aRec.tied_linker_type;
		/* END - output preparation */
		align_destroy(aRec.paa);
		align_destroy(aRec.pla);
		
		
		// Process Read/2
		memset(&aRec, 0, sizeof(tag_align_t));
		//aRec.adapter_type = ADAPTER_NONE; aRec.paa = 0;
		//aRec.linker_type = LINKER_NONE; aRec.pla = 0; aRec.tied_linker_type = LINKER_NONE;
		
		s_profile *p2 = 0;
		readLen = w->seqs[i<<1|1].l_seq;
		maskLen = readLen / 2;
		//int32_t maskLen = w->opt->halfLinkerA.l / 2;
		num = (int8_t*)malloc(readLen+1);
		for (m = 0; m < readLen; ++m) num[m] = nt_table[(int)w->seqs[i<<1|1].seq[m]];
		
		/* Adapter */
		p2 = ssw_init(num, readLen, w->opt->adapter_mata, w->opt->nNT, 2);
		
		// TODO: this should be saved and processed centrally by the caller
		// test linkers alignment
		s_align *resNTrc5p = 0; const tag_t *adapterNTrc5p = &(w->opt->align_rcAdapter5);
		resNTrc5p = ssw_align (p2, adapterNTrc5p->tag_num, adapterNTrc5p->tag_bufsize, w->opt->adapter_gap_open, w->opt->adapter_gap_extension, w->opt->adapter_aln_flag, w->opt->adapter_filter, 0, maskLen);
		aRec.paa = resNTrc5p;
		if (resNTrc5p->score1>=w->opt->filterAdapter5) {
			aRec.adapter_type = ADAPTER_5Prc;
			if (0!=nexteraRCAdapter5Overlap3prime(aRec.paa, w->opt, &w->seqs[i<<1|0])) {
				// okie, we have a prefix suffix case
				aRec.adapter_type = ADAPTER_5Prc;
			}
		}
		init_destroy(p2);
		/* END - Adapter */
		
		/* Linker */
		// TODO: might have to use a half width to cater for possible boundary
		p2 = ssw_init(num, (ADAPTER_NONE==aRec.adapter_type) ? readLen : aRec.paa->read_begin1, w->opt->linker_mata, w->opt->nNT, 2);
		
		
		s_align *resflAA2 = 0; const tag_t *flAA2 = &(w->opt->align_flAA);
		resflAA2 = ssw_align (p2, flAA2->tag_num, flAA2->tag_bufsize, w->opt->linker_gap_open, w->opt->linker_gap_extension, w->opt->linker_aln_flag, w->opt->linker_filter, 0, maskLen);
		aRec.linker_type = LINKER_AA; aRec.pla = resflAA2;
		
		if (aRec.pla->score1>=w->opt->filterFL) {
			// okie, we have a winner for full linker
		} else {
			// full linker is NOT a good solution
			
			if (0!=flOverlap5prime(aRec.pla, w->opt, &w->seqs[i<<1|1])) {
				// okie, we have a prefix fix case
			} else if (0!=flOverlap3prime(aRec.pla, w->opt, &w->seqs[i<<1|1])) {
				// okie, we have a suffix fix case
			} else {
				align_destroy(aRec.pla);
				aRec.linker_type = LINKER_NONE; aRec.pla = 0; aRec.tied_linker_type = LINKER_NONE;
				
				// test the half linker!
				s_align *reshlA2 = 0; const tag_t *hlA2 = &(w->opt->align_hlA);
				reshlA2 = ssw_align (p2, hlA2->tag_num, hlA2->tag_bufsize, w->opt->linker_gap_open, w->opt->linker_gap_extension, w->opt->linker_aln_flag, w->opt->linker_filter, 0, maskLen);
				
				s_align *reshlA2rc = 0; const tag_t *hlA2rc = &(w->opt->align_hlrcA);
				reshlA2rc = ssw_align (p2, hlA2rc->tag_num, hlA2rc->tag_bufsize, w->opt->linker_gap_open, w->opt->linker_gap_extension, w->opt->linker_aln_flag, w->opt->linker_filter, 0, maskLen);
				
				if (reshlA2->score1 > reshlA2rc->score1) {
					aRec.linker_type = LINKER_A; aRec.pla = reshlA2;
					align_destroy(reshlA2rc);
				} else if (reshlA2->score1 < reshlA2rc->score1) {
					aRec.linker_type = LINKER_RC_A; aRec.pla = reshlA2rc;
					align_destroy(reshlA2);
				} else {
					// a tie, does it matter?
					aRec.linker_type = LINKER_A; aRec.pla = reshlA2; aRec.tied_linker_type = LINKER_RC_A;
					align_destroy(reshlA2rc);
				}
				
				if (aRec.pla->score1>=w->opt->filterHL) {
					// okie, we have a winner for half linker
				} else {
					if (0!=hlOverlap5prime(aRec.pla, w->opt, &w->seqs[i<<1|1])) {
						// okie, we have a prefix fix case
					} else if (0!=hlOverlap3prime(aRec.pla, w->opt, &w->seqs[i<<1|1])) {
						// okie, we have a suffix fix case
					} else {
						// half linker is NOT a good solution
						//align_destroy(aRec.pla);
						aRec.linker_type = LINKER_NONE;
						//aRec.pla = 0;
						aRec.tied_linker_type = LINKER_NONE;
					}
				}
			}
		}
		
		// release profile for Read/2
		init_destroy(p2);
		free(num);
		/* END - Linker */
		
		
		/* prepare the record for binary compression writing */
		rec = &(w->cpus[i<<1|1]);
		// adapter 1
		rec->adapter_type = aRec.adapter_type;
		pa = aRec.paa;
		rec->adapter_score1 = pa->score1;
		rec->adapter_ref_begin1 = pa->ref_begin1;
		rec->adapter_ref_end1 = pa->ref_end1;
		rec->adapter_read_begin1 = pa->read_begin1;
		rec->adapter_read_end1 = pa->read_end1;
		// linker 1
		rec->linker_type = aRec.linker_type;
		pa = aRec.pla;
		rec->linker_score1 = pa->score1;
		rec->linker_ref_begin1 = pa->ref_begin1;
		rec->linker_ref_end1 = pa->ref_end1;
		rec->linker_read_begin1 = pa->read_begin1;
		rec->linker_read_end1 = pa->read_end1;
		rec->tied_linker_type = aRec.tied_linker_type;
		/* END - output preparation */
		
		align_destroy(aRec.paa);
		align_destroy(aRec.pla);
		
		// perform the tags
		processPairedTag (i + w->n_processed + 1,
						  &(w->cpus[i<<1|0]), &(w->seqs[i<<1|0]), &(w->cpu_fqs[i<<1|0]),
						  &(w->cpus[i<<1|1]), &(w->seqs[i<<1|1]), &(w->cpu_fqs[i<<1|1]),
						  w->opt->minTagLen, w->opt->cpuflag & (CPU_FASTQ|CPU_RAW|CPU_RC_READ), CPU_TAGTYPE_LMP);
	}
}

void lmp_process_seqs(const lmp_opt_t *opt, bseq1_t *seqs, cpu_record_t *cpus, cpu_fq_t *cpu_fqs, int64_t n_processed, int n)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	w.opt = opt;
	w.seqs = seqs;
	w.cpus = cpus;
	w.cpu_fqs = cpu_fqs;
	//w.n_processed = n_processed;
	w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	kt_for(opt->n_threads, lmp_worker, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
}

int main_lmp(int argc, char *argv[])
{
	lmp_opt_t *opt;
	int fd, fd2, i, c, n;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	void *ko = 0, *ko2 = 0;
	int64_t n_processed = 0;

	gzFile outfd;
	gzFile outfds[CPU_OUTPUT_COUNT] = {0};
	
	double t_diff;

	opt = lmp_opt_init();
	while ((c = getopt(argc, argv, "RFWLNCSQp5:3:A:B:t:T:O:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'n') opt->T = atoi(optarg), opt->T = opt->T > 0 ? opt->T : 21;
		else if (c == 'A') { ks_set(optarg, &(opt->halfLinkerA)); }
		else if (c == '5') { ks_set(optarg, &(opt->adapter5)); }
		else if (c == '3') { ks_set(optarg, &(opt->adapter3)); }
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'L') opt->cpuflag |= CPU_LABEL;
		else if (c == 'N') opt->cpuflag |= CPU_NAME;
		else if (c == 'S') opt->cpuflag |= CPU_SEQ;
		else if (c == 'Q') opt->cpuflag |= CPU_QUALITY;
		else if (c == 'C') opt->cpuflag |= CPU_COMPRESS;
		else if (c == 'W') opt->cpuflag |= CPU_RC_READ;
		else if (c == 'F') opt->cpuflag &= ~CPU_FASTQ;
		else if (c == 'R') opt->cpuflag |= CPU_RAW;
		else if (c == 'T') opt->minTagLen = atoi(optarg), opt->minTagLen = opt->minTagLen >= 0 ? opt->minTagLen : 20;
		else if (c == 'O') { ks_set(optarg, &(opt->outputPrefix)); }
		else {
			lmp_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu tag [options] <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -A STR     half linker A [%s]\n", opt->halfLinkerA.s);
		fprintf(stderr, "       -5 STR     5' adapter [%s]\n", opt->adapter5.s);
		fprintf(stderr, "       -3 STR     3' adapter [%s]\n", opt->adapter3.s);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p         first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -O STR     output prefix [%s]\n", opt->outputPrefix.s);
		fprintf(stderr, "       -n INT     minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -T INT     min. tag len to be written to .fastq\n");
		fprintf(stderr, "       -C         compress .fastq\n");
		fprintf(stderr, "       -W         reverse complement the read for .fastq\n");
		fprintf(stderr, "       -F         skip .fastq writing\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -L         write label instead of numeric encoding\n");
		fprintf(stderr, "       -N         write read name\n");
		fprintf(stderr, "       -S         write sequence\n");
		fprintf(stderr, "       -Q         write quality\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		lmp_opt_terminate(opt);
		free(opt);
		return 1;
	}

	ko = kopen(argv[optind], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind]);
		
		lmp_opt_terminate(opt);
		free(opt);
		return 1;
	}
	fp = gzdopen(fd, "r");
	gzbuffer(fp, G_GZIP_BUFFER_SIZE);
	ks = kseq_init(fp);
	if (optind + 1 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file will be ignored.\n", __func__);
		} else {
			ko2 = kopen(argv[optind + 1], &fd2);
			if (ko2 == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
				
				lmp_opt_terminate(opt);
				free(opt);
				return 1;
			}
			fp2 = gzdopen(fd2, "r");
			gzbuffer(fp2, G_GZIP_BUFFER_SIZE);
			ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}
	
	initLmp(opt); // TODO: is this necessary?
	initLmpAlignmentParameters(opt);
	// TODO: should initialize shared resources across threads
	
	if (opt->cpuflag & CPU_FASTQ) {
		if (init_CPTags_Ouputs(opt->cpuflag, opt->outputPrefix.s, outfds)) {
			terminate_CPTags_Ouputs(outfds);
			
			lmp_opt_terminate(opt);
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
	
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.cpu", opt->outputPrefix.s);
		//const char *mode = (opt->cpuflag & CPU_COMPRESS) ? "w6" : "wT"; //"w0";
		const char *mode = (0==(opt->cpuflag & CPU_DEBUG_MASK)) ? "w6" : "wT";
		outfd = gzopen (filename.s, mode);
		if (outfd == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing: %s\n", __func__, filename.s, strerror (errno));
			lmp_opt_terminate(opt);
			free(opt);
			return 1;
		}
		free(filename.s);
		gzbuffer(outfd, G_GZIP_BUFFER_SIZE);
	}
	
	gzprintf(outfd, "%d\t%s\tCPLMP\t%s\n", sizeof(cpu_record_t), CPU_RECORD_VERSION, (opt->cpuflag & CPU_DEBUG_MASK) ? "Text" : "Binary");
	{
		gzprintf(outfd, "COMMAND:");
		for (i = 0; i < argc; ++i)
			gzprintf(outfd, " %s", argv[i]);
		gzprintf(outfd, "\n");
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
		cpu_record_t *cpus = calloc(n, sizeof(cpu_record_t));
		cpu_fq_t *cpu_fqs = calloc(n, sizeof(cpu_fq_t));
		lmp_process_seqs(opt, seqs, cpus, cpu_fqs, n_processed, n);
		
		
		// OUTPUT:
		if (0==(opt->cpuflag & CPU_DEBUG_MASK)) {
			
			int nWritten = gzwrite(outfd, cpus, n*sizeof(cpu_record_t));
			nWritten /= sizeof(cpu_record_t);
			if (nWritten != n) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to write %d cpu records. %d written instead.\n", __func__, n, nWritten);
			}
		} else {
			for (i = 0; i < n; i+=2) {
				// rcAdapt3<start,end>,read/1<start,end>,<score>,<best>
				// rcAdapt5<start,end>,read/2<start,end>,<score>,<best>
				// linker1<start,end>,read/1<start,end>,<score>,<best>,<tied>
				// linker2<start,end>,read/2<start,end>,<score>,<best>,<tied>
				// read name
				// <read1>,<read2>
				// <quality1>,<quality2>
				
				/* adapter */
				cpu_record_t *pCPU1 = &(cpus[i]);
				cpu_record_t *pCPU2 = &(cpus[i+1]);
				
				if (opt->cpuflag & CPU_LABEL) {
					gzprintf(outfd, "%d\t%d\t%d\t%d\t%d", pCPU1->adapter_ref_begin1, pCPU1->adapter_ref_end1, pCPU1->adapter_read_begin1, pCPU1->adapter_read_end1, pCPU1->adapter_score1);
					gzprintf(outfd, "\t%s", adapterTypeToStr(pCPU1->adapter_type));
					
					gzprintf(outfd, "\t%d\t%d\t%d\t%d\t%d", pCPU2->adapter_ref_begin1, pCPU2->adapter_ref_end1, pCPU2->adapter_read_begin1, pCPU2->adapter_read_end1, pCPU2->adapter_score1);
					gzprintf(outfd, "\t%s", adapterTypeToStr(pCPU2->adapter_type));
					/* END - adapter */
					
					/* linker */
					gzprintf(outfd, "\t%d\t%d\t%d\t%d\t%d", pCPU1->linker_ref_begin1, pCPU1->linker_ref_end1, pCPU1->linker_read_begin1, pCPU1->linker_read_end1, pCPU1->linker_score1);
					gzprintf(outfd, "\t%s\t%s", linkerTypeToStr(pCPU1->linker_type), linkerTypeToStr(pCPU1->tied_linker_type));
					gzprintf(outfd, "\t%d\t%d\t%d\t%d\t%d", pCPU2->linker_ref_begin1, pCPU2->linker_ref_end1, pCPU2->linker_read_begin1, pCPU2->linker_read_end1, pCPU2->linker_score1);
					gzprintf(outfd, "\t%s\t%s", linkerTypeToStr(pCPU2->linker_type), linkerTypeToStr(pCPU2->tied_linker_type));
					/* END - linker */
					
					gzprintf(outfd, "\t%s\t%s", linkerClassToStr(pCPU1->finalLinkerType), outputClassToStr(pCPU1->classification));
				} else {
					gzprintf(outfd, "%d\t%d\t%d\t%d\t%d", pCPU1->adapter_ref_begin1, pCPU1->adapter_ref_end1, pCPU1->adapter_read_begin1, pCPU1->adapter_read_end1, pCPU1->adapter_score1);
					gzprintf(outfd, "\t%#x", pCPU1->adapter_type);
					
					gzprintf(outfd, "\t%d\t%d\t%d\t%d\t%d", pCPU2->adapter_ref_begin1, pCPU2->adapter_ref_end1, pCPU2->adapter_read_begin1, pCPU2->adapter_read_end1, pCPU2->adapter_score1);
					gzprintf(outfd, "\t%#x", pCPU2->adapter_type);
					/* END - adapter */
					
					/* linker */
					gzprintf(outfd, "\t%d\t%d\t%d\t%d\t%d", pCPU1->linker_ref_begin1, pCPU1->linker_ref_end1, pCPU1->linker_read_begin1, pCPU1->linker_read_end1, pCPU1->linker_score1);
					gzprintf(outfd, "\t%#x\t%#x", pCPU1->linker_type, pCPU1->linker_type);
					gzprintf(outfd, "\t%d\t%d\t%d\t%d\t%d", pCPU2->linker_ref_begin1, pCPU2->linker_ref_end1, pCPU2->linker_read_begin1, pCPU2->linker_read_end1, pCPU2->linker_score1);
					gzprintf(outfd, "\t%#x\t%#x", pCPU2->linker_type, pCPU2->linker_type);
					/* END - linker */
					
					gzprintf(outfd, "\t%#x\t%#x", pCPU1->finalLinkerType, pCPU1->classification);
				}

				gzprintf(outfd, "\t%d\t%d\t%d\t%d", pCPU1->len[0], pCPU1->len[1], pCPU2->len[0], pCPU2->len[1]);
				
				gzprintf(outfd, "\t%s", seqs[i].name);
				
				if (opt->cpuflag & CPU_SEQ) {
					gzprintf(outfd, "\t%s\t%s", seqs[i].seq, seqs[i+1].seq);
				}
				if (opt->cpuflag & CPU_QUALITY) {
					gzprintf(outfd, "\t%s\t%s", seqs[i].qual, seqs[i+1].qual);
				}
				
				gzprintf(outfd, "\n");
			}
		}
		
		// additional output
		if (opt->cpuflag & CPU_FASTQ) {
			for (i = 0; i < n; ++i)
				if (cpu_fqs[i].fq.l>0) gzwrite(outfds[cpus[i].classification], cpu_fqs[i].fq.s, cpu_fqs[i].fq.l);
		}
		// END - OUTPUT
		
		free(cpus);
		
		for (i = 0; i < n; ++i) {
			//err_fputs(seqs[i].sam, stdout);
			free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual);
			free(cpu_fqs[i].fq.s);
			//free(seqs[i].sam);
		}
		free(cpu_fqs);
		free(seqs);
		
		n_processed += n;
		if (bwa_verbose >= 3) {
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s] %lld sequences processed, %.0f seq/sec, %.2f min..\n", __func__, n_processed, 1.0*n_processed/t_diff, t_diff/60.0);
		}
	}
	
	int status = 0;
	
	if (gzflush(outfd, Z_FINISH)) {
		fprintf(stderr, "[E::%s] fail to flush stdout for gzip stream: %s\n", __func__, strerror (errno));
		status = 1;
	}
	// TODO: for stdout, we do NOT call gzclose
	
	terminate_CPTags_Ouputs(outfds);
	
	if (gzclose(outfd)) {
		fprintf(stderr, "[E::%s] fail to close CPTags data stream: %s\n", __func__, strerror (errno));
	}

	lmp_opt_terminate(opt);
	free(opt);
	kseq_destroy(ks);
	err_gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	
	return status;
}


