/*******************************************************************************
 ** PENDING:
 1) determine if we have a mutual mode or hybrid mode
    mutual mode: for tagLen <= x (short tag), use bwa aln logic
                 for tagLen >x (non-short tag), use bwa mem logic
 2) determine the effect of chunksize and speed
*******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "kstring.h"
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"

/* bwa aln */
#include "bwtaln.h"
#include "bwtgap.h"
#define MAX_ALN_FIRST_SHORT_TAG_LEN 25
#define MAX_ALN_SHORT_TAG_LEN 69
/* END - bwa aln */
/* bwa se */
#include "bntseq.h"
/* END - bwa se */

#define WRITE_SAM_LOGIC_TAG


#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include "kseq.h"
#include "utils.h"
KSEQ_DECLARE(gzFile)


extern double G_t_real;

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

/***************************
 * this is really bwamem.c + bwtaln.c + bwase.c
 ***************************/

static const bntseq_t *global_bns = 0; // for debugging only

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

/***************************
 * Collection SA invervals *
 ***************************/

//#define intv_lt(a, b) ((a).info < (b).info)
//KSORT_INIT(mem_intv, bwtintv_t, intv_lt)

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

static smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
	return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	const mem_pestat_t *pes;
	const gap_opt_t *optAln; /* bwa aln */
	int n_occ; /* bwa se */
	smem_aux_t **aux;
	gap_stack_t **stacks; /* bwa aln */
	bseq1_t *seqs;
	mem_alnreg_v *regs;
	int64_t n_processed;
	bwa_seq_t *seqsAln; /* bwa aln */
} worker_t;

//bwamem.c
static inline int get_rlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; ++k) {
		int op = cigar[k]&0xf;
		if (op == 0 || op == 2)
			l += cigar[k]>>4;
	}
	return l;
}

//bwamem.c
void memaln_aln2sam(const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_, int softclip_all
					, int X0, int X0Score, int X1, int X1Score) // for chia-pet
{
	int i, l_name;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert
	
	if (m_) mtmp = *m_, m = &mtmp;
	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand
	
	// print up to CIGAR
	l_name = strlen(s->name);
	ks_resize(str, str->l + s->l_seq + l_name + (s->qual? s->l_seq : 0) + 20);
	kputsn(s->name, l_name, str); kputc('\t', str); // QNAME
	kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
	if (p->rid >= 0) { // with coordinate
		kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
		kputl(p->pos + 1, str); kputc('\t', str); // POS
		kputw(p->mapq, str); kputc('\t', str); // MAPQ
		if (p->n_cigar) { // aligned
			for (i = 0; i < p->n_cigar; ++i) {
				int c = p->cigar[i]&0xf;
				if (!softclip_all && (c == 3 || c == 4))
					c = which? 4 : 3; // use hard clipping for supplementary alignments
				kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
			}
		} else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
	} else kputsn("*\t0\t0\t*", 7, str); // without coordinte
	kputc('\t', str);
	
	// print the mate position if applicable
	if (m && m->rid >= 0) {
		if (p->rid == m->rid) kputc('=', str);
		else kputs(bns->anns[m->rid].name, str);
		kputc('\t', str);
		kputl(m->pos + 1, str); kputc('\t', str);
		if (p->rid == m->rid) {
			int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
			int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
			if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
			else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
		} else kputc('0', str);
	} else kputsn("*\t0\t0", 5, str);
	kputc('\t', str);
	
	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, str);
	} else if (!p->is_rev) { // the forward strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !softclip_all) { // have cigar && not the primary alignment && not softclip all
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	} else { // the reverse strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !softclip_all) {
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	}
	
	// print optional tags
	if (p->n_cigar) {
		if (X0>=0) { kputsn("\tXT:A:", 6, str); kputc((1==X0)?'U':'R', str); }
		kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
		// additional tag for chia-pet
		if (X0>=0) { kputsn("\tX0:i:", 6, str); kputw(X0, str); kputsn("\tY0:i:", 6, str); kputw(X0Score, str); }
		if (X1>=0) { kputsn("\tX1:i:", 6, str); kputw(X1, str); kputsn("\tY1:i:", 6, str); kputw(X1Score, str); }
		// END - additional tag for chia-pet
		kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
	}
	if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
	if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
	if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
			kputsn("\tSA:Z:", 6, str);
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (list[i].flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
				kputs(bns->anns[r->rid].name, str); kputc(',', str);
				kputl(r->pos+1, str); kputc(',', str);
				kputc("+-"[r->is_rev], str); kputc(',', str);
				for (k = 0; k < r->n_cigar; ++k) {
					kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
				}
				kputc(',', str); kputw(r->mapq, str);
				kputc(',', str); kputw(r->NM, str);
				kputc(';', str);
			}
		}
	}
	if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str); }
	if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
	
#ifdef WRITE_SAM_LOGIC_TAG
	kputs("\tXL:Z:mem", str);
#endif
	kputc('\n', str);
}

void memaln_regX0X1count (const mem_opt_t *opt, mem_alnreg_v *a, int entryK, int *X0, int *X0Score, int *X1, int *X1Score)
{
	int k;
	int x1ScoreIndex=-1;
	
	*X0 = 0; *X0Score = a->a[entryK].score;
	*X1 = 0; *X1Score = 0;
	
	for (k = entryK; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		if (*X0Score==p->score) {
			*X0 += 1;
		} else if (p->score<*X0Score) {
			x1ScoreIndex = k;
			*X1Score = p->score;
			break;
		}
	}
	if (-1!=x1ScoreIndex) {
		for (k = x1ScoreIndex; k < a->n; ++k) {
			mem_alnreg_t *p = &a->a[k];
			if (*X1Score==p->score) {
				*X1 += 1;
			} else {
				break;
			}
		}
	}
}

// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible
void memaln_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m)
{
	extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, mem_alnreg_v *a, int l_query, const char *query);
	
	kstring_t str;
	kvec_t(mem_aln_t) aa;
	int k;
	char **XA = 0;
	int entryK = -1;
	
	if (!(opt->flag & MEM_F_ALL))
		XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
	kv_init(aa);
	str.l = str.m = 0; str.s = 0;
	for (k = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;
		if (p->secondary >= 0 && !(opt->flag&MEM_F_ALL)) continue;
		if (p->secondary >= 0 && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
		if (-1==entryK) entryK = k;
		q = kv_pushp(mem_aln_t, aa);
		*q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
		assert(q->rid >= 0); // this should not happen with the new code
		q->XA = XA? XA[k] : 0;
		q->flag |= extra_flag; // flag secondary
		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
		if (k && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
		if (k && q->mapq > aa.a[0].mapq) q->mapq = aa.a[0].mapq;
		if (bwa_verbose >= 4) printf("* memaln_reg2sam_se: %s score=%d, truesc=%d, sub=%d, csub=%d, sub_n=%d, w=%d, n_cmp=%d\n", s->name, p->score, p->truesc, p->sub, p->csub, p->sub_n, p->w, p->n_comp);
	}
	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
		mem_aln_t t;
		t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
		t.flag |= extra_flag;
		memaln_aln2sam(bns, &str, s, 1, &t, 0, m, opt->flag&MEM_F_SOFTCLIP, -1, -1, -1, -1);
	} else {
		int X0, X0Score, X1, X1Score;
		memaln_regX0X1count (opt, a, entryK, &X0, &X0Score, &X1, &X1Score);
		for (k = 0; k < aa.n; ++k)
			memaln_aln2sam(bns, &str, s, aa.n, aa.a, k, m, opt->flag&MEM_F_SOFTCLIP, 0==k?X0:-1, 0==k?X0Score:-1, 0==k?X1:-1, 0==k?X1Score:-1);
		for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
		free(aa.a);
	}
	s->sam = str.s;
	if (XA) {
		for (k = 0; k < a->n; ++k) free(XA[k]);
		free(XA);
	}
}

/* bwa aln */
// adapted from bwtaln.c: bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
// original form loop thru' every elements skip by n_thread
// new form compute a single element so that looping can be done in different ways
// gap_stack initialization is done by caller as this only handle a single element, also we know the max length anyway "MAX_SHORT_TAG_LEN"
void bwa_cal_sa_reg_gap_single(const bwt_t *bwt, bwa_seq_t *p, const gap_opt_t *opt,
							   gap_stack_t *stack
							   )
{
	//bwtaln.c
	extern int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);
	
	int j, max_l = 0;
	bwt_width_t *w, *seed_w;
	gap_opt_t local_opt = *opt;
	
	seed_w = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	w = 0;
	
	p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
	if (max_l < p->len) {
		max_l = p->len;
		w = (bwt_width_t*)realloc(w, (max_l + 1) * sizeof(bwt_width_t));
		memset(w, 0, (max_l + 1) * sizeof(bwt_width_t));
	}
	bwt_cal_width(bwt, p->len, p->seq, w);
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
	local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
	if (p->len > opt->seed_len)
		bwt_cal_width(bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), seed_w);
	// core function
	for (j = 0; j < p->len; ++j) // we need to complement
		p->seq[j] = p->seq[j] > 3? 4 : 3 - p->seq[j];
	p->aln = bwt_match_gap(bwt, p->len, p->seq, w, p->len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, stack);
	//fprintf(stderr, "mm=%lld,ins=%lld,del=%lld,gapo=%lld\n", p->aln->n_mm, p->aln->n_ins, p->aln->n_del, p->aln->n_gapo);
	// clean up the unused data in the record
	// TOOD: we don't need to clean up; preparing the sam record immediate after this
	free(p->name);
	free(p->seq);
	free(p->rseq);
	free(p->qual);
	p->name = 0; p->seq = p->rseq = p->qual = 0;
	
	free(seed_w); free(w);
}
/* END - bwa aln */

/* TODO: do we perform local clean up instead
void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs)
{
	int i, j;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
		for (j = 0; j < p->n_multi; ++j)
			if (p->multi[j].cigar) free(p->multi[j].cigar);
		free(p->name);
		free(p->seq); free(p->rseq); free(p->qual); free(p->aln); free(p->md); free(p->multi);
		free(p->cigar);
	}
	free(seqs);
}
*/

static void memaln_short_tag_align1_core (const gap_opt_t *optAln, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, const bseq1_t *seq, bwa_seq_t *p, gap_stack_t *stack, int n_occ)
{
	// bwase.c
	extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
	//extern void bwa_cal_pac_pos(const bntseq_t *bns, const char *prefix, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr);
	extern void bwa_cal_pac_pos_core(const bntseq_t *bns, const bwt_t *bwt, bwa_seq_t *seq, const int max_mm, const float fnr);
	extern bwtint_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int ref_len, int *strand);
	extern void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq);
	// bwaseqio.c
	int bwa_trim_read(int trim_qual, bwa_seq_t *p);
	
	int j;
	int trim_qual = 0; // TODO: pass this in or this is already handled by caller?
	
	// TODO: perform bwa aln
	// copied from bwaseqio.c (invoked via bwa_read_seq(..) in bwa aln)
	p->tid = -1; // no assigned to a thread
	p->qual = 0;
	p->full_len = p->clip_len = p->len = seq->l_seq;
	p->seq = (ubyte_t*)calloc(p->full_len, 1);
	for (j = 0; j != p->full_len; ++j)
		p->seq[j] = nst_nt4_table[(int)seq->seq[j]];
	if (seq->qual) { // copy quality
		p->qual = (ubyte_t*)strdup((char*)seq->qual);
		if (trim_qual >= 1) bwa_trim_read(trim_qual, p);
	}
	//XCode analysis: 0 length allocation
	//p->rseq = (ubyte_t*)calloc(p->full_len, 1);
	p->rseq = (ubyte_t*)calloc(p->full_len>0?p->full_len:1, 1);
	memcpy(p->rseq, p->seq, p->len);
	seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
	// TODO: int is_comp = mode&BWA_MODE_COMPREAD;
	int is_comp = optAln->mode&BWA_MODE_COMPREAD;
	seq_reverse(p->len, p->rseq, is_comp);
	
	// TODO: what's the difference between "const bwt_t" * vs "bwt *const"
	bwa_cal_sa_reg_gap_single(bwt, p, optAln, stack); // p-seq got restored
	
	// TODO: perform bwa samse (bwase.c)
	// read alignment
	bwa_aln2seq_core(p->n_aln, p->aln, p, 1, n_occ);
	
	// TODO: [bwa_aln_core] convert to sequence coordinate...
	//bwa_cal_pac_pos(bns, prefix, 1, p, opt.max_diff, w->optAln.fnr); // forward bwt will be destroyed here
	bwa_cal_pac_pos_core(bns, bwt, p, optAln->max_diff, optAln->fnr);
	int strand, n_multi;
	for (j = n_multi = 0; j < p->n_multi; ++j) {
		bwt_multi1_t *q = p->multi + j;
		q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len + q->ref_shift, &strand);
		q->strand = strand;
		if (q->pos != p->pos && q->pos != (bwtint_t)-1)
			p->multi[n_multi++] = *q;
	}
	p->n_multi = n_multi;
	
	// TODO: [bwa_aln_core] refine gapped alignments...
	// BLOCK : repeat as bwa_refined_gapped() needs the data
	p->seq = (ubyte_t*)calloc(p->full_len, 1);
	for (j = 0; j != p->full_len; ++j)
		p->seq[j] = nst_nt4_table[(int)seq->seq[j]];
	if (seq->qual) { // copy quality
		p->qual = (ubyte_t*)strdup((char*)seq->qual);
		if (trim_qual >= 1) bwa_trim_read(trim_qual, p);
	}
	p->rseq = (ubyte_t*)calloc(p->full_len, 1);
	memcpy(p->rseq, p->seq, p->len);
	seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
	//is_comp = optAln->mode&BWA_MODE_COMPREAD;;
	seq_reverse(p->len, p->rseq, is_comp);
	// we need the name for printing sam
	p->name = strdup(seq->name);
	{ // trim /[12]$
		size_t t = strlen(p->name);
		if (t > 2 && p->name[t-2] == '/' && (p->name[t-1] == '1' || p->name[t-1] == '2')) p->name[t-2] = '\0';
	}
	// END-BLOCK
	//bwa_refine_gapped(w->bns, 1, p, 0);
	//seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
	bwa_refine_gapped(bns, 1, p, pac);
}

static void memaln_worker1(void *data, int i, int tid)
{
	// bwamem.c
	extern mem_alnreg_v mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, void *buf);

	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
		if (w->seqs[i].l_seq>MAX_ALN_SHORT_TAG_LEN) {
			// TODO: parameter tuning back to normal
			w->regs[i] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, w->aux[tid]);
		} else if (w->seqs[i].l_seq>MAX_ALN_FIRST_SHORT_TAG_LEN) {
			// length = (MAX_ALN_FIRST_SHORT_TAG_LEN, MAX_ALN_SHORT_TAG_LEN]
			memaln_short_tag_align1_core(w->optAln, w->bwt, w->bns, w->pac, &(w->seqs[i]), &(w->seqsAln[i]), w->stacks[tid], w->n_occ);
			if (BWA_TYPE_NO_MATCH==w->seqsAln[i].type) {
				// TODO: parameter tuning for short read
				/*mem_opt_t *opt = mem_opt_init();
				opt->a = 1; opt->b = 2; opt->T = 15;
				opt->o_del = opt->o_ins = 6;
				opt->e_del = opt->e_ins = 1;
				opt->pen_clip5 = opt->pen_clip3 = 0;
				bwa_fill_scmat(opt->a, opt->b, opt->mat);
				w->regs[i] = mem_align1_core(opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, w->aux[tid]);
				// TODO: one-time initialization
				free(opt);*/
				w->regs[i] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, w->aux[tid]);
			}
		} else {
			// length = (0, MAX_ALN_FIRST_SHORT_TAG_LEN]
			memaln_short_tag_align1_core(w->optAln, w->bwt, w->bns, w->pac, &(w->seqs[i]), &(w->seqsAln[i]), w->stacks[tid], w->n_occ);
			
			// TODO: handle 18bp and 19bp specific
			// 18bp = 18 exact match & unique, else pointless
			// 19bp = allow 1 mismatch & unique, else pointless
		}
	} else {
		// TODO: not implemented!!!
		if (bwa_verbose >= 4) printf("=====> Processing read '%s'/1 <=====\n", w->seqs[i<<1|0].name);
		w->regs[i<<1|0] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i<<1|0].l_seq, w->seqs[i<<1|0].seq, w->aux[tid]);
		if (bwa_verbose >= 4) printf("=====> Processing read '%s'/2 <=====\n", w->seqs[i<<1|1].name);
		w->regs[i<<1|1] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i<<1|1].l_seq, w->seqs[i<<1|1].seq, w->aux[tid]);
	}
}

// from bwase.c: bwa_print_seq
void memaln_print_seq(kstring_t *s, bwa_seq_t *seq) {
	char buffer[4096];
	const int bsz = sizeof(buffer);
	int i, j, l;
	
	if (seq->strand == 0) {
		for (i = 0; i < seq->full_len; i += bsz) {
			l = seq->full_len - i > bsz ? bsz : seq->full_len - i;
			for (j = 0; j < l; j++) buffer[j] = "ACGTN"[seq->seq[i + j]];
			//err_fwrite(buffer, 1, l, stream);
			kputsn(buffer, l, s);
		}
	} else {
		for (i = seq->full_len - 1; i >= 0; i -= bsz) {
			l = i + 1 > bsz ? bsz : i + 1;
			for (j = 0; j < l; j++) buffer[j] = "TGCAN"[seq->seq[i - j]];
			//err_fwrite(buffer, 1, l, stream);
			kputsn(buffer, l, s);
		}
	}
}

// from bwase.c: bwa_print_sam1
void memaln_bwaseq2sam_se(const bntseq_t *bns, bwa_seq_t *p, bseq1_t *s, int mode, int max_top2)
{
	int i, j;
	kstring_t str = {0,0,0};
	
	if (BWA_TYPE_NO_MATCH==p->type) {
		// this read has no match
		
		int flag = p->extra_flag | SAM_FSU;
		ksprintf(&str, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
		
		//bwa_print_seq(stdout, p);
		memaln_print_seq(&str, p);
		
		kputc('\t', &str);
		if (p->qual) {
			if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
			kputs((const char*)p->qual, &str);
		} else kputc('*', &str);
		if (bwa_rg_id[0]) ksprintf(&str, "\tRG:Z:%s", bwa_rg_id);
		if (p->bc[0]) ksprintf(&str, "\tBC:Z:%s", p->bc);
		if (p->clip_len < p->full_len) ksprintf(&str, "\tXC:i:%d", p->clip_len);
		kputc('\n', &str);
	} else {
		// TODO: len or full_len?
		if ((18==p->full_len && 0==p->nm) || (19==p->full_len && 2>p->nm) || p->full_len>=20) {
			int seqid, nn, flag = p->extra_flag;
			char XT;
			
			j = pos_end(p) - p->pos; // j is the length of the reference in the alignment
			
			// get seqid
			nn = bns_cnt_ambi(bns, p->pos, j, &seqid);
			if (p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
				flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences
			
			// update flag and print it
			if (p->strand) flag |= SAM_FSR;
			ksprintf(&str, "%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
			ksprintf(&str, "%d\t%d\t", (int)(p->pos - bns->anns[seqid].offset + 1), p->mapQ);
			
			// print CIGAR
			if (p->cigar) {
				for (j = 0; j != p->n_cigar; ++j)
					ksprintf(&str, "%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
			} else if (p->type == BWA_TYPE_NO_MATCH) kputc('*', &str);
			else ksprintf(&str, "%dM", p->len);
			
			ksprintf(&str, "\t*\t0\t0\t");
			
			// print sequence and quality
			//bwa_print_seq(stdout, p);
			memaln_print_seq(&str, p);
			
			kputc('\t', &str);
			if (p->qual) {
				if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
				ksprintf(&str, "%s", p->qual);
			} else kputc('*', &str);
			
			if (bwa_rg_id[0]) ksprintf(&str, "\tRG:Z:%s", bwa_rg_id);
			if (p->bc[0]) ksprintf(&str, "\tBC:Z:%s", p->bc);
			if (p->clip_len < p->full_len) ksprintf(&str, "\tXC:i:%d", p->clip_len);
			
			// calculate XT tag
			XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';
			// print tags
			ksprintf(&str, "\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD)? "NM" : "CM", p->nm);
			if (nn) ksprintf(&str, "\tXN:i:%d", nn);
			if (p->type != BWA_TYPE_MATESW) { // X0 and X1 are not available for this type of alignment
				ksprintf(&str, "\tX0:i:%d", p->c1);
				if (p->c1 <= max_top2) ksprintf(&str, "\tX1:i:%d", p->c2);
			}
			ksprintf(&str, "\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
			if (p->md) ksprintf(&str, "\tMD:Z:%s", p->md);
			// print multiple hits
			if (p->n_multi) {
				ksprintf(&str, "\tXA:Z:");
				for (i = 0; i < p->n_multi; ++i) {
					bwt_multi1_t *q = p->multi + i;
					int k;
					j = pos_end_multi(q, p->len) - q->pos;
					nn = bns_cnt_ambi(bns, q->pos, j, &seqid);
					ksprintf(&str, "%s,%c%d,", bns->anns[seqid].name, q->strand? '-' : '+',
							 (int)(q->pos - bns->anns[seqid].offset + 1));
					if (q->cigar) {
						for (k = 0; k < q->n_cigar; ++k)
							ksprintf(&str, "%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
					} else ksprintf(&str, "%dM", p->len);
					ksprintf(&str, ",%d;", q->gap + q->mm);
				}
			}
			
#ifdef WRITE_SAM_LOGIC_TAG
			kputs("\tXL:Z:aln", &str);
#endif
			kputc('\n', &str);
		} else {
			// no a valid mapping anyway; treated as no match
			int flag = p->extra_flag | SAM_FSU;
			ksprintf(&str, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
			
			//bwa_print_seq(stdout, p);
			memaln_print_seq(&str, p);
			
			kputc('\t', &str);
			if (p->qual) {
				if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
				kputs((const char*)p->qual, &str);
			} else kputc('*', &str);
			if (bwa_rg_id[0]) ksprintf(&str, "\tRG:Z:%s", bwa_rg_id);
			if (p->bc[0]) ksprintf(&str, "\tBC:Z:%s", p->bc);
			if (p->clip_len < p->full_len) ksprintf(&str, "\tXC:i:%d", p->clip_len);
			kputc('\n', &str);
		}
	}
	s->sam = str.s;
}


static void memaln_worker2(void *data, int i, int tid)
{
	extern int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);
	extern void mem_reg2ovlp(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a);
	// bwamem.c
	extern void mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
	// bwamem.c BUT will have to be changed to our version
	//extern void mem_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);

	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (bwa_verbose >= 4) printf("=====> Finalizing read '%s' <=====\n", w->seqs[i].name);
#if 1
		if (w->seqs[i].l_seq>MAX_ALN_SHORT_TAG_LEN) {
			if (w->opt->flag & MEM_F_ALN_REG) {
				mem_reg2ovlp(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i]);
			} else {
				mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
				memaln_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
			}
		} else if (w->seqs[i].l_seq>MAX_ALN_FIRST_SHORT_TAG_LEN) {
			// length = (MAX_ALN_FIRST_SHORT_TAG_LEN, MAX_ALN_SHORT_TAG_LEN]
			if (w->regs[i].n>0) {
				if (w->opt->flag & MEM_F_ALN_REG) {
					mem_reg2ovlp(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i]);
				} else {
					// TODO: parameter tuning for short read
					/*mem_opt_t *opt = mem_opt_init();
					opt->a = 1; opt->b = 2; opt->T = 15;
					opt->o_del = opt->o_ins = 6;
					opt->e_del = opt->e_ins = 1;
					bwa_fill_scmat(opt->a, opt->b, opt->mat);
					mem_mark_primary_se(opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
					memaln_reg2sam_se(opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
					// TODO: one-time initialization
					free(opt);*/
					mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
					memaln_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
				}
			} else {
				// report bwa aln result
				// [bwa_aln_core] print alignments...; should be worker2 ???
				// bwa_print_sam1(w->bns, p, 0, w->optAln->mode, w->optAln->max_top2);
				// load sam records for external printing
				memaln_bwaseq2sam_se(w->bns, &(w->seqsAln[i]), &w->seqs[i], w->optAln->mode, w->optAln->max_top2);
				// we do not release current single element as the code has the assumption to free the array
				// bwa_free_read_seq(1, p);
			}
		} else {
			// length = (0, MAX_ALN_FIRST_SHORT_TAG_LEN]
			
			// report bwa aln result
			// [bwa_aln_core] print alignments...; should be worker2 ???
			// bwa_print_sam1(w->bns, p, 0, w->optAln->mode, w->optAln->max_top2);
			// load sam records for external printing
			memaln_bwaseq2sam_se(w->bns, &(w->seqsAln[i]), &w->seqs[i], w->optAln->mode, w->optAln->max_top2);
			// we do not release current single element as the code has the assumption to free the array
			// bwa_free_read_seq(1, p);
		}
#else
		if (w->regs[i].n>0) {
			if (w->opt->flag & MEM_F_ALN_REG) {
				mem_reg2ovlp(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i]);
			} else {
				mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
				memaln_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
			}
		} else {
			// report bwa aln result
			// [bwa_aln_core] print alignments...; should be worker2 ???
			// bwa_print_sam1(w->bns, p, 0, w->optAln->mode, w->optAln->max_top2);
			// load sam records for external printing
			memaln_bwaseq2sam_se(w->bns, &(w->seqsAln[i]), &w->seqs[i], w->optAln->mode, w->optAln->max_top2);
			// we do not release current single element as the code has the assumption to free the array
			// bwa_free_read_seq(1, p);
		}
#endif
		free(w->regs[i].a);
	} else {
		// TODO: not implemented!!!
		if (bwa_verbose >= 4) printf("=====> Finalizing read pair '%s' <=====\n", w->seqs[i<<1|0].name);
		mem_sam_pe(w->opt, w->bns, w->pac, w->pes, (w->n_processed>>1) + i, &w->seqs[i<<1], &w->regs[i<<1]);
		free(w->regs[i<<1|0].a); free(w->regs[i<<1|1].a);
	}
}

void memaln_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0,
	/* bwa aln */ const gap_opt_t *optAln, /* bwa se */ int n_occ)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	mem_alnreg_v *regs;
	mem_pestat_t pes[4];
	double ctime, rtime;
	int i;
	gap_opt_t localOptAln;
	
	ctime = cputime(); rtime = realtime();
	global_bns = bns;
	//regs = malloc(n * sizeof(mem_alnreg_v));
	regs = calloc(n, sizeof(mem_alnreg_v));
	w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
	w.seqs = seqs; w.regs = regs; w.n_processed = n_processed; // TODO:
	//TODO: w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	w.pes = &pes[0];
	w.optAln = optAln; /* bwa aln */
	w.n_occ = n_occ; /* bwa se */
	w.aux = malloc(opt->n_threads * sizeof(smem_aux_t));
	/* bwa aln */
	w.stacks = malloc(opt->n_threads * sizeof(gap_stack_t)); // TODO: check if this is allocated correctly!!!
	localOptAln = *optAln;
	if (optAln->fnr > 0.0) localOptAln.max_diff = bwa_cal_maxdiff(MAX_ALN_SHORT_TAG_LEN, BWA_AVG_ERR, optAln->fnr);
	if (localOptAln.max_diff < localOptAln.max_gapo) localOptAln.max_gapo = localOptAln.max_diff;
	w.seqsAln = calloc(n, sizeof(bwa_seq_t));
	/* END - bwa aln */
	
	for (i = 0; i < opt->n_threads; ++i) {
		w.aux[i] = smem_aux_init();
		/* bwa aln */
		// initiate priority stack
		w.stacks[i] = gap_init_stack(localOptAln.max_diff, localOptAln.max_gapo, localOptAln.max_gape, &localOptAln);
		/* END - bwa aln */
	}
	kt_for(opt->n_threads, memaln_worker1, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
	for (i = 0; i < opt->n_threads; ++i) {
		smem_aux_destroy(w.aux[i]);
		gap_destroy_stack(w.stacks[i]); /* bwa aln */
	}
	free(w.aux);
	free(w.stacks); /* bwa aln */
	if (opt->flag&MEM_F_PE) { // infer insert sizes if not provided
		if (pes0) memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
		else mem_pestat(opt, bns->l_pac, n, regs, pes); // otherwise, infer the insert size distribution from data
	}
	kt_for(opt->n_threads, memaln_worker2, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // generate alignment
	/* bwa se */
	bwa_free_read_seq(n, w.seqsAln);
	/* END - bwa se */
	free(regs);
	if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}

int main_memaln(int argc, char *argv[])
{
	/* bwa se */
	extern void bwase_initialize();
	
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, n, copy_comment = 0;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	bwaidx_t *idx;
	char *p, *rg_line = 0;
	const char *mode = 0;
	void *ko = 0, *ko2 = 0;
	int64_t n_processed = 0;
	mem_pestat_t pes[4], *pes0 = 0;
	
	/* bwa aln */
	gap_opt_t *optAln;
	// bwt_t *bwt;
	/* END - bwa aln */
	
	int n_occ = 3; /* bwa se */
	
	double t_diff;
	double t_timepointIO, t_timepointProcessing;
	double t_diffIO, t_diffProcessing;
	
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;
	
	opt = mem_opt_init();
	/* bwa aln */
	optAln = gap_init_opt();
	optAln->seed_len = 20; optAln->max_seed_diff = 2;
	/* END - bwa aln */

	memset(&opt0, 0, sizeof(mem_opt_t));
	while ((c = getopt(argc, argv, "epaFMCSPHYn:k:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == 'x') mode = optarg;
		else if (c == 'n') n_occ = atoi(optarg); /* bwa se */
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		//else if (c == 'p') opt->flag |= MEM_F_PE; // WCH: to be treated as SE
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'e') opt->flag |= MEM_F_SELF_OVLP;
		else if (c == 'F') opt->flag |= MEM_F_ALN_REG;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 'h') opt->max_hits = atoi(optarg), opt0.max_hits = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'C') copy_comment = 1;
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		} else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = (int32_t) strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = (int32_t) strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = (int32_t) strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = (int32_t) strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = (int32_t) strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = (int32_t) strtol(p+1, &p, 10);
		} else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 'I') { // specify the insert size distribution
			pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
						__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else {
			free(opt);
			free(optAln);  /* bwa aln */
			return 1;
		}
	}
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 1 >= argc || optind + 3 < argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu memaln [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		//		fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
		fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "       -e            discard full-length exact matches\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n", opt->pen_unpaired);
		fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
		fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A2 -B5 -O2 -E1 -L0\n");
		fprintf(stderr, "                     pbread: -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -N25 -FeaD.001\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -h INT        if there are <INT hits with score >80%% of the max score, output all in XA [%d]\n", opt->max_hits);
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		free(optAln);  /* bwa aln */
		return 1;
	}
	
	if (mode) {
		if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "pbread1") == 0 || strcmp(mode, "pbread") == 0) {
			if (!opt0.a) opt->a = 2, opt0.a = 1;
			update_a(opt, &opt0);
			if (!opt0.o_del) opt->o_del = 2;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 2;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 5;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
			if (strcmp(mode, "pbread1") == 0 || strcmp(mode, "pbread") == 0) {
				opt->flag |= MEM_F_ALL | MEM_F_SELF_OVLP | MEM_F_ALN_REG;
				if (!opt0.max_occ) opt->max_occ = 1000;
				if (!opt0.min_seed_len) opt->min_seed_len = 13;
				if (!opt0.max_chain_extend) opt->max_chain_extend = 25;
				if (opt0.drop_ratio == 0.) opt->drop_ratio = .001;
			} else {
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			free(opt);
			free(optAln);  /* bwa aln */
			return 1; /*// FIXME memory leak*/
		}
	} else update_a(opt, &opt0);
	//	if (opt->T < opt->min_HSP_score) opt->T = opt->min_HSP_score; // TODO: tie ->T to MEM_HSP_COEF
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	
	 /* bwa aln */
	if (optAln->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, optAln->fnr);
			if (l != k) fprintf(stderr, "[M::%s] %dbp reads: max_diff = %d\n", __func__, i, l);
			k = l;
		}
	}
	/* END - bwa aln */
	
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL )) == 0) {
	//if ((idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS )) == 0) {
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] BWT, BNS and PAC data loading failed for %s\n", __func__, argv[optind]);
		}
		free(opt);
		free(optAln);  /* bwa aln */
		return 1; /*// FIXME: memory leak*/
	}
	if (bwa_verbose >= 3) {
		t_diff = realtime() - G_t_real;
		fprintf(stderr, "[M::%s] BWT, BNS and PAC data loaded in %.2f min..\n", __func__, t_diff/60.0);
	}
	
	/* bwa aln */
	// WCH: we piggyback on bwa_idx_load for what bwa aln needs
	// bwt = idx->bwt;
	/* END - bwa aln */
	/* bwa se */
	// initialization
	bwase_initialize();
	//bns = bns_restore(prefix);
	srand48(idx->bns->seed);
	/* END - bwa se */
	
	ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		return 1;
	}
	fp = gzdopen(fd, "r");
	ks = kseq_init(fp);
	if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file will be ignored.\n", __func__);
		} else {
			ko2 = kopen(argv[optind + 2], &fd2);
			if (ko2 == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
				return 1;
			}
			fp2 = gzdopen(fd2, "r");
			ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}
	if (!(opt->flag & MEM_F_ALN_REG))
		bwa_print_sam_hdr(idx->bns, rg_line);
	t_timepointProcessing = realtime();
	while ((seqs = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
		t_timepointIO = realtime();
		t_diffIO = t_timepointIO - t_timepointProcessing;
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
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, n, (long)size);
		memaln_process_seqs(opt, idx->bwt, idx->bns, idx->pac, n_processed, n, seqs, pes0,
							/* bwa aln */ optAln, /* bwa se */ n_occ);
		t_timepointProcessing = realtime();
		t_diffProcessing = t_timepointProcessing - t_timepointIO;
		
		// TODO: write out .bam instead for .sam
		for (i = 0; i < n; ++i) {
			if (seqs[i].sam) err_fputs(seqs[i].sam, stdout);
		}
		t_timepointIO = realtime();
		t_diffIO += (t_timepointIO - t_timepointProcessing);
		for (i = 0; i < n; ++i) {
			free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
		}
		free(seqs);
		
		n_processed += n;
		if (bwa_verbose >= 3) {
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s] %lld sequences processed, %.0f seq/sec, %.2f min, i/o %.2f sec, processing %.2f sec..\n", __func__, n_processed, 1.0*n_processed/t_diff, t_diff/60.0, t_diffIO, t_diffProcessing);
		}
		
		t_timepointProcessing = realtime();
	}
	
	free(opt);
	free(optAln);  /* bwa aln */
	bwa_idx_destroy(idx);
	kseq_destroy(ks);
	err_gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	return 0;
}

