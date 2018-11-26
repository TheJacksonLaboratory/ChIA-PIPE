#ifndef CPU_SPAN_H
#define CPU_SPAN_H

#include "bam.h"

static inline int is_overlap(int32_t s1, int32_t e1, int32_t s2, int32_t e2) {
	if (e1<s2) {
		return 0;
	} else if (e2<s1) {
		return 0;
	}
	return -1;
}

static int getIntraSpan(int32_t s1, int32_t e1, int32_t s2, int32_t e2) {
	int32_t pos_1 = (s1 < s2) ? s1 : s2;
	int32_t pos_2 = (e1 > e2) ? e1 : e2;
	int32_t tlen = ((pos_1>pos_2) ? pos_1 - pos_2 : pos_2 - pos_1) + 1;
	return tlen;
}

static int getIntraSpan_MidPoint(int32_t s1, int32_t e1, int32_t s2, int32_t e2) {
	int32_t center_1 = (s1+e1)/2;
	int32_t center_2 = (s2+e2)/2;
	int32_t tlen = ((center_1>center_2) ? center_1 - center_2 : center_2 - center_1) + 1;
	return tlen;
}

#define PF_STATES_PASS	0
#define PF_STATES_FAIL	1

typedef struct {
	int minMAPQ;
	float minFractionalMAPQ;
	int maxX1;
} pairqc_opt_t;

static int getPairQCStateByMapQ (const bam1_t *bam, const bam1_t *mbam, const pairqc_opt_t *opt) {
	return (bam->core.qual>=opt->minMAPQ && mbam->core.qual>=opt->minMAPQ) ? PF_STATES_PASS : PF_STATES_FAIL;
}


static int getQCStateByUniqueLoci (const bam1_t *bam, const pairqc_opt_t *opt) {
	// aln: {1==X0, 0==X1} strictly for now... will have to determine if there is additional case to relax
	// mem: {1==X0, 0==X1} || {1==X0, X1>0, Y1/Y0>=0.95}
	
	uint8_t *xl = bam_aux_get(bam, G_MAPPROGRAM_TAG);
	if (0!=xl) {
		char *a_xl = bam_aux2Z(xl);
		if (0==strcmp(a_xl, G_MAPPROGRAM_TAG_BWAALN)) {
			uint8_t *hitTag = bam_aux_get(bam, G_TOPHIT_TAG);
			int32_t x0 = bam_aux2i(hitTag);
			if (1!=x0) return PF_STATES_FAIL;
			
			hitTag = bam_aux_get(bam, G_ALTERNATIVEHIT_TAG);
			int32_t x1 = bam_aux2i(hitTag);
			if (x1>opt->maxX1) return PF_STATES_FAIL;
			
		} else if (0==strcmp(a_xl, G_MAPPROGRAM_TAG_BWAMEM)) {
			uint8_t *hitTag = bam_aux_get(bam, G_TOPHIT_TAG);
			int32_t x0 = bam_aux2i(hitTag);
			if (1!=x0) return PF_STATES_FAIL;
			
			hitTag = bam_aux_get(bam, G_ALTERNATIVEHIT_TAG);
			int32_t x1 = bam_aux2i(hitTag);
			if (x1>opt->maxX1) {
				// check that Y1/Y0>=0.95
				hitTag = bam_aux_get(bam, G_ALTERNATIVEHIT_SCORE_TAG);
				int32_t y1 = bam_aux2i(hitTag);
				if (y1>0) {
					hitTag = bam_aux_get(bam, G_TOPHIT_SCORE_TAG);
					int32_t y0 = bam_aux2i(hitTag);
					float fractionalScore = (y1*1.0f) / (y0*1.0f);
					if (!(fractionalScore>opt->minFractionalMAPQ)) return PF_STATES_FAIL;
				}
			}
			
		} else {
			fprintf(stderr, "[E::%s] Unrecognized mapper %s in %s tag for \"%s\".\n", __func__, a_xl, G_MAPPROGRAM_TAG, bam1_qname(bam));
			return PF_STATES_FAIL;
		}
	} else {
		// data has not been processed by CPU!!!
		// TODO: we might have to a lenient mode for passing through
		//       but in this case, we do not know the scoring mechanism if there is no alignment program specified
		//if (bwa_verbose >= 1)
		fprintf(stderr, "[E::%s] %s tag unavailable for \"%s\". Make sure that mapper produce %s tag like CPU..\n", __func__, G_MAPPROGRAM_TAG, bam1_qname(bam), G_MAPPROGRAM_TAG);
		return PF_STATES_FAIL;
	}
	
	
	return PF_STATES_PASS;
}

static int getPairQCStateByUniqueLoci (const bam1_t *bam, const bam1_t *mbam, const pairqc_opt_t *opt) {
	// aln: {1==X0, 0==X1} strictly for now... will have to determine if there is additional case to relax
	// mem: {1==X0, 0==X1} || {1==X0, X1>0, Y1/Y0>=0.95}
	
	int pairQCState = getQCStateByUniqueLoci (bam, opt);
	if (PF_STATES_PASS==pairQCState) {
		pairQCState = getQCStateByUniqueLoci (mbam, opt);
	}
	
	return pairQCState;
}

#endif
