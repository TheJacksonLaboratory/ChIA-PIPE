#ifndef CPU_TAG_H
#define CPU_TAG_H

#include <errno.h>

#include "kstring.h"

#include "zlib.h"

#include "cpuadapter.h"
#include "cpulinker.h"
#include "cputagcommon.h"

//#define CPU_OUTPUT_NONE        0x0000
//#define CPU_OUTPUT_TIED        0x0001
//#define CPU_OUTPUT_CONFLICT    0x0002
#define CPU_OUTPUT_HL          0x0003
#define CPU_OUTPUT_FL_C_2TAGS  0x0004
#define CPU_OUTPUT_FL_C_1TAG   0x0005
#define CPU_OUTPUT_FL_NC_2TAGS 0x0006
#define CPU_OUTPUT_FL_NC_1TAG  0x0007
#define CPU_OUTPUT_COUNT       0x0008

static char outputClassBuffer[11];
static const char *outputClassToStr(int16_t t) {
	if (CPU_OUTPUT_NONE==t) return "none";
	else if (CPU_OUTPUT_TIED==t) return "tied";
	else if (CPU_OUTPUT_CONFLICT==t) return "conflict";
	else if (CPU_OUTPUT_HL==t) return "half";
	else if (CPU_OUTPUT_FL_C_2TAGS==t) return "full2tc";
	else if (CPU_OUTPUT_FL_C_1TAG==t) return "full1tc";
	else if (CPU_OUTPUT_FL_NC_2TAGS==t) return "full2tnc";
	else if (CPU_OUTPUT_FL_NC_1TAG==t) return "full1tnc";
	else {
		//return "?";
		sprintf(outputClassBuffer, "%#x", t);
		return outputClassBuffer;
	}
}

static char shortOutputClassBuffer[11];
static const char *shortOutputClassToStr(int16_t t) {
	if (CPU_OUTPUT_NONE==t) return "none";
	else if (CPU_OUTPUT_TIED==t) return "tied";
	else if (CPU_OUTPUT_CONFLICT==t) return "confl";
	else if (CPU_OUTPUT_HL==t) return "hl";
	else if (CPU_OUTPUT_FL_C_2TAGS==t) return "fl2tc";
	else if (CPU_OUTPUT_FL_C_1TAG==t) return "fl1tc";
	else if (CPU_OUTPUT_FL_NC_2TAGS==t) return "fl2tnc";
	else if (CPU_OUTPUT_FL_NC_1TAG==t) return "fl1tnc";
	else {
		//return "?";
		sprintf(shortOutputClassBuffer, "%#x", t);
		return shortOutputClassBuffer;
	}
}

static inline int isFLFLConsistent(int l1, int l2) {
	if (LINKER_AA==l1) {
		return (LINKER_AA==l2) ? LINKER_AA : LINKER_NONE;
	} else if (LINKER_BB==l1) {
		return (LINKER_BB==l2) ? LINKER_BB : LINKER_NONE;
	} else if (LINKER_AB==l1) {
		return (LINKER_BA==l2) ? LINKER_AB : LINKER_NONE;
	} else if (LINKER_BA==l1) {
		return (LINKER_AB==l2) ? LINKER_BA : LINKER_NONE;
	}
	
	return LINKER_NONE;
}

static inline int isFLHLConsistent(int l1, int l2) {
	if (LINKER_AA==l1) {
		return (LINKER_A==l2) ? LINKER_AA : LINKER_NONE;
	} else if (LINKER_BB==l1) {
		return (LINKER_B==l2) ? LINKER_BB : LINKER_NONE;
	} else if (LINKER_AB==l1) {
		return (LINKER_B==l2) ? LINKER_AB : LINKER_NONE;
	} else if (LINKER_BA==l1) {
		return (LINKER_A==l2) ? LINKER_BA : LINKER_NONE;
	}
	
	return LINKER_NONE;
}

static inline int isHLFLConsistent(int l1, int l2) {
#if 1
	return isFLHLConsistent(l2, l1);
#else
	if (LINKER_AA==l2) {
		return (LINKER_A==l1) ? LINKER_AA : LINKER_NONE;
	} else if (LINKER_BB==l2) {
		return (LINKER_B==l1) ? LINKER_BB : LINKER_NONE;
	} else if (LINKER_AB==l2) {
		return (LINKER_B==l1) ? LINKER_AB : LINKER_NONE;
	} else if (LINKER_BA==l2) {
		return (LINKER_A==l1) ? LINKER_BA : LINKER_NONE;
	}
	
	return LINKER_NONE;
#endif
}

static inline int isHLHLConsistent(int l1, int l2) {
	if (LINKER_A==l1) {
		return (LINKER_A==l2) ? LINKER_AA : LINKER_NONE;
	} else if (LINKER_B==l1) {
		return (LINKER_B==l2) ? LINKER_BB : LINKER_NONE;
	} else if (LINKER_A==l1) {
		return (LINKER_B==l2) ? LINKER_AB : LINKER_NONE;
	} else if (LINKER_B==l1) {
		return (LINKER_A==l2) ? LINKER_BA : LINKER_NONE;
	}
	
	return LINKER_NONE;
}

// TODO:
// 1) need quality trimming to be included
// 2) need to record how many tags survived for mapping (inject into the id?! better for mapping results processing later)
// 3) currently keep all reads which output still have to check. can be better by keeping array of tag(s) to write
// 4) check if we have linker at the 5' end of a read?

/*
// pairing version of tag extraction
static inline void extractFullFullTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeTwoTagsPosition(pr1, minTagLen);
	computeTwoTagsPosition(pr2, minTagLen);
	
	if (flag & CPU_RAW) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
	} else if (flag & CPU_FASTQ) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		uint16_t totalTags = pr1->tags_n + pr2->tags_n;
		buildTwoTagsFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildTwoTagsFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}

static inline void extractFullHalfTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeTwoTagsPosition(pr1, minTagLen);
	computeTwoTagsPosition(pr2, minTagLen);
	
	if (flag & CPU_RAW) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
	} else if (flag & CPU_FASTQ) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		uint16_t totalTags = pr1->tags_n + pr2->tags_n;
		buildTwoTagsFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildTwoTagsFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}

static inline void extractFullNoneTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeTwoTagsPosition(pr1, minTagLen);
	computeSingleTagPosition(pr2, minTagLen);
	
	if (flag & CPU_RAW) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
	} else if (flag & CPU_FASTQ) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		uint16_t totalTags = pr1->tags_n + pr2->tags_n;
		buildTwoTagsFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildSingleTagFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}

static inline void extractHalfFullTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeTwoTagsPosition(pr1, minTagLen);
	computeTwoTagsPosition(pr2, minTagLen);
	
	if (flag & CPU_RAW) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
	} else if (flag & CPU_FASTQ) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		uint16_t totalTags = pr1->tags_n + pr2->tags_n;
		buildTwoTagsFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildTwoTagsFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}

static inline void extractHalfHalfTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeTwoTagsPosition(pr1, minTagLen);
	computeTwoTagsPosition(pr2, minTagLen);
	
	if (flag & CPU_RAW) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
	} else if (flag & CPU_FASTQ) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		uint16_t totalTags = pr1->tags_n + pr2->tags_n;
		buildTwoTagsFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildTwoTagsFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}

static inline void extractHalfNoneTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeTwoTagsPosition(pr1, minTagLen);
	computeSingleTagPosition(pr2, minTagLen);
	
	if (flag & CPU_RAW) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
	} else if (flag & CPU_FASTQ) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		uint16_t totalTags = pr1->tags_n + pr2->tags_n;
		buildTwoTagsFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildSingleTagFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}

static inline void extractNoneHalfTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeSingleTagPosition(pr1, minTagLen);
	computeTwoTagsPosition(pr2, minTagLen);
	
	if (flag & CPU_RAW) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
	} else if (flag & CPU_FASTQ) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		uint16_t totalTags = pr1->tags_n + pr2->tags_n;
		buildSingleTagFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildTwoTagsFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}

static inline void extractNoneFullTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeSingleTagPosition(pr1, minTagLen);
	computeTwoTagsPosition(pr2, minTagLen);
	
	if (flag & CPU_RAW) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
	} else if (flag & CPU_FASTQ) {
		pfq1->pairId = pfq2->pairId = pid;
		pfq1->readId = 1; pfq2->readId = 2;
		uint16_t totalTags = pr1->tags_n + pr2->tags_n;
		buildSingleTagFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildTwoTagsFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}
*/

// Single-end read version (only R/1)
static inline void processSingleTag
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, int minTagLen, int flag, int tagFamily)
{
	if (NULL!=pSeq1) pr1->readLen = pSeq1->l_seq;
	
	if (LINKER_NONE!=pr1->tied_linker_type) {
		// DECISION: write R/1 into .tied.fastq.gz
		computeSingleTagPosition(pr1, minTagLen);
		
		pr1->finalLinkerType = LINKER_TIED;
		pr1->classification = CPU_OUTPUT_TIED;
		
		if (flag & CPU_RAW) {
			pfq1->pairId = pid;
			pfq1->readId = 1;
			buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
		} else if (flag & CPU_FASTQ) {
			pfq1->pairId = pid;
			pfq1->readId = 1;
			buildSingleTagFq(pr1, pSeq1, pfq1, pr1->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
		}
	} else {
		// at this point, R/1 does not have tied linker type detected
		if (isFullLinker(pr1->linker_type)) {
			// FL(R1) && FL(R2)
			// DECISION: write R/1 and R/2 into .FL.fastq.gz if consistent
			// DECISION: write R/1 and R/2 into .conflict.fastq.gz if inconsistent
			//extractFullFullTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
			computeTwoTagsPosition(pr1, minTagLen);
			
			pr1->finalLinkerType = pr1->linker_type; // there is only a single read, thus will be the same linker
			// OUTPUT:
			//   Tag1(R/1), FL, Tag2(R/1)
			if (isNonChimericFullLinker(pr1->finalLinkerType)) {
				pr1->classification = (isBothTagsPresentR1(pr1, minTagLen) ? CPU_OUTPUT_FL_NC_2TAGS : CPU_OUTPUT_FL_NC_1TAG);
			} else {
				pr1->classification = (isBothTagsPresentR1(pr1, minTagLen) ? CPU_OUTPUT_FL_C_2TAGS : CPU_OUTPUT_FL_C_1TAG);
			}
			
			if (flag & CPU_RAW) {
				pfq1->pairId = pid;
				pfq1->readId = 1;
				buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
			} else if (flag & CPU_FASTQ) {
				pfq1->pairId = pid;
				pfq1->readId = 1;
				buildTwoTagsFq(pr1, pSeq1, pfq1, pr1->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
			}
			
		} else if (isHalfLinker(pr1->linker_type)) {
			// HL(R1)
			// DECISION: write R/1 into .HL.fastq.gz
			
			// OUTPUT:
			//   Tag1(R/1), HL, Tag2(R/1)
			computeTwoTagsPosition(pr1, minTagLen);
			
			pr1->finalLinkerType = pr1->linker_type;
			pr1->classification = CPU_OUTPUT_HL;
			
			if (flag & CPU_RAW) {
				pfq1->pairId = pid;
				pfq1->readId = 1;
				buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
			} else if (flag & CPU_FASTQ) {
				pfq1->pairId = pid;
				pfq1->readId = 1;
				buildTwoTagsFq(pr1, pSeq1, pfq1, pr1->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
			}
		} else {
			// None(R1)
			// DECISION: write R/1 into .none.fastq.gz
			// nDecision = CPU_OUTPUT_NONE;
			// OUTPUT:
			//   R/1
			computeSingleTagPosition(pr1, minTagLen);
			
			pr1->finalLinkerType = LINKER_NONE;
			pr1->classification = CPU_OUTPUT_NONE;
			
			if (flag & CPU_RAW) {
				pfq1->pairId = pid;
				pfq1->readId = 1;
				buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
			} else if (flag & CPU_FASTQ) {
				pfq1->pairId = pid;
				pfq1->readId = 1;
				buildSingleTagFq(pr1, pSeq1, pfq1, pr1->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
			}
		}
	}
	if (flag & (CPU_RAW|CPU_FASTQ)) {
		pfq1->classification = pr1->classification;
	}
}

// Paired-end read version (both R/1 & R/2)
static inline void processPairedTag
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag, int tagFamily)
{
	if (NULL!=pSeq1) pr1->readLen = pSeq1->l_seq;
	if (NULL!=pSeq2) pr2->readLen = pSeq2->l_seq;
	
	if (LINKER_NONE!=pr1->tied_linker_type || LINKER_NONE!=pr2->tied_linker_type) {
		// DECISION: write R/1 and R/2 into .tied.fastq.gz
		// OUTPUT:
		//   Tag1(R/1), <HL,HL>, Tag2(R/1)
		//   Tag1(R/1), <FL,FL>, Tag2(R/1)
		//   Tag1(R/2), <HL,HL>, Tag2(R/2)
		//   Tag1(R/2), <FL,FL>, Tag2(R/2)
		
		//extractNoneNoneTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
		computeSingleTagPosition(pr1, minTagLen);
		computeSingleTagPosition(pr2, minTagLen);
		
		pr1->finalLinkerType = pr2->finalLinkerType = LINKER_TIED;
		pr1->classification = pr2->classification = CPU_OUTPUT_TIED;

		if (flag & CPU_RAW) {
			pfq1->pairId = pfq2->pairId = pid;
			pfq1->readId = 1; pfq2->readId = 2;
			buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
			buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
		} else if (flag & CPU_FASTQ) {
			pfq1->pairId = pfq2->pairId = pid;
			pfq1->readId = 1; pfq2->readId = 2;
			buildSingleTagFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
			buildSingleTagFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
		}
	} else {
		// at this point, both R/1 and R/2 do not have tied linker type detected
		
		if (isFullLinker(pr1->linker_type)) {
			if (isFullLinker(pr2->linker_type)) {
				// FL(R1) && FL(R2)
				// DECISION: write R/1 and R/2 into .FL.fastq.gz if consistent
				// DECISION: write R/1 and R/2 into .conflict.fastq.gz if inconsistent
				//extractFullFullTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeTwoTagsPosition(pr1, minTagLen);
				computeTwoTagsPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = isFLFLConsistent(pr1->linker_type, pr2->linker_type);
				if (LINKER_NONE!=pr1->finalLinkerType) {
					// OUTPUT:
					//   Tag1(R/1), FL, Tag2(R/1)
					//   Tag1(R/2), FL, Tag2(R/2)
					//
					// TODO:
					// BETTER: unsure of the pairing, so keep it to first-class evident
					//   Tag1(R/1) -=- Tag2(R/1)
					//   Tag1(R/2) -=- Tag2(R/2)
					// BEST:
					//       [max(Tag1(R/1), rc(Tag2(R/2)))] -=- FL -=- [max(rc(Tag1(R/2)), Tag2(R/1))]
					//
					//           Tag1(R/1),  FL    , Tag2(R/1)
					//       rc(Tag2(R/2)), rc(FL) , rc(Tag1(R/2))
					//   we expect that there will be overlap between this paired-read
					//
					if (isNonChimericFullLinker(pr1->finalLinkerType)) {
						pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_NC_2TAGS : CPU_OUTPUT_FL_NC_1TAG);
					} else {
						pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_C_2TAGS : CPU_OUTPUT_FL_C_1TAG);
					}
				} else {
					// OUTPUT:
					//   Tag1(R/1), FL, Tag2(R/1)
					//   Tag1(R/2), FL, Tag2(R/2)
					//
					// TODO:
					// BETTER: unsure of the pairing, so keep it to first-class evident
					//   Tag1(R/1) -=- Tag2(R/1)
					//   Tag1(R/2) -=- Tag2(R/2)
					pr1->classification = pr2->classification = CPU_OUTPUT_CONFLICT;
				}
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildTwoTagsFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildTwoTagsFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			} else if (isHalfLinker(pr2->linker_type)) {
				// FL(R1) && HL(R2)
				// DECISION: write R/1 and R/2 into .FL.fastq.gz if consistent
				// DECISION: write R/1 and R/2 into .conflict.fastq.gz if inconsistent
				//extractFullHalfTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeTwoTagsPosition(pr1, minTagLen);
				computeTwoTagsPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = isFLHLConsistent(pr1->linker_type, pr2->linker_type);
				if (LINKER_NONE!=pr1->finalLinkerType) {
					// OUTPUT:
					//   Tag1(R/1), FL, Tag2(R/1)
					//   Tag1(R/2), HL, Tag2(R/2)
					//
					// TODO:
					// BETTER: unsure of the pairing, so keep it to first-class evident
					//   Tag1(R/1) -=- Tag2(R/1)
					//   Tag1(R/2) -=- Tag2(R/2)
					// BEST:
					//       [Tag1(R/1), rc(Tag2(R/2))] -=- FL, rc(HL) -=- [rc(Tag1(R/2)), Tag2(R/1)]
					//
					//       Tag1(R/1), FL        , Tag2(R/1)
					//       rc(Tag2(R/2)), rc(HL), rc(Tag1(R/2))
					//   we expect that there will be overlap between this paired-read
					//
					if (isNonChimericFullLinker(pr1->finalLinkerType)) {
						pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_NC_2TAGS : CPU_OUTPUT_FL_NC_1TAG);
					} else {
						pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_C_2TAGS : CPU_OUTPUT_FL_C_1TAG);
					}
				} else {
					// OUTPUT:
					//   Tag1(R/1), FL, Tag2(R/1)
					//   Tag1(R/2), HL, Tag2(R/2)
					//
					// TODO:
					// BETTER: unsure of the pairing, so keep it to first-class evident
					//   Tag1(R/1) -=- Tag2(R/1)
					//   Tag1(R/2) -=- Tag2(R/2)
					pr1->classification = pr2->classification = CPU_OUTPUT_CONFLICT;
				}
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildTwoTagsFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildTwoTagsFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			} else {
				// FL(R1) && None(R2)
				// DECISION: write R/1 and R/2 into .FL.fastq.gz
				
				// OUTPUT:
				//   Tag1(R/1), FL, Tag2(R/1) - possibly overlap with RC(R/2) or with gap
				//   TODO: function to check if it is overlap or with gap
				//   R/2
				//extractFullNoneTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeTwoTagsPosition(pr1, minTagLen);
				computeSingleTagPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = pr1->linker_type;
				if (isNonChimericFullLinker(pr1->finalLinkerType)) {
					pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_NC_2TAGS : CPU_OUTPUT_FL_NC_1TAG);
				} else {
					pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_C_2TAGS : CPU_OUTPUT_FL_C_1TAG);
				}
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildTwoTagsFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildSingleTagFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			}
		} else if (isHalfLinker(pr1->linker_type)) {
			if (isFullLinker(pr2->linker_type)) {
				// HL(R1) && FL(R2)
				// DECISION: write R/1 and R/2 into .FL.fastq.gz if consistent
				// DECISION: write R/1 and R/2 into .conflict.fastq.gz if inconsistent
				//extractHalfFullTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeTwoTagsPosition(pr1, minTagLen);
				computeTwoTagsPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = isHLFLConsistent(pr1->linker_type, pr2->linker_type);
				if (LINKER_NONE!=pr1->finalLinkerType) {
					// OUTPUT:
					//   Tag1(R/1), HL, Tag2(R/1)
					//   Tag1(R/2), FL, Tag2(R/2)
					//
					// TODO:
					// BETTER: unsure of the pairing, so keep it to first-class evident
					//   Tag1(R/1) -=- Tag2(R/1)
					//   Tag1(R/2) -=- Tag2(R/2)
					// BEST:
					//       [Tag1(R/1), rc(Tag2(R/2))] -=- HL, rc(FL) -=- [rc(Tag1(R/2)), Tag2(R/1)]
					//
					//       Tag1(R/1), HL, Tag2(R/1)
					//   rc(Tag2(R/2)), rc(FL), rc(Tag1(R/2))
					//   we expect that there will be overlap between this paired-read
					//
					if (isNonChimericFullLinker(pr1->finalLinkerType)) {
						pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_NC_2TAGS : CPU_OUTPUT_FL_NC_1TAG);
					} else {
						pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_C_2TAGS : CPU_OUTPUT_FL_C_1TAG);
					}
				} else {
					// OUTPUT:
					//   Tag1(R/1), HL, Tag2(R/1)
					//   Tag1(R/2), FL, Tag2(R/2)
					//
					// TODO:
					// BETTER: unsure of the pairing, so keep it to first-class evident
					//   Tag1(R/1) -=- Tag2(R/1)
					//   Tag1(R/2) -=- Tag2(R/2)
					pr1->classification = pr2->classification = CPU_OUTPUT_CONFLICT;
				}
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildTwoTagsFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildTwoTagsFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			} else if (isHalfLinker(pr2->linker_type)) {
				// HL(R1) && HL(R2)
				// DECISION: write R/1 and R/2 into .FL.fastq.gz if consistent
				// DECISION: write R/1 and R/2 into .conflict.fastq.gz if inconsistent
				//extractHalfHalfTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeTwoTagsPosition(pr1, minTagLen);
				computeTwoTagsPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = isHLHLConsistent(pr1->linker_type, pr2->linker_type);
				if (LINKER_NONE!=pr1->finalLinkerType) {
					// OUTPUT:
					//   Tag1(R/1), HL, Tag2(R/1)
					//   Tag1(R/2), HL, Tag2(R/2)
					// TODO:
					//       Tag1(R/1), HL    , Tag2(R/1)
					//           rc(Tag2(R/2)), rc(HL), rc(Tag1(R/2))
					//   we expect that there will be overlap between this paired-read
					// BETTER:
					//       Tag1(R/1) -=- HL, rc(HL) -=- rc(Tag1(R/2))
					//
					if (isNonChimericFullLinker(pr1->finalLinkerType)) {
						pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_NC_2TAGS : CPU_OUTPUT_FL_NC_1TAG);
					} else {
						pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_C_2TAGS : CPU_OUTPUT_FL_C_1TAG);
					}
				} else {
					// OUTPUT:
					//   Tag1(R/1), HL, Tag2(R/1)
					//   Tag1(R/2), HL, Tag2(R/2)
					// BETTER:
					//   Tag1(R/1) -=- Tag2(R/1)
					//   Tag1(R/2) -=- Tag2(R/2)
					pr1->classification = pr2->classification = CPU_OUTPUT_CONFLICT;
				}
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildTwoTagsFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildTwoTagsFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			} else {
				// HL(R1) && None(R2)
				// DECISION: write R/1 and R/2 into .HL.fastq.gz
				
				// OUTPUT:
				//   Tag1(R/1), HL, Tag2(R/1) - possibly overlap with RC(R/2) or with gap
				//   TODO: function to check if it is overlap or with gap
				//   R/2
				//extractHalfNoneTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeTwoTagsPosition(pr1, minTagLen);
				computeSingleTagPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = pr1->linker_type;
				pr1->classification = pr2->classification = CPU_OUTPUT_HL;
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildTwoTagsFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildSingleTagFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			}
		} else {
			if (isFullLinker(pr2->linker_type)) {
				// None(R1) && FL(R2)
				// DECISION: write R/1 and R/2 into .FL.fastq.gz
				
				// OUTPUT:
				//   R/1
				//   Tag1(R/2), FL, Tag2(R/2) - possibly overlap with RC(R/1) or with gap
				//   TODO: function to check if it is overlap or with gap
				//extractNoneFullTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeSingleTagPosition(pr1, minTagLen);
				computeTwoTagsPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = pr2->linker_type;
				if (isNonChimericFullLinker(pr2->linker_type)) {
					pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_NC_2TAGS : CPU_OUTPUT_FL_NC_1TAG);
				} else {
					pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_FL_C_2TAGS : CPU_OUTPUT_FL_C_1TAG);
				}
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildSingleTagFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildTwoTagsFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			} else if (isHalfLinker(pr2->linker_type)) {
				// None(R1) && HL(R2)
				// DECISION: write R/1 and R/2 into .HL.fastq.gz
				
				// OUTPUT:
				//   R/1
				//   Tag1(R/2), HL, Tag2(R/2) - possibly overlap with RC(R/1) or with gap
				//   TODO: function to check if it is overlap or with gap
				//extractNoneHalfTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeSingleTagPosition(pr1, minTagLen);
				computeTwoTagsPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = pr2->linker_type;
				pr1->classification = pr2->classification = CPU_OUTPUT_HL;
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildSingleTagFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildTwoTagsFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			} else {
				// None(R1) && None(R2)
				// DECISION: write R/1 and R/2 into .none.fastq.gz
				// nDecision = CPU_OUTPUT_NONE;
				
				// OUTPUT:
				//   R/1
				//   R/2
				//extractNoneNoneTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeSingleTagPosition(pr1, minTagLen);
				computeSingleTagPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = LINKER_NONE;
				pr1->classification = pr2->classification = CPU_OUTPUT_NONE;
				
				if (flag & CPU_RAW) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildRawFq(pr1, pSeq1, pfq1, flag & CPU_RC_READ);
					buildRawFq(pr2, pSeq2, pfq2, flag & CPU_RC_READ);
				} else if (flag & CPU_FASTQ) {
					pfq1->pairId = pfq2->pairId = pid;
					pfq1->readId = 1; pfq2->readId = 2;
					buildSingleTagFq(pr1, pSeq1, pfq1, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
					buildSingleTagFq(pr2, pSeq2, pfq2, pr1->tags_n + pr2->tags_n, minTagLen, flag & CPU_RC_READ, tagFamily);
				}
			}
		}
		
	}
	if (flag & (CPU_RAW|CPU_FASTQ)) {
		pfq1->classification = pfq2->classification = pr1->classification;
	}
}

static inline int init_CPTags_Ouputs (int flag, const char *outPrefix, gzFile *outfds) {
	const char *mode = (flag & CPU_COMPRESS) ? "w6" : "wT"; //"w0";
	const char *outSuffix = (flag & CPU_COMPRESS) ? ".gz" : "";
	
	// CPU_OUTPUT_NONE
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.none.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_NONE] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_NONE] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_NONE+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_NONE], G_GZIP_BUFFER_SIZE);
	}
	
	// CPU_OUTPUT_TIED
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.tied.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_TIED] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_TIED] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_TIED+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_TIED], G_GZIP_BUFFER_SIZE);
	}
	
	// CPU_OUTPUT_CONFLICT
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.conflict.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_CONFLICT] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_CONFLICT] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_CONFLICT+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_CONFLICT], G_GZIP_BUFFER_SIZE);
	}
	
	// CPU_OUTPUT_HL
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.halflinker.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_HL] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_HL] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_HL+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_HL], G_GZIP_BUFFER_SIZE);
	}
	
	// CPU_OUTPUT_FL_C_2TAGS
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.fulllinker.chimeric.paired.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_FL_C_2TAGS] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_FL_C_2TAGS] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_FL_C_2TAGS+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_FL_C_2TAGS], G_GZIP_BUFFER_SIZE);
	}
	
	// CPU_OUTPUT_FL_C_1TAG
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.fulllinker.chimeric.single.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_FL_C_1TAG] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_FL_C_1TAG] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_FL_C_1TAG+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_FL_C_1TAG], G_GZIP_BUFFER_SIZE);
	}
	
	// CPU_OUTPUT_FL_NC_2TAGS
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.FullLinker.NonChimeric.paired.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_FL_NC_2TAGS] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_FL_NC_2TAGS] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_FL_NC_2TAGS+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_FL_NC_2TAGS], G_GZIP_BUFFER_SIZE);
	}
	
	// CPU_OUTPUT_FL_NC_1TAG
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.FullLinker.NonChimeric.single.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_FL_NC_1TAG] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_FL_NC_1TAG] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_FL_NC_1TAG+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_FL_NC_1TAG], G_GZIP_BUFFER_SIZE);
	}
	
	return 0;
}

static inline int terminate_CPTags_Ouputs (gzFile *outfds) {
	// TODO: always synchronized with #defined
	int nFailed = 0;
	if (outfds[CPU_OUTPUT_NONE]) {
		if  (gzclose(outfds[CPU_OUTPUT_NONE])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-no-linker gzip stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_TIED]) {
		if (gzclose(outfds[CPU_OUTPUT_TIED])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-tied stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_CONFLICT]) {
		if (gzclose(outfds[CPU_OUTPUT_CONFLICT])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-conflict stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_HL]) {
		if (gzclose(outfds[CPU_OUTPUT_HL])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-half-linker stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_FL_C_2TAGS]) {
		if (gzclose(outfds[CPU_OUTPUT_FL_C_2TAGS])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-full-linker-chimeric-paired stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_FL_C_1TAG]) {
		if (gzclose(outfds[CPU_OUTPUT_FL_C_1TAG])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-full-linker-chimeric-single stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_FL_NC_2TAGS]) {
		if (gzclose(outfds[CPU_OUTPUT_FL_NC_2TAGS])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-full-linker-non-chimeric-paired stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_FL_NC_1TAG]) {
		if (gzclose(outfds[CPU_OUTPUT_FL_NC_1TAG])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-full-linker-non-chimeric-single stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	
	return nFailed;
}

#endif
