#ifndef CPU_STAG_H
#define CPU_STAG_H

#include <errno.h>

#include "kstring.h"

#include "zlib.h"

#include "cpuadapter.h"
#include "cpulinker.h"
#include "cputagcommon.h"

//#define CPU_OUTPUT_NONE        0x0000
//#define CPU_OUTPUT_TIED        0x0001
//#define CPU_OUTPUT_CONFLICT    0x0002
#define CPU_OUTPUT_SL_2TAGS    0x0003
#define CPU_OUTPUT_SL_1TAG     0x0004
#define CPU_OUTPUT_SL_COUNT    0x0005

#define G_CP_LinkerSingle "ACGCGATATCTTATCTGACT"

// assume linkerA and rcLinkerA
#define LINKER_RC_SINGLE_PRIV    0x0001
#define LINKER_SINGLE_PRIV       0x0002

#define LINKER_RC_SINGLE    (LINKER_RC_SINGLE_PRIV|LINKER_FULL)
#define LINKER_SINGLE       (LINKER_SINGLE_PRIV|LINKER_FULL)


static char outputSLClassBuffer[11];
static const char *outputSLClassToStr(int16_t t) {
	if (CPU_OUTPUT_NONE==t) return "none";
	else if (CPU_OUTPUT_TIED==t) return "tied";
	else if (CPU_OUTPUT_CONFLICT==t) return "conflict";
	else if (CPU_OUTPUT_SL_2TAGS==t) return "sl2t";
	else if (CPU_OUTPUT_SL_1TAG==t) return "sl1t";
	else {
		//return "?";
		sprintf(outputSLClassBuffer, "%#x", t);
		return outputSLClassBuffer;
	}
}

static char slinkerTypeBuffer[11];
static const char *slinkerTypeToStr(int16_t t) {
	if (LINKER_NONE==t) return ".";
	else if (LINKER_SINGLE==(LINKER_SINGLE&t)) return "SL";
	else if (LINKER_RC_SINGLE==(LINKER_RC_SINGLE&t)) return "ls";
	else {
		//return "?";
		sprintf(slinkerTypeBuffer, "%#x", t);
		return slinkerTypeBuffer;
	}
}

static inline int isSLSLConsistent(int l1, int l2) {
	if (LINKER_SINGLE==l1) {
		return (LINKER_RC_SINGLE==l2) ? LINKER_SINGLE : LINKER_NONE;
	} else if (LINKER_RC_SINGLE==l1) {
		return (LINKER_SINGLE==l2) ? LINKER_RC_SINGLE : LINKER_NONE;
	}
	return LINKER_NONE;
}

// TODO:
// 1) need quality trimming to be included
// 2) need to record how many tags survived for mapping (inject into the id?! better for mapping results processing later)
// 3) currently keep all reads which output still have to check. can be better by keeping array of tag(s) to write
// 4) check if we have linker at the 5' end of a read?

// pairing version of tag extraction
/*
static inline void extractSLSLTags
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

static inline void extractSLNoneTags
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

static inline void extractNoneSLTags
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

static inline void processPairedSTag
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
				//extractSLSLTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeTwoTagsPosition(pr1, minTagLen);
				computeTwoTagsPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = isSLSLConsistent(pr1->linker_type, pr2->linker_type);
				if (LINKER_NONE!=pr1->finalLinkerType) {
					pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_SL_2TAGS : CPU_OUTPUT_SL_1TAG);
				} else {
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
				//extractSLNoneTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeTwoTagsPosition(pr1, minTagLen);
				computeSingleTagPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = pr1->linker_type;
				pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_SL_2TAGS : CPU_OUTPUT_SL_1TAG);
				
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
				//extractNoneSLTags (pid, pr1, pSeq1, pfq1, pr2, pSeq2, pfq2, minTagLen, flag);
				computeSingleTagPosition(pr1, minTagLen);
				computeTwoTagsPosition(pr2, minTagLen);
				
				pr1->finalLinkerType = pr2->finalLinkerType = pr2->linker_type;
				pr1->classification = pr2->classification = (isBothTagsPresent(pr1, pr2, minTagLen) ? CPU_OUTPUT_SL_2TAGS : CPU_OUTPUT_SL_1TAG);
				
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

static inline int init_CPSTags_Ouputs (int flag, const char *outPrefix, gzFile *outfds) {
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
	
	// CPU_OUTPUT_SL_2TAGS
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.singlelinker.paired.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_SL_2TAGS] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_SL_2TAGS] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_SL_2TAGS+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_SL_2TAGS], G_GZIP_BUFFER_SIZE);
	}
	
	// CPU_OUTPUT_SL_1TAG
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.singlelinker.single.fastq%s", outPrefix, outSuffix);
		outfds[CPU_OUTPUT_SL_1TAG] = gzopen (filename.s, mode);
		if (outfds[CPU_OUTPUT_SL_1TAG] == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for gzip writing: %s\n", __func__, filename.s, strerror (errno));
			free(filename.s);
			return CPU_OUTPUT_SL_1TAG+1;
		}
		free(filename.s);
		gzbuffer(outfds[CPU_OUTPUT_SL_1TAG], G_GZIP_BUFFER_SIZE);
	}
	return 0;
}

static inline int terminate_CPSTags_Ouputs (gzFile *outfds) {
	// TODO: always synchronized with #defined
	int nFailed = 0;
	if (outfds[CPU_OUTPUT_NONE]) {
		if  (gzclose(outfds[CPU_OUTPUT_NONE])) {
			fprintf(stderr, "[E::%s] fail to close CPSTags-no-linker gzip stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_TIED]) {
		if (gzclose(outfds[CPU_OUTPUT_TIED])) {
			fprintf(stderr, "[E::%s] fail to close CPSTags-tied stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_CONFLICT]) {
		if (gzclose(outfds[CPU_OUTPUT_CONFLICT])) {
			fprintf(stderr, "[E::%s] fail to close CPSTags-conflict stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_SL_2TAGS]) {
		if (gzclose(outfds[CPU_OUTPUT_SL_2TAGS])) {
			fprintf(stderr, "[E::%s] fail to close CPSTags-linker-paired stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	if (outfds[CPU_OUTPUT_SL_1TAG]) {
		if (gzclose(outfds[CPU_OUTPUT_SL_1TAG])) {
			fprintf(stderr, "[E::%s] fail to close CPTags-linker-single stream: %s\n", __func__, strerror (errno));
			nFailed++;
		}
	}
	return nFailed;
}

#endif
