#ifndef CPU_TAG_COMMON_H
#define CPU_TAG_COMMON_H

#include <errno.h>

#include "kstring.h"

#include "zlib.h"

// available in >=v1.2.4
#define G_GZIP_BUFFER_SIZE (128*1024)

#define CPU_RECORD_VERSION "CPU-0.1.0a"

#define CPU_OUTPUT_PREFIX     "cptags"

#define CPU_LABEL      0x0001
#define CPU_NAME       0x0002
#define CPU_SEQ        0x0004
#define CPU_QUALITY    0x0008
#define CPU_DEBUG_MASK (CPU_LABEL|CPU_NAME|CPU_SEQ|CPU_QUALITY)

#define CPU_COMPRESS   0x0010
#define CPU_FASTQ      0x0020
#define CPU_RAW        0x0040
#define CPU_RC_READ    0x1000


#define CPU_OUTPUT_NONE          0x0000
#define CPU_OUTPUT_TIED          0x0001
#define CPU_OUTPUT_CONFLICT      0x0002
#define CPU_OUTPUT_COMMON_COUNT  0x0003
//#define CPU_OUTPUT_HL          0x0003
//#define CPU_OUTPUT_FL_C_2TAGS  0x0004
//#define CPU_OUTPUT_FL_C_1TAG   0x0005
//#define CPU_OUTPUT_FL_NC_2TAGS 0x0006
//#define CPU_OUTPUT_FL_NC_1TAG  0x0007
//#define CPU_OUTPUT_COUNT       0x0008


#define CPU_TAGTYPE_DILINKER     (0<<5)
#define CPU_TAGTYPE_LMP          (1<<5)
#define CPU_TAGTYPE_SILINKER     (2<<5)
#define CPU_TAGTYPE_EXTENDED     (7<<5)


#include "cpuadapter.h"
#include "cpulinker.h"

typedef struct {
	long pairId;
	uint8_t readId;
	uint8_t  classification;
	kstring_t fq;
} cpu_fq_t;

typedef struct {
	int16_t adapter_type;
	uint16_t adapter_score1;
	int16_t adapter_ref_begin1;
	int16_t adapter_ref_end1;
	int16_t	adapter_read_begin1;
	int16_t adapter_read_end1;
	
	int16_t linker_type;
	uint16_t linker_score1;
	int16_t linker_ref_begin1;
	int16_t linker_ref_end1;
	int16_t	linker_read_begin1;
	int16_t linker_read_end1;
	int16_t tied_linker_type;
	
	uint16_t finalLinkerType;
	uint8_t  classification;
	uint8_t  tags_n;
	uint16_t readLen;
	uint16_t start[2];
	uint16_t len[2];
} cpu_record_t;

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/* This table is used to transform nucleotide letters into numbers. */
extern int8_t nt_table[128];

static void reverse_comple(const char* seq, char* rc, size_t end) {
	//size_t end = strlen(seq);
	size_t start = 0;
	static const int8_t rc_table[128] = {
		78, 78, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,
		78, 78, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,
		78, 78, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,
		78, 78, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,
		78, 84, 78, 71,  78, 78, 78, 67,  78, 78, 78, 78,  78, 78, 78, 78,
		78, 78, 78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78,
		78, 84, 78, 71,  78, 78, 78, 67,  78, 78, 78, 78,  78, 78, 78, 78,
		78, 78, 78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78
	};
	rc[end] = '\0';
	-- end;
	while (LIKELY(start < end)) {
		rc[start] = (char)rc_table[(int8_t)seq[end]];
		rc[end] = (char)rc_table[(int8_t)seq[start]];
		++ start;
		-- end;
	}
	if (start == end) rc[start] = (char)rc_table[(int8_t)seq[start]];
}

// TODO:
// 1) need quality trimming to be included
// 2) need to record how many tags survived for mapping (inject into the id?! better for mapping results processing later)
// 3) currently keep all reads which output still have to check. can be better by keeping array of tag(s) to write
// 4) check if we have linker at the 5' end of a read?

// pairing version of tag extraction

static inline int isBothTagsPresentR1(cpu_record_t *pr1, int minTagLen) {
	int nStatus = 0;
	if (pr1->len[0]>=minTagLen) nStatus |= 0x02;
	if (pr1->len[1]>=minTagLen) nStatus |= 0x01;

	return (0x03==nStatus);
}

static inline int isBothTagsPresent(cpu_record_t *pr1, cpu_record_t *pr2, int minTagLen) {
#if 0
	if (pr1->len[0]>=minTagLen && pr2->len[0]>=minTagLen)
		return 1;
	if (pr1->len[1]>=minTagLen && pr2->len[1]>=minTagLen)
		return 2;
	return 0;
#else
	// WCH: 2015-03-10 bugfix
	int nStatus = 0;
	if (pr1->len[0]>=minTagLen) nStatus |= 0x02;
	if (pr1->len[1]>=minTagLen) nStatus |= 0x01;
	if (pr2->len[0]>=minTagLen) nStatus |= 0x01;
	if (pr2->len[1]>=minTagLen) nStatus |= 0x02;
	
	return (0x03==nStatus);
#endif
}

static inline void computeSingleTagPosition (cpu_record_t *pr, int minTagLen)
{
	pr->start[0] = 0;
	int nLen = (ADAPTER_NONE==pr->adapter_type) ? pr->readLen : (pr->adapter_read_begin1);
	pr->len[0] = nLen;
	if (pr->len[0]>=minTagLen) pr->tags_n++;
	
	pr->start[1] = 0;
	pr->len[1] = 0;
}

static inline void computeTwoTagsPosition (cpu_record_t *pr, int minTagLen)
{
	pr->tags_n = 0;
	
	pr->start[0] = 0;
	pr->len[0] = pr->linker_read_begin1;
	if (pr->len[0]>=minTagLen) pr->tags_n++;
	
	pr->start[1] = pr->linker_read_end1 + 1;
	int nEnd = (ADAPTER_NONE==pr->adapter_type) ? pr->readLen : (pr->adapter_read_begin1);
	int nLen = nEnd - pr->start[1];
	pr->len[1] = (nLen<0) ? 0 : nLen;
	if (pr->len[1]>=minTagLen) pr->tags_n++;
}

static inline void strrev(char *s)
{
	if (!s || !(*s)) return;
	char *e = s + strlen(s) - 1;
	while( e > s ) {
		char t = *s;
		*s = *e;
		*e = t;
		s++;
		e--;
	}
}

/*
 // this is the original verbose read id which is missing the tile, x, y information for de-duplication
static inline void buildReadId (const cpu_record_t *pr, const cpu_fq_t *pfq, uint16_t totalTags, uint16_t numTags, uint16_t tagId, const char *position, char *szBuffer)
{
	sprintf(szBuffer, "%d:%ld:%d:%d:%d:%d:%d:%d:%d", totalTags, pfq->pairId, pfq->readId, numTags, tagId, pr->linker_type, pr->tied_linker_type, pr->finalLinkerType, pr->classification);
}
*/

static inline const char *getLanePosition (const char *name)
{
	size_t len = strlen(name);
	const char *position = name + len;
	uint8_t numColon = 0;
	while (position >= name) {
		if (':' == *position) {
			numColon++;
			if (3==numColon) return position;
		}
		position--;
	}
	return "";
}

typedef struct {
	long pairId;
	//num2
	uint8_t  tagFamily;
	uint8_t  totalTags;
	uint8_t  readId;
	uint8_t  numTags;
	uint8_t  tagId;
	
	//num3
	uint8_t tied_linker_type;
	uint8_t linker_type;
	
	//num4
	uint8_t classification;
	uint8_t finalLinkerType;
	
	//position
	int32_t tile;
	int32_t x;
	int32_t y;
} cpu_readid_t;

static inline int buildPairReadId (const cpu_readid_t *pRead1, const cpu_readid_t *pRead2, char *szQName)
{
	// see void buildReadId(...) for details
	if (pRead1->pairId==pRead2->pairId && pRead1->tile==pRead2->tile && pRead1->x==pRead2->x && pRead1->y==pRead2->y) {
		
#if 0
		uint16_t num2 = (pRead2->tagFamily | (pRead2->totalTags-1)<<3 | ((pRead2->readId-1)<<2) | ((pRead2->numTags-1)<<1) | (pRead2->tagId-1)) << 8;
		num2 |= (pRead1->tagFamily | (pRead1->totalTags-1)<<3 | ((pRead1->readId-1)<<2) | ((pRead1->numTags-1)<<1) | (pRead1->tagId-1));
#else
		uint16_t num2 = (pRead2->tagFamily | (pRead2->totalTags)<<3 | ((pRead2->readId)<<2) | ((pRead2->numTags)<<1) | (pRead2->tagId)) << 8;
		num2 |= (pRead1->tagFamily | (pRead1->totalTags)<<3 | ((pRead1->readId)<<2) | ((pRead1->numTags)<<1) | (pRead1->tagId));
#endif	
		uint32_t num3 = (pRead2->tied_linker_type << 8 | pRead2->linker_type); num3 = num3 << 16;
		num3 |= (pRead1->tied_linker_type << 8 | pRead1->linker_type);
		
		uint32_t num4 = (pRead2->classification << 8 | pRead2->finalLinkerType); num4 = num4 << 16;
		num4 |= (pRead1->classification << 8 | pRead1->finalLinkerType);
		
		return sprintf(szQName, "%ld:%u:%u:%u:%d:%d:%d", pRead1->pairId, num2, num3, num4, pRead1->tile, pRead1->x, pRead1->y);
	} else {
		szQName[0] = '\0';
		return 0;
	}
}

static inline int parseReadId (const char *readId, cpu_readid_t *pRead, int isRead2)
{
	long pairId=0;
	int num2=0;
	int num3=0;
	int num4=0;
	int tile=0;
	int x=0;
	int y=0;
	
	// TODO: this is taking much longer than necessary
	//       implement a faster version
#if 0
	int ret = sscanf(readId, "%ld:%u:%u:%u:%d:%d:%d", &pairId, &num2, &num3, &num4, &tile, &x, &y);
	if (7!=ret) {
		return ret;
	}
#else
	int ret = 0;
	const char *pch = readId;
	while (*pch!='\0' && *pch!=':') {
		if (isdigit(*pch)) { pairId *= 10; pairId += (*pch - '0'); }
		else break;
		pch++;
	}
	if (':' == *pch) {
		ret++; // 1
		
		pch++;
		while (*pch!='\0' && *pch!=':') {
			if (isdigit(*pch)) { num2 *= 10; num2 += (*pch - '0'); }
			else break;
			pch++;
		}
		
		if (':' == *pch) {
			ret++; // 2
			
			pch++;
			while (*pch!='\0' && *pch!=':') {
				if (isdigit(*pch)) { num3 *= 10; num3 += (*pch - '0'); }
				else break;
				pch++;
			}
			
			if (':' == *pch) {
				ret++; // 3
				
				pch++;
				while (*pch!='\0' && *pch!=':') {
					if (isdigit(*pch)) { num4 *= 10; num4 += (*pch - '0'); }
					else break;
					pch++;
				}
				
				if (':' == *pch) {
					ret++; // 4
					
					pch++;
					while (*pch!='\0' && *pch!=':') {
						if (isdigit(*pch)) { tile *= 10; tile += (*pch - '0'); }
						else break;
						pch++;
					}
					
					if (':' == *pch) {
						ret++; // 5
						
						pch++;
						while (*pch!='\0' && *pch!=':') {
							if (isdigit(*pch)) { x *= 10; x += (*pch - '0'); }
							else break;
							pch++;
						}
						
						if (':' == *pch) {
							ret++; // 6
							
							pch++;
							while (*pch!='\0' && *pch!=':') {
								if (isdigit(*pch)) { y *= 10; y += (*pch - '0'); }
								else break;
								pch++;
							}
							
							if ('\0' == *pch) {
								ret++; // 7
							}
						}
					}
				}
			}
		}
	}
	
	if (7!=ret) {
		return ret;
	}
#endif
	
	// pairid
	pRead->pairId = pairId;

	if (0!=isRead2) {
		num2 = num2 >> 8;
		num3 = num3 >> 16;
		num4 = num4 >> 16;
	}
	//num2
	pRead->tagFamily = (num2 & 0xe0);
#if 0
	pRead->totalTags = ((num2>>3) & 0x03) + 1;
	pRead->readId = ((num2>>2) & 0x01) + 1;
	pRead->numTags = ((num2>>1) & 0x01) + 1;
	pRead->tagId = (num2 & 0x01) + 1;
#else
	pRead->totalTags = ((num2>>3) & 0x03);
	pRead->readId = ((num2>>2) & 0x01);
	pRead->numTags = ((num2>>1) & 0x01);
	pRead->tagId = (num2 & 0x01);
#endif
	//num3
	pRead->tied_linker_type = (num3 >> 8)& 0xFF;
	pRead->linker_type = num3 & 0xFF;
	
	//num4
	pRead->classification = (num4 >> 8)& 0xFF;
	pRead->finalLinkerType = num4 & 0xFF;
	
	//position
	pRead->tile = tile;
	pRead->x = x;
	pRead->y = y;
	
	return 0;
}

static inline void buildReadId (const cpu_record_t *pr, const cpu_fq_t *pfq, int tagFamily, uint16_t totalTags, uint16_t numTags, uint16_t tagId, const char *position, char *szBuffer)
{
	//@<#tags-in-pairs>:<pair-id>:<read1/2>:<#Tag>:<TagId>:<linkerType>:<tiedLinkertype>:<finalLinkerType>:<classification>
	// ---> transformed --->
	//@<pair-id>
	//:[5-bits as dec]total(2b),readid(1b),totalTgs(1b),tagid(1b)
	//:[16-bits as dec]tied(8b)linker(8b)
	//:[16-bits as dec]class(8b)final linker(8b)
	//:tile:x:y
	uint8_t num2 = tagFamily | (totalTags-1)<<3 | ((pfq->readId-1)<<2) | ((numTags-1)<<1) | (tagId-1);
	uint16_t num3 = pr->tied_linker_type << 8 | pr->linker_type;
	uint16_t num4 = pr->classification << 8 | pr->finalLinkerType;
	sprintf(szBuffer, "%ld:%d:%d:%d%s", pfq->pairId, num2, num3, num4, position);
}

static inline void buildSingleTagFq (const cpu_record_t *pr, const bseq1_t *pSeq, cpu_fq_t *pfq, uint16_t totalTags, int minTagLen, int rcRead, int tagFamily)
{
	//if (!genFastq) return;
	
	ks_set(NULL, &(pfq->fq));
	
	if (pr->len[0]>=minTagLen) {
		// build the full fq
		
		const char *position = getLanePosition (pSeq->name);
		//@<#tags-in-pairs>:<pair-id>:<read1/2>:<#Tag>:<TagId>:<linkerType>:<tiedLinkertype>:<finalLinkerType>:<classification>
		char szBuffer[1024];
		buildReadId (pr, pfq, tagFamily, totalTags, 1, 1, position, szBuffer);
		//sprintf(szBuffer, "%d:%ld:%d:%d:%d:%d:%d:%d:%d", totalTags, pfq->pairId, pfq->readId, 1, 1, pr->linker_type, pr->tied_linker_type, pr->finalLinkerType, pr->classification);
		
		// id
		kputc('@', &(pfq->fq));
		kputs(szBuffer, &(pfq->fq));
		kputc('\n', &(pfq->fq));
		
		// sequence
		char *seq = malloc(pr->len[0]+1);
		if (!rcRead) {
			strncpy(seq, pSeq->seq, pr->len[0]); seq[pr->len[0]] = '\0';
		} else {
			reverse_comple(pSeq->seq, seq, pr->len[0]);
		}
		kputs(seq, &(pfq->fq)); free(seq);
		kputc('\n', &(pfq->fq));
		
		// separator
		kputs("+\n", &(pfq->fq));
		
		// quality
		char *qual = malloc(pr->len[0]+1);
		strncpy(qual, pSeq->qual, pr->len[0]); qual[pr->len[0]] = '\0';
		if (rcRead) strrev(qual);
		kputs(qual, &(pfq->fq)); free(qual);
		kputc('\n', &(pfq->fq));
	}
}

static inline void buildTwoTagsFq (const cpu_record_t *pr, const bseq1_t *pSeq, cpu_fq_t *pfq, uint16_t totalTags, int minTagLen, int rcRead, int tagFamily)
{
	//if (!genFastq) return;
	
	ks_set(NULL, &(pfq->fq));
	const char *position = NULL;
	
	if (pr->len[0]>=minTagLen) {
		// build the full fq
		
		position = getLanePosition (pSeq->name);
		//@<#tags-in-pairs>:<pair-id>:<read1/2>:<#Tag>:<TagId>:<linkerType>:<tiedLinkertype>:<finalLinkerType>:<classification>
		char szBuffer[1024];
		buildReadId (pr, pfq, tagFamily, totalTags, 2, 1, position, szBuffer);
		//sprintf(szBuffer, "%d:%ld:%d:%d:%d:%d:%d:%d:%d", totalTags, pfq->pairId, pfq->readId, 2, 1, pr->linker_type, pr->tied_linker_type, pr->finalLinkerType, pr->classification);
		
		// id
		kputc('@', &(pfq->fq));
		kputs(szBuffer, &(pfq->fq));
		kputc('\n', &(pfq->fq));
		
		// sequence
		char *seq = malloc(pr->len[0]+1);
		if (!rcRead) {
			strncpy(seq, pSeq->seq, pr->len[0]); seq[pr->len[0]] = '\0';
		} else {
			reverse_comple(pSeq->seq, seq, pr->len[0]);
		}
		kputs(seq, &(pfq->fq)); free(seq);
		kputc('\n', &(pfq->fq));
		
		// separator
		kputs("+\n", &(pfq->fq));
		
		// quality
		char *qual = malloc(pr->len[0]+1);
		strncpy(qual, pSeq->qual, pr->len[0]); qual[pr->len[0]] = '\0';
		if (rcRead) strrev(qual);
		kputs(qual, &(pfq->fq)); free(qual);
		kputc('\n', &(pfq->fq));
	}
	if (pr->len[1]>=minTagLen) {
		// build the full fq
		
		if (!position) position = getLanePosition (pSeq->name);
		//@<#tags-in-pairs>:<pair-id>:<read1/2>:<#Tag>:<TagId>:<linkerType>:<tiedLinkertype>:<finalLinkerType>:<classification>
		char szBuffer[1024];
		buildReadId (pr, pfq, tagFamily, totalTags, 2, 2, position, szBuffer);
		//sprintf(szBuffer, "%d:%ld:%d:%d:%d:%d:%d:%d:%d", totalTags, pfq->pairId, pfq->readId, 2, 2, pr->linker_type, pr->tied_linker_type, pr->finalLinkerType, pr->classification);
		
		// id
		kputc('@', &(pfq->fq));
		kputs(szBuffer, &(pfq->fq));
		kputc('\n', &(pfq->fq));
		
		// sequence
		char *seq = malloc(pr->len[1]+1);
#ifdef AS_IS_R1_R2
		if (!rcRead) {
			strncpy(seq, pSeq->seq+pr->start[1], pr->len[1]); seq[pr->len[1]] = '\0';
		} else {
			reverse_comple(pSeq->seq+pr->start[1], seq, pr->len[1]);
		}
#else
		// we wish to convert R/1 tag2 as R/2 tag1, R/2 tag2 as R/1 tag1
		if (!rcRead) {
			reverse_comple(pSeq->seq+pr->start[1], seq, pr->len[1]);
		} else {
			strncpy(seq, pSeq->seq+pr->start[1], pr->len[1]); seq[pr->len[1]] = '\0';
		}
#endif
		kputs(seq, &(pfq->fq)); free(seq);
		kputc('\n', &(pfq->fq));
		
		// separator
		kputs("+\n", &(pfq->fq));
		
		// quality
		char *qual = malloc(pr->len[1]+1);
		strncpy(qual, pSeq->qual+pr->start[1], pr->len[1]); qual[pr->len[1]] = '\0';
#ifdef AS_IS_R1_R2
		if (rcRead) strrev(qual);
#else
		if (!rcRead) strrev(qual);
#endif
		kputs(qual, &(pfq->fq)); free(qual);
		kputc('\n', &(pfq->fq));
	}
}

static inline void buildRawFq (const cpu_record_t *pr, const bseq1_t *pSeq, cpu_fq_t *pfq, int rcRead)
{
	//if (!genFastq) return;
	
	ks_set(NULL, &(pfq->fq));
	// id
	kputc('@', &(pfq->fq));
	kputs(pSeq->name, &(pfq->fq));
	
	if (pSeq->comment) {
		kputc(' ', &(pfq->fq));
		kputs(pSeq->comment, &(pfq->fq));
	}
		
	// let's record the decision after the name
	char szBuffer[1024];
	// TODO: we do not have these two decisions made
	//sprintf(szBuffer, " %s %s", linkerClassToStr(pr->finalLinkerType), outputClassToStr(pr->classification));
	//kputs(szBuffer, &(pfq->fq));
	
	if (ADAPTER_NONE!=pr->adapter_type) {
		sprintf(szBuffer, " %s(%d,%d,%d)", adapterTypeToStr(pr->adapter_type), pr->adapter_read_begin1, pr->adapter_read_end1, pr->adapter_score1);
		kputs(szBuffer, &(pfq->fq));
	}
	
	if (LINKER_NONE!=pr->linker_type) {
		sprintf(szBuffer, " %s(%d,%d,%d)", linkerTypeToStr(pr->linker_type), pr->linker_read_begin1, pr->linker_read_end1, pr->linker_score1);
		kputs(szBuffer, &(pfq->fq));
		if (LINKER_NONE!=pr->tied_linker_type) {
			sprintf(szBuffer, " tied(%s)", linkerTypeToStr(pr->tied_linker_type));
			kputs(szBuffer, &(pfq->fq));
		}
	} else {
		sprintf(szBuffer, " NoLinker");
		kputs(szBuffer, &(pfq->fq));
	}
	
	kputc('\n', &(pfq->fq));
	
	// sequence
	kputs(pSeq->seq, &(pfq->fq));
	kputc('\n', &(pfq->fq));
	
	// separator
	kputs("+\n", &(pfq->fq));
	
	// quality
	kputs(pSeq->qual, &(pfq->fq));
	kputc('\n', &(pfq->fq));
}

/*
static inline void extractNoneNoneTags
(int64_t pid, cpu_record_t *pr1, const bseq1_t *pSeq1, cpu_fq_t *pfq1, cpu_record_t *pr2, const bseq1_t *pSeq2, cpu_fq_t *pfq2, int minTagLen, int flag)
{
	computeSingleTagPosition(pr1, minTagLen);
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
		buildSingleTagFq(pr1, pSeq1, pfq1, totalTags, minTagLen, flag & CPU_RC_READ);
		buildSingleTagFq(pr2, pSeq2, pfq2, totalTags, minTagLen, flag & CPU_RC_READ);
	}
}
*/

#endif
