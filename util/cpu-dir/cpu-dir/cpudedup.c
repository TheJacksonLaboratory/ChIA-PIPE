// TODO:
// 1. experimental code for de-duplication
// 2. output bedpe data columns in text file for exploration
// 3. perform sorting and marking of duplicated PETs
// 4. re-read the .bam file and either mark the PET as duplicate or omit PET from output
//

#include "zlib.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h> // for multi-threading
#include <errno.h>
#include <math.h> // for floor
#include "utils.h"
#include "kseq.h"
#include "kstring.h"
#include "bwa.h"
KSEQ_DECLARE(gzFile)

#include "sam_header.h"
#include "sam.h"
#include "cputagcommon.h"
#include "cputagbam.h"
#include "cpulinker.h"
#include "cputag.h"

//extern void qsort_mt(void *a, size_t n, size_t es, int (*compar)(const void *, const void *), int nThread);

extern double G_t_real;
extern const char *CPU_version;

// TODO: move this to be sharable from cpupair
#define CPU_MAPSTATE_NA         0
#define CPU_MAPSTATE_Ns         1
#define CPU_MAPSTATE_NOT_MAPPED 2
#define CPU_MAPSTATE_UNIQUE     3
#define CPU_MAPSTATE_REPEAT     4
#define CPU_MAPSTATE_COUNT      (CPU_MAPSTATE_REPEAT+1)

#define CPU_PAIR_OUTPUT_UU            0
#define CPU_PAIR_OUTPUT_UxxU          1
#define CPU_PAIR_OUTPUT_xx            2
#define CPU_PAIR_OUTPUT_nn            3
#define CPU_PAIR_OUTPUT_UU_DISCORDANT 4
#define CPU_PAIR_OUTPUT_COUNT         (CPU_PAIR_OUTPUT_UU_DISCORDANT+1)

#define G_MAX_QNAME_LEN 256
// END - move this to be sharable from cpupair

#define CPU_DEDUP_OUTPUT_REPRESENTATIVE  0
#define CPU_DEDUP_OUTPUT_DUPLICATED      1
#define CPU_DEDUP_OUTPUT_COUNT         (CPU_DEDUP_OUTPUT_DUPLICATED+1)

#define CPU_DEDUP_DEBUG_CONVERT_LIST     0x01
#define CPU_DEDUP_DEBUG_READ_LIST        0x02
#define CPU_DEDUP_DEBUG_SORTED_LIST      0x04

#if 0
// sort by chrom1, start1, end1, chrom2, start2, end2, strand1, strand2, final linker, linker1, linker2, score
#define CPU_DEDUP_METADATA_QUAL(a) (((a)>>8)&0x7F)

#define CPU_DEDUP_METADATA_STRAND(a)               (((a)>>31)&0x01)
#define CPU_DEDUP_METADATA_STRAND_SET(a,b)         ((a)|(((b)&0x01)<<31))
#define CPU_DEDUP_METADATA_MSTRAND(a)              (((a)>>30)&0x01)
#define CPU_DEDUP_METADATA_MSTRAND_SET(a,b)        ((a)|(((b)&0x01)<<30))
#define CPU_DEDUP_METADATA_CLASSIFICATION(a)       (((a)>>27)&0x07)
#define CPU_DEDUP_METADATA_CLASSIFICATION_SET(a,b) ((a)|(((b)&0x07)<<27))
#define CPU_DEDUP_METADATA_FINALLINKER(a)          (((a)>>23)&0x0F)
#define CPU_DEDUP_METADATA_FINALLINKER_SET(a,b)    ((a)|(((b)&0x0F)<<23))
#define CPU_DEDUP_METADATA_LINKER(a)               (((a)>>19)&0x0F)
#define CPU_DEDUP_METADATA_LINKER_SET(a,b)         ((a)|(((b)&0x0F)<<19))
#define CPU_DEDUP_METADATA_MLINKER(a)              (((a)>>15)&0x0F)
#define CPU_DEDUP_METADATA_MLINKER_SET(a,b)        ((a)|(((b)&0x0F)<<15))
#define CPU_DEDUP_METADATA_QUAL(a)                 (((a)>>8)&0x7F)
#define CPU_DEDUP_METADATA_QUAL_SET(a,b)           ((a)|(((b)&0x7F)<<8))
#define CPU_DEDUP_METADATA_UNIQ(a)                 (((a)>>7)&0x01)
#define CPU_DEDUP_METADATA_UNIQ_SET(a,b)           ((a)|(((b)&0x01)<<7))
#define CPU_DEDUP_METADATA_MUNIQ(a)                (((a)>>6)&0x01)
#define CPU_DEDUP_METADATA_MUNIQ_SET(a,b)          ((a)|(((b)&0x01)<<6))
#define CPU_DEDUP_METADATA_RTID(a)                 (((a)>>4)&0x03)
#define CPU_DEDUP_METADATA_RTID_SET(a,b)           ((a)|(((b)&0x03)<<4))
#define CPU_DEDUP_METADATA_MRTID(a)                (((a)>>2)&0x03)
#define CPU_DEDUP_METADATA_MRTID_SET(a,b)          ((a)|(((b)&0x03)<<2))
#define CPU_DEDUP_METADATA_DUPLICATED(a)           (((a)>>1)&0x01)
#define CPU_DEDUP_METADATA_DUPLICATED_SET(a,b)     ((a)|(((b)&0x01)<<1))
#define CPU_DEDUP_METADATA_SWAP(a)                 ((a)&0x01)
#define CPU_DEDUP_METADATA_SWAP_SET(a,b)           ((a)|((b)&0x01))
#endif

typedef struct {
	int32_t tid;
	int32_t pos;
	int32_t endpos;
	int32_t mtid;
	int32_t mpos;
	int32_t mendpos;
	// TODO: we downgraded qual from 8 bits to 7 bits to keep within 32 bits
	//uint32_t metaData;
	uint32_t strand:1, mstrand:1, classification:3, finalLinkerType:4, linkerType:4, mlinkerType:4, qual:7, uniq:1, muniq:1, rtid:2, mrtid:2, duplicated:1, swap:1;
	// might have to change to uint64_t if supporting more than 4 billion PET tags
	uint32_t lineno; // it is the bam "record id"; in future, we might use this for "direct" bam record access
	uint32_t mlineno; // it is the bam "record id"; in future, we might use this for "direct" bam record access
	long pairId; // TODO: for checking purpose only
} cpu_dedup_record_t;

typedef struct {
	int chunk_size;
	int n_threads;
	
	int selfLigation;
	
	int toGroup;
	int disponoParallel; // TODO: we should have bit based flags!
	int outputSortedList;
	int calculateLibraryComplexity; // TODO: we should have bit based flags!
	
	kstring_t outputPrefix;
} dedup_opt_t;

dedup_opt_t *dedup_opt_init()
{
	dedup_opt_t *o;
	o = calloc(1, sizeof(dedup_opt_t));
	
	o->chunk_size = 1000000;
	o->n_threads = 1;
	
	o->selfLigation = G_SELF_LIGATION;
	
	o->toGroup = 0;
	o->disponoParallel = 0;
	o->outputSortedList = 0;
	o->calculateLibraryComplexity = 1;
	
	memset(&(o->outputPrefix), 0, sizeof(kstring_t));
	return o;
}

void dedup_opt_terminate(dedup_opt_t *o)
{
	free(o->outputPrefix.s);
}

typedef struct {
	const dedup_opt_t *opt;
	
	uint32_t n_processed;
	
	bam1_t *bamRecs;
	cpu_readid_t *readids;
	
	int nTagGroups;
	tag_group_t *tagGroups;
	cpu_dedup_record_t *dedupRecords;
	
	uint32_t nUpdateBam;
} worker_t;

static inline int init_CPDedup_Outputs (const char *outPrefix, char *out_mode, bam_header_t *header, samfile_t **outfds, int selfLigation) {
	// CPU_SPAN_OUTPUT_SELF_LIGATION
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.nr.bam", outPrefix, selfLigation);
		outfds[CPU_DEDUP_OUTPUT_REPRESENTATIVE] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_DEDUP_OUTPUT_REPRESENTATIVE]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_DEDUP_OUTPUT_REPRESENTATIVE+1;
		}
		free(filename.s);
	}
	// CPU_SPAN_OUTPUT_INTRA_LIGATION
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.dup.bam", outPrefix, selfLigation);
		outfds[CPU_DEDUP_OUTPUT_DUPLICATED] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_DEDUP_OUTPUT_DUPLICATED]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_DEDUP_OUTPUT_DUPLICATED+1;
		}
		free(filename.s);
	}
	return 0;
}

static inline int terminate_CPDedup_Outputs (samfile_t **outfds) {
	// TODO: always synchronized with #defined
	int nFailed = 0;
	if (outfds[CPU_DEDUP_OUTPUT_REPRESENTATIVE]) {
		samclose(outfds[CPU_DEDUP_OUTPUT_REPRESENTATIVE]);
	}
	if (outfds[CPU_DEDUP_OUTPUT_DUPLICATED]) {
		samclose(outfds[CPU_DEDUP_OUTPUT_DUPLICATED]);
	}
	
	return nFailed;
}

static void readid_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	// read id parsing does not care about if reads are paried-end
	int nRet = parseReadId(bam1_qname(&(w->bamRecs[i])), &(w->readids[i]), (BAM_FREAD2==(BAM_FREAD2 & w->bamRecs[i].core.flag))?1:0);
	if (0!=nRet) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to parse read id \"%s\". Make sure that read is prepared by CPU.\n", __func__, bam1_qname(&(w->bamRecs[i])));
	}
}

static void dedup_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	// TODO: set up the result in w->dedupRecords[i]
	cpu_dedup_record_t *dedupRec = &(w->dedupRecords[i]);
	// we considered this PET as duplicated by default and non-uniquely mapped
	// this has the same effect as excluding this PET from all downstream analysis
	// thus, we only 'reset' those which are usable for downstream analysis
	dedupRec->duplicated = 0; dedupRec->uniq = 0; dedupRec->muniq = 0;
	
	// we are going to assume that these are paired-end data
	// the data must be name sorted, NOT coorindate-sorted
	// we will also check that the tag is set by CPU
	
	// DECISION: this section is copied and modified from cpupair.c
	//           we likely have to fold both into a single process for code maintenance and efficiency
	//           in the end, we need the pairing and the genomics end position
	
	// TODO: consideration
	//       code UU : it wil be just two entries; #entries=[2,2]
	//       code Ux : [Ux] non-unique in that R/2 anchor has multiple mapping; #entries=[2,inf)
	//       code Ux : [Un] unique, but can have two entries from both R/1 and R/2; #entries=[1,2]
	//       code xU : [xU] non-unique in that R/2 anchor has multiple mapping; #entries=[2,inf)
	//       code xU : [nU] unique, but can have two entries from both R/1 and R/2; #entries=[1,2]
	//       code xx : both anchors have multiple mapping; #entries=[4,inf)
	//       code nn : only two entries to be written, considered as unmapped; #entries=[2,2]
	
	// perform pairing
	w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_nn;
	
	int32_t readTags[2][2] = {{-1,-1},{-1,-1}};
	int32_t readMappedLen[2][2]={{0,0},{0,0}};
	int32_t readEnd[2][2]={{0,0},{0,0}};
	int j, k;
	// STAGE 1: we set up the representative tags in reads for pairing selection
	uint8_t readTagsMapState[2][2] = {{CPU_MAPSTATE_NA,CPU_MAPSTATE_NA},{CPU_MAPSTATE_NA,CPU_MAPSTATE_NA}};
	uint8_t numTagsSet = 0;
	for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
		if (BAM_FSECONDARY==(w->bamRecs[k].core.flag & BAM_FSECONDARY)) {
			// TODO: we ignore sceondary alignment for now
		} else {
			if (-1==readTags[w->readids[k].readId][w->readids[k].tagId]) {
				readTags[w->readids[k].readId][w->readids[k].tagId] = k;
				bam1_t *bam = &(w->bamRecs[k]);
				if (BAM_FUNMAP==(bam->core.flag & BAM_FUNMAP)) {
					readTagsMapState[w->readids[k].readId][w->readids[k].tagId] = CPU_MAPSTATE_NOT_MAPPED;
				} else {
					readEnd[w->readids[k].readId][w->readids[k].tagId] = bam_calend(&bam->core, bam1_cigar(bam));
					readMappedLen[w->readids[k].readId][w->readids[k].tagId] = readEnd[w->readids[k].readId][w->readids[k].tagId] - bam->core.pos;
					uint8_t *xt = bam_aux_get(bam, "XT");
					if (0!=xt) {
						char a_xt = bam_aux2A(xt);
						if ('U' == a_xt) readTagsMapState[w->readids[k].readId][w->readids[k].tagId] = CPU_MAPSTATE_UNIQUE;
						else if ('R' == a_xt) readTagsMapState[w->readids[k].readId][w->readids[k].tagId] = CPU_MAPSTATE_REPEAT;
						else if ('N' == a_xt) readTagsMapState[w->readids[k].readId][w->readids[k].tagId] = CPU_MAPSTATE_Ns;
						else {
							// TODO: this is considered fatal, cannot continue
							//if (bwa_verbose >= 1)
							fprintf(stderr, "[E::%s] Unknown XT tag \"%c\" for \"%s\". Make sure that mapper produce XT tag like bwa aln..\n", __func__, a_xt, bam1_qname(bam));
						}
					} else {
						// TODO: this is considered fatal, cannot continue
						//if (bwa_verbose >= 1)
						fprintf(stderr, "[E::%s] XT tag is not recorded for \"%s\". Make sure that mapper produce XT tag like bwa aln..\n", __func__, bam1_qname(bam));
					}
				}
				
				numTagsSet++;
				if (4==numTagsSet) break; // we have all 4 tags set, proceed with next stage
			}
		}
	}
	
	// STAGE 2: encode the mapping
	
	// STAGE 2.1 : proceed anchor#1
	uint8_t anchorDiscordant = 0;
	uint8_t anchor1MapState = CPU_MAPSTATE_NA;
	int anchor1ReadIndex = -1;
	if (CPU_MAPSTATE_UNIQUE==readTagsMapState[0][0]) {
		if (CPU_MAPSTATE_UNIQUE==readTagsMapState[1][1]) {
			// both unique
			if (w->bamRecs[readTags[0][0]].core.tid == w->bamRecs[readTags[1][1]].core.tid) {
				// TODO: we pick the 5'-most read for now as it should has more reliable bases
				//       in future we might consider the length, etc
				if (readMappedLen[0][0]>=readMappedLen[1][1]) {
					anchor1ReadIndex = readTags[0][0]; dedupRec->endpos = readEnd[0][0];
				} else {
					anchor1ReadIndex = readTags[1][1]; dedupRec->endpos = readEnd[1][1];
				}
			} else {
				// different chromosome! thus consider discordant
				// TODO: but we still pick one for processing
				anchorDiscordant |= 0x01;
				if (readMappedLen[0][0]>=readMappedLen[1][1]) {
					anchor1ReadIndex = readTags[0][0]; dedupRec->endpos = readEnd[0][0];
				} else {
					anchor1ReadIndex = readTags[1][1]; dedupRec->endpos = readEnd[1][1];
				}
			}
		} else {
			// only R/1 tag 1 unique
			anchor1ReadIndex = readTags[0][0]; dedupRec->endpos = readEnd[0][0];
		}
		anchor1MapState = CPU_MAPSTATE_UNIQUE;
	} else if (CPU_MAPSTATE_UNIQUE==readTagsMapState[1][1]) {
		// only R/2 tag 2 unique
		anchor1ReadIndex = readTags[1][1]; dedupRec->endpos = readEnd[1][1];
		anchor1MapState = CPU_MAPSTATE_UNIQUE;
	} else {
		// other cases
		if (CPU_MAPSTATE_REPEAT==readTagsMapState[0][0] || CPU_MAPSTATE_REPEAT==readTagsMapState[1][1]) {
			// either repeat, considered as repeat
			anchor1MapState = CPU_MAPSTATE_REPEAT;
			// we pick the read tag with longer mapped region
			anchor1ReadIndex = (readMappedLen[0][0]>=readMappedLen[1][1]) ? readTags[0][0] : readTags[1][1];
		} else {
			// not mapped, too many Ns, and not available
			anchor1MapState = CPU_MAPSTATE_NOT_MAPPED;
			// both unmapped, we keep the copy with the longer sequence, or there could really be no data
			if (CPU_MAPSTATE_NA==readTagsMapState[0][0]) {
				if (CPU_MAPSTATE_NA==readTagsMapState[1][1]) {
					// it will remain as -1
					//anchor1ReadIndex = -1;
				} else {
					anchor1ReadIndex = readTags[1][1]; dedupRec->endpos = readEnd[1][1];
				}
			} else {
				if (CPU_MAPSTATE_NA==readTagsMapState[1][1]) {
					anchor1ReadIndex = readTags[0][0]; dedupRec->endpos = readEnd[0][0];
				} else {
					if (w->bamRecs[readTags[0][0]].core.l_qseq>=w->bamRecs[readTags[1][1]].core.l_qseq) {
						anchor1ReadIndex = readTags[0][0]; dedupRec->endpos = readEnd[0][0];
					} else {
						anchor1ReadIndex = readTags[1][1]; dedupRec->endpos = readEnd[1][1];
					}
				}
			}
		}
	}
	
	// STATE 2.2 : proceed anchor#2
	uint8_t anchor2MapState = CPU_MAPSTATE_NA;
	int anchor2ReadIndex = -1;
	if (CPU_MAPSTATE_UNIQUE==readTagsMapState[1][0]) {
		if (CPU_MAPSTATE_UNIQUE==readTagsMapState[0][1]) {
			// both unique
			if (w->bamRecs[readTags[1][0]].core.tid == w->bamRecs[readTags[0][1]].core.tid) {
				// TODO: we pick the 5'-most read for now as it should has more reliable bases
				//       in future we might consider the length, etc
				if (readMappedLen[1][0]>=readMappedLen[0][1]) {
					anchor2ReadIndex = readTags[1][0]; dedupRec->mendpos = readEnd[1][0];
				} else {
					anchor2ReadIndex = readTags[0][1]; dedupRec->mendpos = readEnd[0][1];
				}
			} else {
				// different chromosome!
				// TODO: but we still pick one for processing
				anchorDiscordant |= 0x02;
				if (readMappedLen[1][0]>=readMappedLen[0][1]) {
					anchor2ReadIndex = readTags[1][0]; dedupRec->mendpos = readEnd[1][0];
				} else {
					anchor2ReadIndex = readTags[0][1]; dedupRec->mendpos = readEnd[0][1];
				}
			}
		} else {
			// only R/2 tag 1 unique
			anchor2ReadIndex = readTags[1][0]; dedupRec->mendpos = readEnd[1][0];
		}
		anchor2MapState = CPU_MAPSTATE_UNIQUE;
	} else if (CPU_MAPSTATE_UNIQUE==readTagsMapState[0][1]) {
		// only R/1 tag 2 unique
		anchor2ReadIndex = readTags[0][1]; dedupRec->mendpos = readEnd[0][1];
		anchor2MapState = CPU_MAPSTATE_UNIQUE;
	} else {
		// other cases
		if (CPU_MAPSTATE_REPEAT==readTagsMapState[1][0] || CPU_MAPSTATE_REPEAT==readTagsMapState[0][1]) {
			// either repeat, considered as repeat
			anchor2MapState = CPU_MAPSTATE_REPEAT;
			// we pick the read tag with longer mapped region
			if (readMappedLen[1][0]>=readMappedLen[0][1]) {
				anchor2ReadIndex = readTags[1][0]; dedupRec->mendpos = readEnd[1][0];
			} else {
				anchor2ReadIndex = readTags[0][1]; dedupRec->mendpos = readEnd[0][1];
			}
		} else {
			// not mapped, too many Ns, and not available
			anchor2MapState = CPU_MAPSTATE_NOT_MAPPED;
			// both unmapped, we keep the copy with the longer sequence, or there could really be no data
			if (CPU_MAPSTATE_NA==readTagsMapState[1][0]) {
				if (CPU_MAPSTATE_NA==readTagsMapState[0][1]) {
					// it will remain as -1
					//anchor2ReadIndex = -1;
				} else {
					anchor2ReadIndex = readTags[0][1]; dedupRec->mendpos = readEnd[0][1];
				}
			} else {
				if (CPU_MAPSTATE_NA==readTagsMapState[0][1]) {
					anchor2ReadIndex = readTags[1][0]; dedupRec->mendpos = readEnd[1][0];
				} else {
					if (w->bamRecs[readTags[1][0]].core.l_qseq>=w->bamRecs[readTags[0][1]].core.l_qseq) {
						anchor2ReadIndex = readTags[1][0]; dedupRec->mendpos = readEnd[1][0];
					} else {
						anchor2ReadIndex = readTags[0][1]; dedupRec->mendpos = readEnd[0][1];
					}
				}
			}
		}
	}
	
	// TODO: there might be a better logics to handle discordant to get more useful information
	
	// order the R/1 & R/2 with smaller genomic coordinate on the left, larger on the right
	// if only R/1 or R/2 is present, we move that to anchor1 position
	const bam1_core_t *c1 = NULL;
	const cpu_readid_t *ri1 = NULL;
	const bam1_core_t *c2 = NULL;
	const cpu_readid_t *ri2 = NULL;
	
	if (-1==anchor1ReadIndex) {
		if (-1==anchor2ReadIndex) {
			// no data, unlikely to happen
			fprintf(stderr, "[E::%s] Cannot locate bam index for both Read/1 and Read/2 bulk id#%d\n", __func__, i);
			return ;
		} else {
			// anchor1 NA, anchor2 avail
			dedupRec->swap = 1;
			anchor1ReadIndex = anchor2ReadIndex;
			anchor2ReadIndex = -1;
			dedupRec->endpos = dedupRec->mendpos;
			dedupRec->mendpos = 0;
			
			uint8_t anchorMapState = anchor1MapState; anchor1MapState = anchor2MapState; anchor2MapState = anchorMapState;
			
			c1 = &(w->bamRecs[anchor1ReadIndex].core);
			ri1 = &(w->readids[anchor1ReadIndex]);
			
			// adjust the deduplication record
			dedupRec->classification = ri1->classification & 0x07;
			dedupRec->finalLinkerType = ri1->finalLinkerType & 0x0F;
			dedupRec->qual = (c1->qual>0x7F) ? 0x7F : c1->qual;
			dedupRec->pairId = ri1->pairId;
			
			dedupRec->tid = c1->tid; dedupRec->pos = c1->pos;
			dedupRec->strand = (BAM_FREVERSE==(c1->flag&BAM_FREVERSE)) ? 1 : 0;
			dedupRec->linkerType = ri1->linker_type & 0x0F;
			dedupRec->uniq = (CPU_MAPSTATE_UNIQUE==anchor1MapState) ? 1 : 0;
			dedupRec->rtid = 2*ri1->readId+ri1->tagId+1;
			dedupRec->lineno = w->n_processed + anchor1ReadIndex; // TODO: 32 bits or 64 bits?
			
			// there is no need to set the value if we want default zeros
			dedupRec->mtid = -1;
			/*
			 dedupRec->mtid = 0; dedupRec->pos = 0;
			 dedupRec->mstrand = 0;
			 dedupRec->mlinkerType = 0;
			 dedupRec->muniq = 0;
			 dedupRec->mrtid = 0;
			 */
			
			w->tagGroups[i].bamRecs[0] = &(w->bamRecs[anchor1ReadIndex]);
		}
	} else {
		if (-1==anchor2ReadIndex) {
			// anchor1 avail, anchor2 NA
			// no swapping needed, proceed as-is
			c1 = &(w->bamRecs[anchor1ReadIndex].core);
			ri1 = &(w->readids[anchor1ReadIndex]);
			
			// adjust the deduplication record
			dedupRec->classification = ri1->classification & 0x07;
			dedupRec->finalLinkerType = ri1->finalLinkerType & 0x0F;
			dedupRec->qual = (c1->qual>0x7F) ? 0x7F : c1->qual;
			dedupRec->pairId = ri1->pairId;
			
			dedupRec->tid = c1->tid; dedupRec->pos = c1->pos;
			dedupRec->strand = (BAM_FREVERSE==(c1->flag&BAM_FREVERSE)) ? 1 : 0;
			dedupRec->linkerType = ri1->linker_type & 0x0F;
			dedupRec->uniq = (CPU_MAPSTATE_UNIQUE==anchor1MapState) ? 1 : 0;
			dedupRec->rtid = 2*ri1->readId+ri1->tagId+1;
			dedupRec->lineno = w->n_processed + anchor1ReadIndex; // TODO: 32 bits or 64 bits?
			
			// there is no need to set the value if we want default zeros
			dedupRec->mtid = -1;
			/*
			dedupRec->mtid = 0; dedupRec->pos = 0;
			dedupRec->mstrand = 0;
			dedupRec->mlinkerType = 0;
			dedupRec->muniq = 0;
			dedupRec->mrtid = 0;
			 */
			
			w->tagGroups[i].bamRecs[0] = &(w->bamRecs[anchor1ReadIndex]);
			
		} else {
			// anchor1 avail, anchor2 avail
			c1 = &(w->bamRecs[anchor1ReadIndex].core);
			ri1 = &(w->readids[anchor1ReadIndex]);
			c2 = &(w->bamRecs[anchor2ReadIndex].core);
			ri2 = &(w->readids[anchor2ReadIndex]);
			
			if (c2->tid<c1->tid) dedupRec->swap = 1;
			else if (c2->tid==c1->tid) {
				if (c2->pos<c1->pos) dedupRec->swap = 1;
				else if (c2->pos==c1->pos) {
					if (dedupRec->mendpos<dedupRec->endpos) dedupRec->swap = 1;
				}
			}

			if (dedupRec->swap) {
				const bam1_core_t *c = c1; c1 = c2; c2 = c;
				const cpu_readid_t *ri = ri1; ri1 = ri2; ri2 = ri;
				int32_t endpos = dedupRec->endpos; dedupRec->endpos = dedupRec->mendpos; dedupRec->mendpos = endpos;
				dedupRec->lineno = w->n_processed + anchor2ReadIndex; // TODO: 32 bits or 64 bits?
				dedupRec->mlineno = w->n_processed + anchor1ReadIndex; // TODO: 32 bits or 64 bits?
				
				uint8_t anchorMapState = anchor1MapState; anchor1MapState = anchor2MapState; anchor2MapState = anchorMapState;
			} else {
				dedupRec->lineno = w->n_processed + anchor1ReadIndex; // TODO: 32 bits or 64 bits?
				dedupRec->mlineno = w->n_processed + anchor2ReadIndex; // TODO: 32 bits or 64 bits?
			}
			
			// adjust the deduplication record
			dedupRec->classification = ri1->classification & 0x07;
			dedupRec->finalLinkerType = ri1->finalLinkerType & 0x0F;
			// we keep the lower of the two scores so that we do not over-reported in the case of repeats
			uint8_t qual = (c1->qual <= c2->qual) ? c1->qual : c2->qual;
			dedupRec->qual = (qual>0x7F) ? 0x7F : qual;
			dedupRec->pairId = ri1->pairId;
			
			dedupRec->tid = c1->tid; dedupRec->pos = c1->pos;
			dedupRec->strand = (BAM_FREVERSE==(c1->flag&BAM_FREVERSE)) ? 1 : 0;
			dedupRec->linkerType = ri1->linker_type & 0x0F;
			dedupRec->uniq = (CPU_MAPSTATE_UNIQUE==anchor1MapState) ? 1 : 0;
			dedupRec->rtid = 2*ri1->readId+ri1->tagId+1;
			
			dedupRec->mtid = c2->tid; dedupRec->mpos = c2->pos;
			dedupRec->mstrand = (BAM_FREVERSE==(c2->flag&BAM_FREVERSE)) ? 1 : 0;
			dedupRec->mlinkerType = ri2->linker_type & 0x0F;
			dedupRec->muniq = (CPU_MAPSTATE_UNIQUE==anchor2MapState) ? 1 : 0;
			dedupRec->mrtid = 2*ri2->readId+ri2->tagId+1;

			w->tagGroups[i].bamRecs[0] = &(w->bamRecs[anchor1ReadIndex]);
			w->tagGroups[i].bamRecs[1] = &(w->bamRecs[anchor2ReadIndex]);
			
		}
	}

	
	// STATE 3: record interaction-pair
	if (CPU_MAPSTATE_UNIQUE==anchor1MapState && CPU_MAPSTATE_UNIQUE==anchor2MapState) {
		// [UU] both unique
		if (0==anchorDiscordant) {
			w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_UU;
			
			// TODO: we don't really need new read id, unless we want to merge pairing with deduplication
			// TODO: introduce a YT:Z tag so that we know the decision
			//G_PAIRING_TAG_UU
			if (w->nUpdateBam) {
				// let's copy most of the structure
				w->tagGroups[i].releaseBamRecs = 1;
				w->tagGroups[i].bamRecs[0] = bam_dup1(&(w->bamRecs[anchor1ReadIndex]));
				w->tagGroups[i].bamRecs[1] = bam_dup1(&(w->bamRecs[anchor2ReadIndex]));
				bam1_t *anchor1Bam = w->tagGroups[i].bamRecs[0];
				bam1_t *anchor2Bam = w->tagGroups[i].bamRecs[1];
				
				// change the read name! for /1 and /2 !!!
				// combined the ids for debugging purposes
				char szQName[G_MAX_QNAME_LEN];
				int nRet = buildPairReadId (&(w->readids[anchor1ReadIndex]), &(w->readids[anchor2ReadIndex]), szQName);
				if (nRet<=0) {
					// TODO: this is considered fatal, cannot continue
					//if (bwa_verbose >= 1)
					fprintf(stderr, "[E::%s] Unmatching paired read ids \"%s\" and \"%s\".\n", __func__, bam1_qname(anchor1Bam), bam1_qname(anchor2Bam));
				} else {
					// adjust the auxiliary data, copy the new qname, and all the rest of the data
#if 1
					// TODO: checking if this is the source of memory leaks
					// handle anchor1 QNAME
					int lengthExcludeQName = anchor1Bam->data_len - anchor1Bam->core.l_qname;
					int newLength = lengthExcludeQName + nRet + 1;
					int bufferLength = newLength; kroundup32(bufferLength);
					uint8_t *data = (uint8_t *)malloc(bufferLength);
					
					if (NULL==data) {
						// TODO: this is considered fatal, cannot continue
						//if (bwa_verbose >= 1)
						fprintf(stderr, "[E::%s] Cannot allocate memory for new data, anchor 1 read ids \"%s\".\n", __func__, bam1_qname(anchor1Bam));
					} else {
						// copy new QName
						memcpy(data, szQName, nRet + 1);
						// copy all remaining data
						memcpy(data + nRet + 1, bam1_cigar(anchor1Bam), lengthExcludeQName);
						anchor1Bam->core.l_qname = nRet + 1;
						anchor1Bam->data_len = newLength;
						anchor1Bam->m_data = bufferLength;
						free(anchor1Bam->data);
						anchor1Bam->data = data;
					}
					
					// handle anchor2 QNAME
					lengthExcludeQName = anchor2Bam->data_len - anchor2Bam->core.l_qname;
					newLength = lengthExcludeQName + nRet + 1;
					bufferLength = newLength; kroundup32(bufferLength);
					data = (uint8_t *)malloc(bufferLength);
					
					if (NULL==data) {
						// TODO: this is considered fatal, cannot continue
						//if (bwa_verbose >= 1)
						fprintf(stderr, "[E::%s] Cannot allocate memory for new data, anchor 2 read ids \"%s\".\n", __func__, bam1_qname(anchor2Bam));
					} else {
						// copy new QName
						memcpy(data, szQName, nRet + 1);
						// copy all remaining data
						memcpy(data + nRet + 1, bam1_cigar(anchor2Bam), lengthExcludeQName);
						anchor2Bam->core.l_qname = nRet + 1;
						anchor2Bam->data_len = newLength;
						anchor2Bam->m_data = bufferLength;
						free(anchor2Bam->data);
						anchor2Bam->data = data;
					}
#endif
				}
				
				// FLAG
				anchor1Bam->core.flag |= BAM_FPAIRED; //FLAG, read paired
				if (BAM_FREVERSE & anchor2Bam->core.flag) anchor1Bam->core.flag |= BAM_FMREVERSE;
				anchor1Bam->core.flag |= BAM_FREAD1;
				
				anchor2Bam->core.flag |= BAM_FPAIRED; //FLAG, read paired
				if (BAM_FREVERSE & anchor1Bam->core.flag) anchor2Bam->core.flag |= BAM_FMREVERSE;
				anchor2Bam->core.flag |= BAM_FREAD2;
				
				// TODO: need to adjust the MAPQ as per a pair-end, rather than single-end
				// anchor1Bam->core.qual; #MAPQ:
				// anchor2Bam->core.qual; #MAPQ:
				
				anchor1Bam->core.mtid = anchor2Bam->core.tid; //RNEXT
				anchor1Bam->core.mpos = anchor2Bam->core.pos; //PNEXT
				anchor2Bam->core.mtid = anchor1Bam->core.tid; //RNEXT
				anchor2Bam->core.mpos = anchor1Bam->core.pos; //PNEXT
				
				if (anchor1Bam->core.tid == anchor2Bam->core.tid) {
					// intra-chromosomal, calculate the template span
					uint32_t pos1_1 = anchor1Bam->core.pos;
					uint32_t pos1_2 = bam_calend(&(anchor1Bam->core), bam1_cigar(anchor1Bam));
					uint32_t pos2_1 = anchor2Bam->core.pos;
					uint32_t pos2_2 = bam_calend(&(anchor2Bam->core), bam1_cigar(anchor2Bam));
					uint32_t pos_1 = (pos1_1 < pos2_1) ? pos1_1 : pos2_1;
					uint32_t pos_2 = (pos1_2 > pos2_2) ? pos1_2 : pos2_2;
					int32_t tlen = ((pos_1>pos_2) ? pos_1 - pos_2 : pos_2 - pos_1) + 1;
					if (anchor1Bam->core.pos<anchor2Bam->core.pos) {
						anchor1Bam->core.isize = tlen;
						anchor2Bam->core.isize = -1*tlen;
					} else {
						anchor1Bam->core.isize = -1*tlen;
						anchor2Bam->core.isize = tlen;
					}
				} else {
					// inter-chromosomal, use 0 as template length
					anchor1Bam->core.isize = anchor2Bam->core.isize = 0; //TLEN
				}
				
				// introduce a YT:Z tag so that we know the decision
				uint8_t *yt = bam_aux_get(anchor1Bam, G_PAIRING_TAG);
				if (0!=yt) bam_aux_del(anchor1Bam, yt); // have to remove
				bam_aux_append(anchor1Bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_UU);
				yt = bam_aux_get(anchor2Bam, G_PAIRING_TAG);
				if (0!=yt) bam_aux_del(anchor2Bam, yt); // have to remove
				bam_aux_append(anchor2Bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_UU);
			}
			
		} else {
			// discordant for one anchor or both the anchors
			// we cannot use this, but we keep them in separate file for investigation
			// introduce a YT:Z tag so that we know the decision
			w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_UU_DISCORDANT;
			
			//G_PAIRING_TAG_UU_Discordant
			// TODO: we exclude discordant from deduplication and uniqueness consideration
			//       i.e. it is considered duplicated and non-unique
			// NOT necessary if we have chosen the best mapping!!!
			if (w->nUpdateBam) {
				for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
					bam1_t *bam = &(w->bamRecs[k]);
					uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
					if (0!=yt) bam_aux_del(bam, yt); // have to remove
					bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_UU_Discordant);
				}
			}
		}
	} else if (CPU_MAPSTATE_UNIQUE==anchor1MapState || CPU_MAPSTATE_UNIQUE==anchor2MapState) {
		// [UxxU] only 1 unique, output all the tags
		// introduce a YT:Z tag so that we know the decision
		w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_UxxU;
		
		// TODO: how do we tell it is Ux or xU later? we need to know which is unique!
		//       consider removing "duplicated"
		if (w->nUpdateBam) {
			if (CPU_MAPSTATE_UNIQUE==anchor1MapState) {
				//G_PAIRING_TAG_Ux
				if (CPU_MAPSTATE_REPEAT==anchor2MapState) {
					// U-Repeat
				} else {
					// U-Unmapped/Ns
				}
				for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
					bam1_t *bam = &(w->bamRecs[k]);
					uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
					if (0!=yt) bam_aux_del(bam, yt); // have to remove
					bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_Ux);
				}
			} else {
				//G_PAIRING_TAG_xU
				if (CPU_MAPSTATE_REPEAT==anchor1MapState) {
					// Repeat-U
				} else {
					// Unmapped/Ns-U
				}
				for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
					bam1_t *bam = &(w->bamRecs[k]);
					uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
					if (0!=yt) bam_aux_del(bam, yt); // have to remove
					bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_xU);
				}
			}
		}
		
	} else {
		// [xx], [nn] separate into two classes, output all the tags
		if (CPU_MAPSTATE_REPEAT==anchor1MapState || CPU_MAPSTATE_REPEAT==anchor2MapState) {
			// [xx]
			// introduce a YT:Z tag so that we know the decision
			w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_xx;
			
			// G_PAIRING_TAG_xx
			if (w->nUpdateBam) {
				for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
					bam1_t *bam = &(w->bamRecs[k]);
					uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
					if (0!=yt) bam_aux_del(bam, yt); // have to remove
					bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_xx);
				}
			}
		} else {
			// [nn]
			// introduce a YT:Z tag so that we know the decision
			w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_nn;
			
			// G_PAIRING_TAG_nn
			// we exclude both unmapped from deduplication and uniqueness consideration
			// i.e. it is considered duplicated and non-unique
			if (w->nUpdateBam) {
				for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
					bam1_t *bam = &(w->bamRecs[k]);
					uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
					if (0!=yt) bam_aux_del(bam, yt); // have to remove
					bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_nn);
				}
			}
		}
	}
}

void dedup_process_alns(const dedup_opt_t *opt, uint32_t n_processed, int n, bam1_t *bamRecs, cpu_readid_t *readids, uint32_t *pNumTagGroups, tag_group_t *tagGroups, cpu_dedup_record_t *dedupRecords, CPUBamBuffer_t* bamsBuffer, uint32_t nUpdateBam)
{
	int i;
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	if (n<=0) return;
	
	w.opt = opt;
	w.n_processed = n_processed; //TODO:
	//TODO: w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;

	w.bamRecs = bamRecs;
	w.readids = readids;
	w.nUpdateBam = nUpdateBam;
	
	// parse the read id
	kt_for(opt->n_threads, readid_worker, &w, n);
	
	// tag grouping
	// TODO: parallelize later
	int numTagGroups = 0;
	long pairId = readids[0].pairId;
	tagGroups[numTagGroups].readIndex = 0; tagGroups[numTagGroups].n_reads = 1;
	int currTagGroup = numTagGroups; numTagGroups++;
	for(i=1; i<n; ++i) {
		if (pairId == readids[i].pairId) {
			// check for overflow
			if (127==tagGroups[currTagGroup].n_reads) {
				fprintf(stderr, "[E::%s] There are >127 tag records for paid id %lu", __func__, pairId);
			}
			tagGroups[currTagGroup].n_reads++;
		} else {
			pairId = readids[i].pairId;
			tagGroups[numTagGroups].readIndex = i; tagGroups[numTagGroups].n_reads = 1;
			currTagGroup++; numTagGroups++;
		}
	}
	
	if (!bamsBuffer->lastBlock) {
		// this is not the last block, and thus we cannot guarantee the completeness of the last tagGroup!
		// let's cache the last block for the next iteration
		numTagGroups--;
		add_CPUBamBuffer(&(bamRecs[tagGroups[numTagGroups].readIndex]), tagGroups[numTagGroups].n_reads, bamsBuffer);
	}
	*pNumTagGroups = numTagGroups;
	
	// pairing selection
	w.nTagGroups = numTagGroups;
	w.tagGroups = tagGroups;
	w.dedupRecords = dedupRecords;
	kt_for(opt->n_threads, dedup_worker, &w, numTagGroups);
	
	// TODO: accumulate all the span information
	for(i=0; i<opt->n_threads; ++i) {
	}
}

int updateDedupTotalRecords(const char *dedupKeyFilename, uint32_t n_pairProcessed) {
	const char *mode = "r+b";
	FILE *fd = fopen(dedupKeyFilename, mode);
	if (fd == NULL) {
		fprintf(stderr, "[E::%s] fail to open '%s' for update: %s\n", __func__, dedupKeyFilename, strerror (errno));
		return errno;
	}
	
	// TODO: header processing synchronization
	fseek(fd, 3*sizeof(uint32_t), SEEK_SET);
	err_fwrite(&n_pairProcessed, sizeof(n_pairProcessed), 1, fd);// filler for total number of records
	fclose(fd);
	
	return 0;
}

// TODO:
cpu_dedup_record_t *readDedupKeys(const char *dedupKeyFilename, uint32_t *nDedupKeys) {
	*nDedupKeys = 0;
	cpu_dedup_record_t *dedupKeys = NULL;
	
	// open file and decide on amount of storage to allocate
	const char *mode = "rb";
	FILE *fd = fopen(dedupKeyFilename, mode);
	if (fd == NULL) {
		fprintf(stderr, "[E::%s] fail to open '%s' for update: %s\n", __func__, dedupKeyFilename, strerror (errno));
		return NULL;
	}
	
	uint32_t nValue = 0;
	err_fread_noeof(&nValue, sizeof(nValue), 1, fd); // cpu version
	err_fread_noeof(&nValue, sizeof(nValue), 1, fd); // sort state: unsorted=0, chunk sorted=1, fully sorted=2
	err_fread_noeof(&nValue, sizeof(nValue), 1, fd); // record size
	err_fread_noeof(&nValue, sizeof(nValue), 1, fd); // filler for total number of records
	*nDedupKeys = nValue;
	if (*nDedupKeys) {
		dedupKeys = (cpu_dedup_record_t *) calloc(*nDedupKeys, sizeof(cpu_dedup_record_t));
		if (dedupKeys) {
			// chunk read the data
			
			// TODO: more error checking needed
			//       what if it is truncated? what if there are more records than allocated?
			uint32_t nRead = 0;
			while (nRead<*nDedupKeys) {
				uint32_t nChunkRecs = 0;
				err_fread_noeof(&nChunkRecs, sizeof(nChunkRecs), 1, fd); // number of records in chunk
				uint32_t nSorted = 0;
				err_fread_noeof(&nSorted, sizeof(nSorted), 1, fd); // chunk sorted?
				if ((nRead+nChunkRecs)>*nDedupKeys) {
					fprintf(stderr, "[E::%s] There are more records than recorded. Specified: %u, Reading: %u.\n", __func__, *nDedupKeys, nRead+nChunkRecs);
					*nDedupKeys = nRead;
					break;
				}
				err_fread_noeof(&(dedupKeys[nRead]), sizeof(cpu_dedup_record_t), nChunkRecs, fd);
				nRead += nChunkRecs;
			}
		} else {
			fprintf(stderr, "[E::%s] out of memory to read %u records from '%s'e: %s\n", __func__, *nDedupKeys, dedupKeyFilename, strerror (errno));
		}
	}
	
	// close file and return prep'd data
	fclose(fd);
	return dedupKeys;
}

int cmpPET(const void *a, const void *b) {
	const cpu_dedup_record_t *ia = (const cpu_dedup_record_t *)a;
	const cpu_dedup_record_t *ib = (const cpu_dedup_record_t *)b;
	
	// sort by chrom1, start1, end1,
	int32_t i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	// follow by chrom2, start2, end2
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;

	// strand1, strand2, final linker, linker1, linker2, score
	// TODO: there are 8x2 comparison here which can be simplified to 4 or 2
#if 0
	i32comp = ((CPU_DEDUP_METADATA_SORT_KEY(ia->metaData)))-(CPU_DEDUP_METADATA_SORT_KEY(ib->metaData))); if (0!=i32comp) return i32comp;
	i32comp = ((CPU_DEDUP_METADATA_QUAL(ib->metaData)))-(CPU_DEDUP_METADATA_QUAL(ia->metaData)));
	return i32comp;
#else
	i32comp = ia->strand - ib->strand; if (0!=i32comp) return i32comp;
	i32comp = ia->mstrand - ib->mstrand; if (0!=i32comp) return i32comp;
	i32comp = ia->classification - ib->classification; if (0!=i32comp) return i32comp;
	i32comp = ia->finalLinkerType - ib->finalLinkerType; if (0!=i32comp) return i32comp;
	i32comp = ia->linkerType - ib->linkerType; if (0!=i32comp) return i32comp;
	i32comp = ia->mlinkerType - ib->mlinkerType; if (0!=i32comp) return i32comp;
	// TODO: qual is higher first!
	//i32comp = ia->qual - ib->qual; if (0!=i32comp) return i32comp;
	i32comp = ib->qual - ia->qual; if (0!=i32comp) return i32comp;
	i32comp = ia->lineno - ib->lineno; //if (!i32comp) return i32comp;
	return i32comp;
#endif
}

int isDuplicate(cpu_dedup_record_t *ia, cpu_dedup_record_t *ib) {
	// sort by chrom1, start1, end1,
	int32_t i32comp = ia->tid - ib->tid; if (0!=i32comp) return 0;
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return 0;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return 0;
	
	// follow by chrom2, start2, end2
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return 0;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return 0;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return 0;
	
	// strand1, strand2, final linker, linker1, linker2
	// TODO: optimized by reducing the number of comparisons
	i32comp = ia->strand - ib->strand; if (0!=i32comp) return 0;
	i32comp = ia->mstrand - ib->mstrand; if (0!=i32comp) return 0;
	i32comp = ia->classification - ib->classification; if (0!=i32comp) return 0;
	i32comp = ia->finalLinkerType - ib->finalLinkerType; if (0!=i32comp) return 0;
	i32comp = ia->linkerType - ib->linkerType; if (0!=i32comp) return 0;
	i32comp = ia->mlinkerType - ib->mlinkerType; if (0!=i32comp) return 0;
	return 1;
}

static inline void mark_duplicate_line(uint32_t *b, uint32_t p) {
	uint32_t i32b = p >> 5;
	uint8_t o32b = p & 0x1F ;
	b[i32b] |= (1 << o32b);
}

static inline int is_line_duplicated(uint32_t *b, uint32_t p) {
	uint32_t i32b = p >> 5;
	uint8_t o32b = p & 0x1F ;
	return (b[i32b] & (1 << o32b));
}

// TODO: we should get a running total of number of duplicates
//       alternatively, parallelized bit counting
uint32_t *markDuplicates(cpu_dedup_record_t *dedupKeys, int nDedupKeys, uint32_t n_processed, uint32_t *pUniquePETs, uint32_t *pUnMappedPETs) {
	int i;
	
	uint32_t nUniquePETs = 0;
	uint32_t nUnMappedPETs = 0;
	
	// allocate the bits string
	if (0==n_processed) return NULL;
	
	uint32_t n32bSize = n_processed >> 5;
	if (n_processed & 0x1F) n32bSize++;
	uint32_t *duplicates = (uint32_t *) calloc(n32bSize, sizeof(uint32_t));
	if (!duplicates) return NULL;
	
	cpu_dedup_record_t *pRec = &dedupKeys[0];
	if (-1==pRec->tid && -1==pRec->mtid) {
		pRec->duplicated = 1;
		nUnMappedPETs++;
	} else {
		nUniquePETs++;
	}
	// pRec->duplicated = 0; // first record can never be duplicated!
	for (i=1; i<nDedupKeys; ++i) {
		// if both anchors are not mapped
		if (-1==dedupKeys[i].tid && -1==dedupKeys[i].mtid) {
			dedupKeys[i].duplicated = 1;
			nUnMappedPETs++;
		} else {
			if (isDuplicate(pRec, &dedupKeys[i])) {
				dedupKeys[i].duplicated = 1;
				
				// TODO: record the duplicated entries
				if (-1!=dedupKeys[i].tid) mark_duplicate_line(duplicates,dedupKeys[i].lineno);
				if (-1!=dedupKeys[i].mtid) mark_duplicate_line(duplicates,dedupKeys[i].mlineno);
			} else {
				nUniquePETs++;
			}
		}
		pRec = &dedupKeys[i];
	}
	
	*pUniquePETs = nUniquePETs;
	*pUnMappedPETs = nUnMappedPETs;
	
	return duplicates;
}


int main_dedup(int argc, char *argv[])
{
	int c, n, i, j, k;
	dedup_opt_t *opt;
	
	bam1_t *bamRecs;
	CPUBamBuffer_t bamsBuffer;
	
	double t_diff;
	double t_timepointIO, t_timepointProcessing;
	double t_diffIO, t_diffProcessing;
	
	uint32_t n_processed = 0;
	uint32_t n_pairProcessed = 0;
	
	int is_bamin = 1;
	char in_mode[5];
	char *fn_list = 0;
	samfile_t *in = 0;
	
	int ret = 0;
	
	int compress_level = -1;
	int is_header = 0;
	int is_bamout = 1;
	char out_mode[5];
	samfile_t *outfds[CPU_DEDUP_OUTPUT_COUNT] = {0};
	
	FILE *outfd = NULL;
	kstring_t dedupKeyFilename = {0,0,0};

	init_CPUBamBuffer(&bamsBuffer);
	
	opt = dedup_opt_init();
	strcpy(in_mode, "r");
	strcpy(out_mode, "w");
	while ((c = getopt(argc, argv, "cdgSs:t:O:l:u:b:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'l') opt->outputSortedList = atoi(optarg), opt->outputSortedList = opt->outputSortedList > 0 ? opt->outputSortedList : 0;
		else if (c == 'c') opt->calculateLibraryComplexity = 1;
		else if (c == 'd') opt->disponoParallel = 1;
		else if (c == 'g') opt->toGroup = 1;
		else if (c == 'S') is_bamin = 0;
		else if (c == 'O') { ks_set(optarg, &(opt->outputPrefix)); }
		else if (c == 's') opt->selfLigation = atoi(optarg), opt->selfLigation = opt->selfLigation > 0 ? opt->selfLigation : G_SELF_LIGATION;
		else {
			dedup_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	
	if (opt->n_threads < 1) opt->n_threads = 1;
	
	if (is_bamin) strcat(in_mode, "b");
	if (is_bamout) strcat(out_mode, "b");
	if (is_header) strcat(out_mode, "h");
	if (compress_level >= 0) {
		char tmp[2];
		tmp[0] = compress_level + '0'; tmp[1] = '\0';
		strcat(out_mode, tmp);
	}
	
	
	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu dedup [options] <in.sam/.bam>\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -d         use parallel sort\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -s INT     self-ligation distance [%d bp]\n", opt->selfLigation);
		fprintf(stderr, "       -S         input sam format\n");
		fprintf(stderr, "       -g         output different categories read to separate file\n");
		fprintf(stderr, "       -l         output sorted PETs to stdout\n");
		fprintf(stderr, "       -c         output library complexity metric\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		dedup_opt_terminate(opt);
		free(opt);
		return 1;
	}
	
	// open file handlers
	if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open \"%s\" for reading.\n", __func__, argv[optind]);
		
		free(fn_list);
		dedup_opt_terminate(opt);
		free(opt);
		
		return 1;
	}
	if (in->header == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to read the header from \"%s\".\n", __func__, argv[optind]);
		
		free(fn_list);
		// close files, free and return
		samclose(in);
		
		dedup_opt_terminate(opt);
		free(opt);
		
		return 1;
	}
	
	// construct the out filename from the in filename
	{
		if (NULL==opt->outputPrefix.s || 0==opt->outputPrefix.l) {
			ks_set(argv[optind], &(opt->outputPrefix));
			char *pch = strrchr(opt->outputPrefix.s, '.');
			if (NULL==pch) {
				// there is no extension, use as-is
			} else {
				if (0==strcmp(pch, ".bam") || 0==strcmp(pch, ".sam")) {
					*pch = '\0';
					opt->outputPrefix.l -= 4;
				} else if (0==strcmp(pch, ".gz")) {
					*pch = '\0';
					opt->outputPrefix.l -= 3;
					char *pch1 = strrchr(opt->outputPrefix.s, '.');
					if (NULL!=pch1) {
						// remove one more extension
						*pch1 = '\0';
						opt->outputPrefix.l = pch1 - opt->outputPrefix.s;
					}
				} else {
					*pch = '\0';
					opt->outputPrefix.l = pch - opt->outputPrefix.s;
				}
			}
		}
	}
	
	if (opt->toGroup) {
		if (init_CPDedup_Outputs (opt->outputPrefix.s, out_mode, in->header, outfds, opt->selfLigation)) {
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open for writing.\n", __func__);
			
			free(fn_list);
			// close files, free and return
			samclose(in);
			
			dedup_opt_terminate(opt);
			free(opt);
			
			return 1;
		}
		
		{
			if (opt->n_threads > 1) samthreads(outfds[CPU_DEDUP_OUTPUT_REPRESENTATIVE], opt->n_threads, 256);
			if (opt->n_threads > 1) samthreads(outfds[CPU_DEDUP_OUTPUT_DUPLICATED], opt->n_threads, 256);
		}
	}

	// we keep the data that we are going to sort in a binary files for single read processing
	// in future, if we have decided on merge sort, we can perform other form of combining the results
	// also, we can already sort the entries before writing them into the binary files
	{
		ksprintf(&dedupKeyFilename, "%s.cpu.dedup", opt->outputPrefix.s);
		
		const char *mode = "wb";
		outfd = fopen(dedupKeyFilename.s, mode);
		if (outfd == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing: %s\n", __func__, dedupKeyFilename.s, strerror (errno));
			if (opt->toGroup) {
				terminate_CPDedup_Outputs(outfds);
			}
			dedup_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	
	// we write header information for future usage
	// TODO: can we safely read text lines from a binary files?
	//       if so, we will like to record the software version and command line
	// TODO: get a unsigned int32 version of the record
	uint32_t nValue = (0x00<<16 | 0x01 <<8 | 0x00); err_fwrite(&nValue, sizeof(nValue), 1, outfd); // cpu version
	nValue = 0; err_fwrite(&nValue, sizeof(nValue), 1, outfd); // sort state: unsorted=0, chunk sorted=1, fully sorted=2
	nValue = sizeof(cpu_dedup_record_t); err_fwrite(&nValue, sizeof(nValue), 1, outfd); // record size
	err_fwrite(&n_pairProcessed, sizeof(n_pairProcessed), 1, outfd);// filler for total number of records
	
	// PROCESSING
	if (argc == optind + 1) { // convert/print the entire file
		if (CPU_DEDUP_DEBUG_CONVERT_LIST==(CPU_DEDUP_DEBUG_CONVERT_LIST&opt->outputSortedList)) {
			fprintf(stdout, "\n\n===\tCONVERTED LIST\n");
		}
		
		//int nNumRequested = opt->chunk_size * opt->n_threads;
		int nNumRequested = opt->chunk_size;
		t_timepointProcessing = realtime();
		while ((bamRecs = readCPUBam(&n, nNumRequested, in, &bamsBuffer)) != 0) {
			t_timepointIO = realtime();
			t_diffIO = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] read %d records..\n", __func__, n-bamsBuffer.unprocessed);
			
			// processing the pairing
#if 0
			fprintf(stderr, "readids allocated with %lu bytes, %.2fMB\n", n*sizeof(cpu_readid_t), n*sizeof(cpu_readid_t)*1.0f/1024.0/1024.0);
			fprintf(stderr, "tagGroups allocated with %lu bytes, %.2fMB\n", n*sizeof(tag_group_t), n*sizeof(tag_group_t)*1.0f/1024.0/1024.0);
			fprintf(stderr, "dedupRecords allocated with %lu bytes, %.2fMB\n", n*sizeof(cpu_dedup_record_t), n*sizeof(cpu_dedup_record_t)*1.0f/1024.0/1024.0);
#endif
			cpu_readid_t *readids = calloc(n, sizeof(cpu_readid_t));
			tag_group_t *tagGroups = calloc(n, sizeof(tag_group_t));
			cpu_dedup_record_t *dedupRecords = calloc(n, sizeof(cpu_dedup_record_t));
			
			uint32_t numTagGroups = 0;
			dedup_process_alns(opt, n_processed, n, bamRecs, readids, &numTagGroups, tagGroups, dedupRecords, &bamsBuffer, 0);
			t_timepointProcessing = realtime();
			t_diffProcessing = t_timepointProcessing - t_timepointIO;
			
			if (opt->toGroup) {
				// TODO: write the bam results
				// TODO: for testing purposes
				// TODO: need to make sure that we have either 1 record (SE) or 2 records (PE)

				if (CPU_DEDUP_DEBUG_CONVERT_LIST==(CPU_DEDUP_DEBUG_CONVERT_LIST&opt->outputSortedList)) {
					for(i=0; i<numTagGroups; ++i) {
						cpu_dedup_record_t *dedupRec = &(dedupRecords[i]);
						
						kstring_t str;
						str.l = str.m = 0; str.s = 0;
						
						// pair id
						kputul(dedupRec->pairId, &str); kputc('\t', &str);
						
						// mapping quality
						kputw(dedupRec->qual, &str); kputc('\t', &str);
						
						// linker 1
						const char *text = shortLinkerTypeToStr(dedupRec->linkerType); kputs(text, &str); kputc('\t', &str);
						// linker 2
						text = shortLinkerTypeToStr(dedupRec->mlinkerType); kputs(text, &str); kputc('\t', &str);
						
						// final linker
						text = shortLinkerTypeToStr(dedupRec->finalLinkerType); kputs(text, &str); kputc('\t', &str);
						
						// classification
						text = shortOutputClassToStr(dedupRec->classification); kputs(text, &str); kputc('\t', &str);
						
						// read 1/tag 1
						kputw(dedupRec->rtid, &str); kputc('\t', &str);
						
						// read 2/tag 2
						kputw(dedupRec->mrtid, &str); kputc('\t', &str);
						
						if (dedupRec->tid < 0) kputsn("*\t0\t0\t+\t", 8, &str);
						else {
							if (in->header) kputs(in->header->target_name[dedupRec->tid] , &str);
							else kputw(dedupRec->tid, &str);
							kputc('\t', &str);
							kputw(dedupRec->pos + 1, &str); kputc('\t', &str);
							kputw(dedupRec->endpos, &str); kputc('\t', &str);
							kputsn(dedupRec->strand ? "-\t" : "+\t", 2, &str);
						}
						
						if (dedupRec->mtid < 0) kputsn("*\t0\t0\t+\t", 8, &str);
						else {
							if (in->header) kputs(in->header->target_name[dedupRec->mtid] , &str);
							else kputw(dedupRec->mtid, &str);
							kputc('\t', &str);
							kputw(dedupRec->mpos + 1, &str); kputc('\t', &str);
							kputw(dedupRec->mendpos, &str); kputc('\t', &str);
							kputsn(dedupRec->mstrand ? "-\t" : "+\t", 2, &str);
						}
						
						// we recompute the insert size just in case it is unavailable
						int32_t iSize = 0;
						if (dedupRec->tid == dedupRec->mtid && -1!=dedupRec->tid) iSize = dedupRec->mendpos - dedupRec->pos;
						kputw(iSize, &str); kputc('\t', &str);
						
						kputsn(bam1_qname(tagGroups[i].bamRecs[0]), tagGroups[i].bamRecs[0]->core.l_qname-1, &str); kputc('\t', &str);
						
						kputw(dedupRec->swap, &str); //kputc('\t', &str);
						
						//TODO: for debugging
						kputc('\t', &str); kputul(dedupRec->lineno, &str);
						kputc('\t', &str); kputul(dedupRec->mlineno, &str);
						kputc('\n', &str);
						
						fprintf(stdout, "%s", str.s);
						free(str.s);
					}
				}
				
			}
			
			// TODO: write the binary data file
			err_fwrite(&numTagGroups, sizeof(numTagGroups), 1, outfd); // write the number of record filler for total number of records
			uint32_t nSorted = 0; err_fwrite(&nSorted, sizeof(nSorted), 1, outfd); // sort state: unsorted=0, chunk sorted=1, fully sorted=2
			err_fwrite(dedupRecords, sizeof(cpu_dedup_record_t), numTagGroups, outfd);
			
			t_timepointIO = realtime();
			t_diffIO += (t_timepointIO - t_timepointProcessing);
			
			
			// clean up
			for(i=0; i<n; ++i) free(bamRecs[i].data);
			free(bamRecs);
			free(readids);
			free(dedupRecords);
			free(tagGroups);
			
			n_processed += (n-bamsBuffer.unprocessed);
			n_pairProcessed += numTagGroups;
			if (bwa_verbose >= 3) {
				t_diff = realtime() - G_t_real;
				fprintf(stderr, "[M::%s] %u tags %u pairs processed, %.0f tags/sec %.0f pairs/sec, %.2f min, i/o %.2f sec, processing %.2f sec..\n", __func__, n_processed, n_pairProcessed, 1.0*n_processed/t_diff, 1.0*n_pairProcessed/t_diff, t_diff/60.0, t_diffIO, t_diffProcessing);
			}
			
			t_timepointProcessing = realtime();
		}
		fclose(outfd);
		
		// need to update the number of records
		t_timepointProcessing = realtime();
		updateDedupTotalRecords(dedupKeyFilename.s, n_pairProcessed);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] update record count, i/o %.2f sec (%.2f min)..\n", __func__, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// we next read all the binary records for sorting
		uint32_t nDedupKeys = 0;
		t_timepointProcessing = realtime();
		cpu_dedup_record_t *dedupKeys = readDedupKeys(dedupKeyFilename.s, &nDedupKeys);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] read %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nDedupKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		if (!dedupKeys) {
			fprintf(stderr, "[E::%s] No duplicated key reads, terminating..\n", __func__);
			// TODO: clean up
		}
		
		if (CPU_DEDUP_DEBUG_READ_LIST==(CPU_DEDUP_DEBUG_READ_LIST&opt->outputSortedList)) {
			// for debugging purposes only
			t_timepointProcessing = realtime();
			fprintf(stdout, "\n\n===\tREAD LIST\n");
			for(i=0; i<nDedupKeys; ++i) {
				
				cpu_dedup_record_t *dedupRec = &(dedupKeys[i]);
				
				kstring_t str;
				str.l = str.m = 0; str.s = 0;
				
				// pair id
				kputul(dedupRec->pairId, &str); kputc('\t', &str);
				
				// mapping quality
				kputw(dedupRec->qual, &str); kputc('\t', &str);
				
				// linker 1
				const char *text = shortLinkerTypeToStr(dedupRec->linkerType); kputs(text, &str); kputc('\t', &str);
				// linker 2
				text = shortLinkerTypeToStr(dedupRec->mlinkerType); kputs(text, &str); kputc('\t', &str);
				
				// final linker
				text = shortLinkerTypeToStr(dedupRec->finalLinkerType); kputs(text, &str); kputc('\t', &str);
				
				// classification
				text = shortOutputClassToStr(dedupRec->classification); kputs(text, &str); kputc('\t', &str);
				
				// read 1/tag 1
				kputw(dedupRec->rtid, &str); kputc('\t', &str);
				
				// read 2/tag 2
				kputw(dedupRec->mrtid, &str); kputc('\t', &str);
				
				if (dedupRec->tid < 0) kputsn("*\t0\t0\t+\t", 8, &str);
				else {
					if (in->header) kputs(in->header->target_name[dedupRec->tid] , &str);
					else kputw(dedupRec->tid, &str);
					kputc('\t', &str);
					kputw(dedupRec->pos + 1, &str); kputc('\t', &str);
					kputw(dedupRec->endpos, &str); kputc('\t', &str);
					kputsn(dedupRec->strand ? "-\t" : "+\t", 2, &str);
				}
				
				if (dedupRec->mtid < 0) kputsn("*\t0\t0\t+\t", 8, &str);
				else {
					if (in->header) kputs(in->header->target_name[dedupRec->mtid] , &str);
					else kputw(dedupRec->mtid, &str);
					kputc('\t', &str);
					kputw(dedupRec->mpos + 1, &str); kputc('\t', &str);
					kputw(dedupRec->mendpos, &str); kputc('\t', &str);
					kputsn(dedupRec->mstrand ? "-\t" : "+\t", 2, &str);
				}
				
				// we recompute the insert size just in case it is unavailable
				int32_t iSize = 0;
				if (dedupRec->tid == dedupRec->mtid && -1!=dedupRec->tid) iSize = dedupRec->mendpos - dedupRec->pos;
				kputw(iSize, &str); kputc('\t', &str);
				
				kputsn("*\t", 2, &str);
				
				kputw(dedupRec->swap, &str); //kputc('\t', &str);
				
				//TODO: for debugging
				kputc('\t', &str); kputul(dedupRec->lineno, &str);
				kputc('\t', &str); kputul(dedupRec->mlineno, &str);
				kputc('\t', &str); kputw(dedupRec->duplicated, &str);
				kputc('\n', &str);
				
				fprintf(stdout, "%s", str.s);
				free(str.s);
			}
			fprintf(stdout, "===\tEND : READ LIST\n");
			fflush(stdout);
			t_timepointIO = realtime();
			t_diffProcessing = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3) {
				fprintf(stderr, "[M::%s] Writing read %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nDedupKeys, t_diffProcessing, t_diffProcessing/60.0);
			}
			// END - for debugging purposes only, remove from production code
		}
		
		// we sort all the records internally
		t_timepointProcessing = realtime();
		if (opt->disponoParallel) {
			//qsort_mt(dedupKeys, nDedupKeys, sizeof(cpu_dedup_record_t), cmpPET, opt->n_threads);
			qsort(dedupKeys, nDedupKeys, sizeof(cpu_dedup_record_t), cmpPET);
		} else {
			qsort(dedupKeys, nDedupKeys, sizeof(cpu_dedup_record_t), cmpPET);
		}
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Sorted %u key records, processing %.2f sec ( %.2f min)..\n", __func__, nDedupKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// TODO: in furture we will use hybrid interal & external sorting with memory consideration
		// we write the sorting back into the same file and update both the sort state and the duplicate status
		// we next prepare a bit setting data represent unique entries as 1 and duplicate as 0 and persist it
		t_timepointProcessing = realtime();
		uint32_t nUniquePETs = 0; uint32_t nUnMappedPETs = 0;
		uint32_t *duplicates = markDuplicates(dedupKeys, nDedupKeys, n_processed, &nUniquePETs, &nUnMappedPETs);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		
		if (opt->calculateLibraryComplexity) {
			fprintf(stdout, "[I::%s] %u PETs, %u unique PETs (%.2f%%), %u unmapped PETs (%.2f%%)\n", __func__, nDedupKeys, nUniquePETs, nUniquePETs*100.0/nDedupKeys,nUnMappedPETs, nUnMappedPETs*100.0/nDedupKeys);
			if (bwa_verbose >= 3) {
				fprintf(stderr, "[I::%s] %u PETs, %u unique PETs (%.2f%%), %u unmapped PETs (%.2f%%)\n", __func__, nDedupKeys, nUniquePETs, nUniquePETs*100.0/nDedupKeys,nUnMappedPETs, nUnMappedPETs*100.0/nDedupKeys);
			}
		}
		
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Marked duplicates amongst %u keys %u records, processing %.2f sec (%.2f min)..\n", __func__, nDedupKeys, n_processed, t_diffProcessing, t_diffProcessing/60.0);
		}
		if (!duplicates) {
			fprintf(stderr, "[E::%s] out of memory to prepare duplicates bitstream for %u pairs processed..\n", __func__, n_pairProcessed);
		}
		
		// TOOD:future
		// writeDedupKeys(dedupKeyFilename.s, dedupKeys, nDedupKeys);

#if 0
		// TODO: for deubgging only
		//       check that those marked duplicated can be crossed checked as duplicated
		for(i=0; i<nDedupKeys; ++i) {
			cpu_dedup_record_t *dedupRec = &(dedupKeys[i]);
			if (0!=dedupRec->duplicated) {
				if (!is_line_duplicated(duplicates, dedupRec->lineno)) {
					if (bwa_verbose >= 1) {
						fprintf(stderr, "[E::%s] Dedup key %u pairid %lu lineno %u marked as duplicated, but not retrieved as such\n", __func__, i, dedupRec->pairId, dedupRec->lineno);
					}
				}
				if (!is_line_duplicated(duplicates, dedupRec->mlineno)) {
					if (bwa_verbose >= 1) {
						fprintf(stderr, "[E::%s] Dedup key %u pairid %lu mlineno %u marked as duplicated, but not retrieved as such\n", __func__, i, dedupRec->pairId, dedupRec->mlineno);
					}
				}
			}
		}
#endif
		
		if (CPU_DEDUP_DEBUG_SORTED_LIST==(CPU_DEDUP_DEBUG_SORTED_LIST&opt->outputSortedList)) {
			// for debugging purposes only, remove from production code
			t_timepointProcessing = realtime();
			fprintf(stdout, "\n\n===\tSORTED LIST\n");
			for(i=0; i<nDedupKeys; ++i) {
				
				cpu_dedup_record_t *dedupRec = &(dedupKeys[i]);
				
				kstring_t str;
				str.l = str.m = 0; str.s = 0;
				
				// pair id
				kputul(dedupRec->pairId, &str); kputc('\t', &str);
				
				// mapping quality
				kputw(dedupRec->qual, &str); kputc('\t', &str);
				
				// linker 1
				const char *text = shortLinkerTypeToStr(dedupRec->linkerType); kputs(text, &str); kputc('\t', &str);
				// linker 2
				text = shortLinkerTypeToStr(dedupRec->mlinkerType); kputs(text, &str); kputc('\t', &str);
				
				// final linker
				text = shortLinkerTypeToStr(dedupRec->finalLinkerType); kputs(text, &str); kputc('\t', &str);
				
				// classification
				text = shortOutputClassToStr(dedupRec->classification); kputs(text, &str); kputc('\t', &str);
				
				// read 1/tag 1
				kputw(dedupRec->rtid, &str); kputc('\t', &str);
				
				// read 2/tag 2
				kputw(dedupRec->mrtid, &str); kputc('\t', &str);
				
				if (dedupRec->tid < 0) kputsn("*\t0\t0\t+\t", 8, &str);
				else {
					if (in->header) kputs(in->header->target_name[dedupRec->tid] , &str);
					else kputw(dedupRec->tid, &str);
					kputc('\t', &str);
					kputw(dedupRec->pos + 1, &str); kputc('\t', &str);
					kputw(dedupRec->endpos, &str); kputc('\t', &str);
					kputsn(dedupRec->strand ? "-\t" : "+\t", 2, &str);
				}
				
				if (dedupRec->mtid < 0) kputsn("*\t0\t0\t+\t", 8, &str);
				else {
					if (in->header) kputs(in->header->target_name[dedupRec->mtid] , &str);
					else kputw(dedupRec->mtid, &str);
					kputc('\t', &str);
					kputw(dedupRec->mpos + 1, &str); kputc('\t', &str);
					kputw(dedupRec->mendpos, &str); kputc('\t', &str);
					kputsn(dedupRec->mstrand ? "-\t" : "+\t", 2, &str);
				}
				
				// we recompute the insert size just in case it is unavailable
				int32_t iSize = 0;
				if (dedupRec->tid == dedupRec->mtid && -1!=dedupRec->tid) iSize = dedupRec->mendpos - dedupRec->pos;
				kputw(iSize, &str); kputc('\t', &str);
				
				kputsn("*\t", 2, &str);
				
				kputw(dedupRec->swap, &str); //kputc('\t', &str);
				
				kputc('\t', &str); kputul(dedupRec->lineno, &str);
				kputc('\t', &str); kputul(dedupRec->mlineno, &str);
				kputc('\t', &str); kputw(dedupRec->duplicated, &str);
				kputc('\n', &str);
				
				fprintf(stdout, "%s", str.s);
				free(str.s);
			}
			fprintf(stdout, "===\tEND : SORTED LIST\n");
			fflush(stdout);
			
			t_timepointIO = realtime();
			t_diffProcessing = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3) {
				fprintf(stderr, "[M::%s] Writing sorted %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nDedupKeys, t_diffProcessing, t_diffProcessing/60.0);
			}
			// TODO: END - for debugging purposes only, remove from production code
		}
		
		// release the sorted records as we no longer need them
		free(dedupKeys);
		
		
		
		//-------------------------------------------------------------------------------------------------------------------
		
		
		// finally, we iterate the bam file and write to either the nr.bam or dup.bam based on the bit string
		// TOOD:future
		// writeNRandDuplicateBams(dedupKeyFilename.s, duplicates, nDedupKeys);
		
		// TODO: option 1: write out as unsorted by iterating the original bam file and use bits to check .nr / .dup
		//                 we then use samtools sort to sort these two output files
		// TODO: option 2: we can perform two sorting in the key sort phase so that we know the order of anchor1 and anchro2
		//                 however, it is unclear how we can try to get the incohesive bam records [need more thoughts]
		
		samclose(in);
		
		if (opt->toGroup) {
			if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open \"%s\" for reading.\n", __func__, argv[optind]);
				
				// TODO: what else do we need to clean up based on the above after looping
				free(fn_list);
				dedup_opt_terminate(opt);
				free(opt);
				
				return 1;
			}
			if (in->header == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to read the header from \"%s\".\n", __func__, argv[optind]);
				
				// TODO: what else do we need to clean up based on the above after looping
				free(fn_list);
				// close files, free and return
				samclose(in);
				
				dedup_opt_terminate(opt);
				free(opt);
				
				return 1;
			}
			
			// we reset the counters
			n_processed = 0;
			n_pairProcessed = 0;
			// TODO: we loop thru' the bam file and pair and write out the entries into .nr and .dup accordingly
			nNumRequested = opt->chunk_size * opt->n_threads;
			t_timepointProcessing = realtime();
			while ((bamRecs = readCPUBam(&n, nNumRequested, in, &bamsBuffer)) != 0) {
				t_timepointIO = realtime();
				t_diffIO = t_timepointIO - t_timepointProcessing;
				if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] read %d records..\n", __func__, n-bamsBuffer.unprocessed);
				
				// processing the pairing
				cpu_readid_t *readids = calloc(n, sizeof(cpu_readid_t));
				tag_group_t *tagGroups = calloc(n, sizeof(tag_group_t));
				cpu_dedup_record_t *dedupRecords = calloc(n, sizeof(cpu_dedup_record_t));
				
				uint32_t numTagGroups = 0;
				dedup_process_alns(opt, n_processed, n, bamRecs, readids, &numTagGroups, tagGroups, dedupRecords, &bamsBuffer, 1);
				t_timepointProcessing = realtime();
				t_diffProcessing = t_timepointProcessing - t_timepointIO;
				
				// TODO: write the bam results
				// TODO: for testing purposes
				// TODO: need to make sure that we have either 1 record (SE) or 2 records (PE)
				for(i=0; i<numTagGroups; ++i) {
					uint32_t lineno = n_processed+tagGroups[i].readIndex;
					// a tagGroup could have more than 1 tag, and they are collectively duplicate or unique
					uint32_t isDuplicated = 0;
					for(j=0, k=tagGroups[i].readIndex; j<tagGroups[i].n_reads; ++j, ++k, ++lineno) {
						if (0!=is_line_duplicated(duplicates, lineno)) {
							isDuplicated = 1; break;
						}
					}
					if (CPU_PAIR_OUTPUT_UU==tagGroups[i].outputClass && 0==isDuplicated) {
						samwrite(outfds[CPU_DEDUP_OUTPUT_REPRESENTATIVE], tagGroups[i].bamRecs[0]);
						samwrite(outfds[CPU_DEDUP_OUTPUT_REPRESENTATIVE], tagGroups[i].bamRecs[1]);
					} else {
						lineno = n_processed+tagGroups[i].readIndex;
						for(j=0, k=tagGroups[i].readIndex; j<tagGroups[i].n_reads; ++j, ++k, ++lineno) {
							if (0==isDuplicated) {
								samwrite(outfds[CPU_DEDUP_OUTPUT_REPRESENTATIVE], &bamRecs[k]);
							} else {
								bamRecs[k].core.flag |= BAM_FDUP;
								samwrite(outfds[CPU_DEDUP_OUTPUT_DUPLICATED], &bamRecs[k]);
							}
						}
					}
				}
				
				t_timepointIO = realtime();
				t_diffIO += (t_timepointIO - t_timepointProcessing);
				
				// clean up
				for(i=0; i<n; ++i) free(bamRecs[i].data);
				free(bamRecs);
				free(readids);
				free(dedupRecords);
				for(i=0; i<numTagGroups; ++i) {
					if (tagGroups[i].releaseBamRecs) {
						bam_destroy1(tagGroups[i].bamRecs[0]);
						bam_destroy1(tagGroups[i].bamRecs[1]);
						tagGroups[i].releaseBamRecs = 0;
					}
				}
				free(tagGroups);
				
				n_processed += (n-bamsBuffer.unprocessed);
				n_pairProcessed += numTagGroups;
				if (bwa_verbose >= 3) {
					t_diff = realtime() - G_t_real;
					fprintf(stderr, "[M::%s] %u tags %u pairs processed, %.0f tags/sec %.0f pairs/sec, %.2f min, i/o %.2f sec, processing %.2f sec..\n", __func__, n_processed, n_pairProcessed, 1.0*n_processed/t_diff, 1.0*n_pairProcessed/t_diff, t_diff/60.0, t_diffIO, t_diffProcessing);
				}
				
				t_timepointProcessing = realtime();
			}
			// close files, free and return
			// TODO: decide if we need additional closing due to re-opening
			samclose(in);
		}
		// END - PROCESSING
		
		free(duplicates);
		
		free(fn_list);
	}
	
	terminate_CPUBamBuffer(&bamsBuffer);
	
#if 0
	// TODO: report the dedup statistics
	int64_t totalInBins = 0;
	// TODO: count those which have been paired
	if (0==totalInBins) {
		if (bwa_verbose >= 2)
			fprintf(stderr, "[W::%s] %lld in bins of %u pairs (%lld tags). Has \"%s\" been paired?\n", __func__, totalInBins, n_pairProcessed, n_processed, argv[optind]);
	}
	fprintf(stdout, "##CPU\t%s\n", CPU_version);
	{
		fprintf(stdout, "##COMMAND\t");
		for (i = 0; i < argc; ++i)
			fprintf(stdout, " %s", argv[i]);
		fprintf(stdout, "\n");
	}
	
	long nTotal = n_pairProcessed;
	if (0==nTotal) nTotal = 1; // prevent division by zero!
	fprintf(stdout, ">>Library information\n");
	fprintf(stdout, "#Measure\tValue\n");
	fprintf(stdout, "Filename\t%s\n", argv[optind]);
	fprintf(stdout, ">>END\n");
	
	// END - report the dedup statistics
#endif
	
	if (opt->toGroup) {
		terminate_CPDedup_Outputs(outfds);
		// TODO: decide if we need additional closing due to re-opening
		// fclose(outfd);
	}
	
	free(dedupKeyFilename.s);
	
	dedup_opt_terminate(opt);
	free(opt);
	
	return ret;
}


