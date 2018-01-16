// TODO:
// 1. experimental code for clustering
// 2. consideration extension size, extension direction (both, 5'->3', 3'->5'
// 3. depending on the timing, we might include the tag count within?
//    e.g. it is fine to re-cluster for different filtering of the iPET filter?
//         but this can really be separated into another call unless we don't want to waste the loaded data
// 4. procedure:
//    a) sort on tag left and number those that overlap as the same cluster left id
//    b) sort on tag right and number those that overlap as the same cluster right id
//    c) sort on cluster left id, cluster right id, then coordinate...
//    d) output clusters with #iPET (with tag count for each anchor)
// 5. sample calls:
//    cpu cluster -t 8 -5 5,-20 -3 5,480 <lib>.FullLinker.NonChimeric.paired.UU.bam
// 6. optimization; there are still many rooms for optimization
//    a) distribute extend into I/O section
//    b) reduce the sort key comparisons
//    c) the processing of sort, process, sort, process, sort, process takes ~50sec, i.e. 20% of the total time
//    d)

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
#include "cpuspan.h"

extern double G_t_real;
extern const char *CPU_version;

#define G_EXTENSION_5P_SOURCE 5
#define G_EXTENSION_5P_SIZE   0
#define G_EXTENSION_3P_SOURCE 5
#define G_EXTENSION_3P_SIZE   500

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

// TODO: optimization: we might not need lineno and mlineno if we do not need bam access
// TODO: swap is being overloaded as 're-check'
typedef struct {
	int32_t tid;
	int32_t pos;
	int32_t endpos;
	int32_t mtid;
	int32_t mpos;
	int32_t mendpos;
	// TODO: we downgraded qual from 8 bits to 7 bits to keep within 32 bits
	//uint32_t metaData;
	uint32_t strand:1, mstrand:1, classification:3, finalLinkerType:4, linkerType:4, mlinkerType:4, qual:7, uniq:1, muniq:1, rtid:2, mrtid:2, usable:1, swap:1;
	// might have to change to uint64_t if supporting more than 4 billion PET tags
	uint32_t lineno; // it is the bam "record id"; in future, we might use this for "direct" bam record access
	uint32_t mlineno; // it is the bam "record id"; in future, we might use this for "direct" bam record access
	int32_t cid; // cluster id
	int32_t mcid; // mate cluster id
	int32_t iPET; // mate cluster id
	int32_t oid; // this is the order id for tag1
	int32_t moid; // this is the order id for tag2
	long pairId; // TODO: for checking purpose only
} cpu_cluster_record_t;

typedef struct {
	int32_t pos;
	int32_t endpos;
	int32_t mpos;
	int32_t mendpos;
	int32_t cid;
	int32_t mcid;
	int32_t iPET;
	int32_t oid;
	int32_t moid;
	int32_t offset; // offset within a bin
	uint32_t lineno; // it is the bam "record id"; in future, we might use this for "direct" bam record access
	long pairId; // TODO: for checking purpose only
} cpu_subcluster_record_t;


typedef struct {
	int source;
	int size;
} extension_opt_t;

typedef struct {
	int chunk_size;
	int n_threads;

	int selfLigation;
	extension_opt_t extension5p;
	extension_opt_t extension3p;
	
	int disponoParallel; // TODO: we should have bit based flags!
	int outputSortedList;
	int cycle; // TODO: we should have bit based flags!
	
	kstring_t outputPrefix;
} cluster_opt_t;

cluster_opt_t *cluster_opt_init()
{
	cluster_opt_t *o;
	o = calloc(1, sizeof(cluster_opt_t));
	
	o->chunk_size = 1000000;
	o->n_threads = 1;
	
	o->selfLigation = G_SELF_LIGATION;
	o->extension5p.source = G_EXTENSION_5P_SOURCE;
	o->extension5p.size = G_EXTENSION_5P_SIZE;
	o->extension3p.source = G_EXTENSION_3P_SOURCE;
	o->extension3p.size = G_EXTENSION_3P_SIZE;
	
	o->disponoParallel = 0;
	o->outputSortedList = 0;
	o->cycle = -1;
	
	memset(&(o->outputPrefix), 0, sizeof(kstring_t));
	return o;
}

void cluster_opt_terminate(cluster_opt_t *o)
{
	free(o->outputPrefix.s);
}

typedef struct {
	const cluster_opt_t *opt;
	
	uint32_t n_processed;
	
	bam1_t *bamRecs;
	cpu_readid_t *readids;
	
	int nTagGroups;
	tag_group_t *tagGroups;
	cpu_cluster_record_t *clusterRecords;

	int32_t nTargets;
	uint32_t *targetLens;
	
	uint32_t nUpdateBam;
} worker_t;

static void readid_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	// read id parsing does not care about if reads are paried-end
	int nRet = parseReadId(bam1_qname(&(w->bamRecs[i])), &(w->readids[i]), (BAM_FREAD2==(BAM_FREAD2 & w->bamRecs[i].core.flag))?1:0);
	if (0!=nRet) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to parse read id \"%s\". Make sure that read is prepared by CPU.\n", __func__, bam1_qname(&(w->bamRecs[i])));
	}
}

static void cluster_pair_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	// TODO: set up the result in w->clusterRecords[i]
	cpu_cluster_record_t *clusterRec = &(w->clusterRecords[i]);
	// we considered this PET as duplicated by default and non-uniquely mapped
	// this has the same effect as excluding this PET from all downstream analysis
	// thus, we only 'reset' those which are usable for downstream analysis
	clusterRec->usable = 0; clusterRec->uniq = 0; clusterRec->muniq = 0;
	clusterRec->iPET = 0;
	
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
					anchor1ReadIndex = readTags[0][0]; clusterRec->endpos = readEnd[0][0];
				} else {
					anchor1ReadIndex = readTags[1][1]; clusterRec->endpos = readEnd[1][1];
				}
			} else {
				// different chromosome! thus consider discordant
				// TODO: but we still pick one for processing
				anchorDiscordant |= 0x01;
				if (readMappedLen[0][0]>=readMappedLen[1][1]) {
					anchor1ReadIndex = readTags[0][0]; clusterRec->endpos = readEnd[0][0];
				} else {
					anchor1ReadIndex = readTags[1][1]; clusterRec->endpos = readEnd[1][1];
				}
			}
		} else {
			// only R/1 tag 1 unique
			anchor1ReadIndex = readTags[0][0]; clusterRec->endpos = readEnd[0][0];
		}
		anchor1MapState = CPU_MAPSTATE_UNIQUE;
	} else if (CPU_MAPSTATE_UNIQUE==readTagsMapState[1][1]) {
		// only R/2 tag 2 unique
		anchor1ReadIndex = readTags[1][1]; clusterRec->endpos = readEnd[1][1];
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
					anchor1ReadIndex = readTags[1][1]; clusterRec->endpos = readEnd[1][1];
				}
			} else {
				if (CPU_MAPSTATE_NA==readTagsMapState[1][1]) {
					anchor1ReadIndex = readTags[0][0]; clusterRec->endpos = readEnd[0][0];
				} else {
					if (w->bamRecs[readTags[0][0]].core.l_qseq>=w->bamRecs[readTags[1][1]].core.l_qseq) {
						anchor1ReadIndex = readTags[0][0]; clusterRec->endpos = readEnd[0][0];
					} else {
						anchor1ReadIndex = readTags[1][1]; clusterRec->endpos = readEnd[1][1];
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
					anchor2ReadIndex = readTags[1][0]; clusterRec->mendpos = readEnd[1][0];
				} else {
					anchor2ReadIndex = readTags[0][1]; clusterRec->mendpos = readEnd[0][1];
				}
			} else {
				// different chromosome!
				// TODO: but we still pick one for processing
				anchorDiscordant |= 0x02;
				if (readMappedLen[1][0]>=readMappedLen[0][1]) {
					anchor2ReadIndex = readTags[1][0]; clusterRec->mendpos = readEnd[1][0];
				} else {
					anchor2ReadIndex = readTags[0][1]; clusterRec->mendpos = readEnd[0][1];
				}
			}
		} else {
			// only R/2 tag 1 unique
			anchor2ReadIndex = readTags[1][0]; clusterRec->mendpos = readEnd[1][0];
		}
		anchor2MapState = CPU_MAPSTATE_UNIQUE;
	} else if (CPU_MAPSTATE_UNIQUE==readTagsMapState[0][1]) {
		// only R/1 tag 2 unique
		anchor2ReadIndex = readTags[0][1]; clusterRec->mendpos = readEnd[0][1];
		anchor2MapState = CPU_MAPSTATE_UNIQUE;
	} else {
		// other cases
		if (CPU_MAPSTATE_REPEAT==readTagsMapState[1][0] || CPU_MAPSTATE_REPEAT==readTagsMapState[0][1]) {
			// either repeat, considered as repeat
			anchor2MapState = CPU_MAPSTATE_REPEAT;
			// we pick the read tag with longer mapped region
			if (readMappedLen[1][0]>=readMappedLen[0][1]) {
				anchor2ReadIndex = readTags[1][0]; clusterRec->mendpos = readEnd[1][0];
			} else {
				anchor2ReadIndex = readTags[0][1]; clusterRec->mendpos = readEnd[0][1];
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
					anchor2ReadIndex = readTags[0][1]; clusterRec->mendpos = readEnd[0][1];
				}
			} else {
				if (CPU_MAPSTATE_NA==readTagsMapState[0][1]) {
					anchor2ReadIndex = readTags[1][0]; clusterRec->mendpos = readEnd[1][0];
				} else {
					if (w->bamRecs[readTags[1][0]].core.l_qseq>=w->bamRecs[readTags[0][1]].core.l_qseq) {
						anchor2ReadIndex = readTags[1][0]; clusterRec->mendpos = readEnd[1][0];
					} else {
						anchor2ReadIndex = readTags[0][1]; clusterRec->mendpos = readEnd[0][1];
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
			clusterRec->swap = 1;
			anchor1ReadIndex = anchor2ReadIndex;
			anchor2ReadIndex = -1;
			clusterRec->endpos = clusterRec->mendpos;
			clusterRec->mendpos = 0;
			
			uint8_t anchorMapState = anchor1MapState; anchor1MapState = anchor2MapState; anchor2MapState = anchorMapState;
			
			c1 = &(w->bamRecs[anchor1ReadIndex].core);
			ri1 = &(w->readids[anchor1ReadIndex]);
			
			// adjust the clusterlication record
			clusterRec->classification = ri1->classification & 0x07;
			clusterRec->finalLinkerType = ri1->finalLinkerType & 0x0F;
			clusterRec->qual = (c1->qual>0x7F) ? 0x7F : c1->qual;
			clusterRec->pairId = ri1->pairId;
			
			clusterRec->tid = c1->tid; clusterRec->pos = c1->pos;
			clusterRec->strand = (BAM_FREVERSE==(c1->flag&BAM_FREVERSE)) ? 1 : 0;
			clusterRec->linkerType = ri1->linker_type & 0x0F;
			clusterRec->uniq = (CPU_MAPSTATE_UNIQUE==anchor1MapState) ? 1 : 0;
			clusterRec->rtid = 2*ri1->readId+ri1->tagId+1;
			clusterRec->lineno = w->n_processed + anchor1ReadIndex; // TODO: 32 bits or 64 bits?
			
			// there is no need to set the value if we want default zeros
			clusterRec->mtid = -1;
			/*
			 clusterRec->mtid = 0; clusterRec->pos = 0;
			 clusterRec->mstrand = 0;
			 clusterRec->mlinkerType = 0;
			 clusterRec->muniq = 0;
			 clusterRec->mrtid = 0;
			 */
			
			w->tagGroups[i].bamRecs[0] = &(w->bamRecs[anchor1ReadIndex]);
		}
	} else {
		if (-1==anchor2ReadIndex) {
			// anchor1 avail, anchor2 NA
			// no swapping needed, proceed as-is
			c1 = &(w->bamRecs[anchor1ReadIndex].core);
			ri1 = &(w->readids[anchor1ReadIndex]);
			
			// adjust the clusterlication record
			clusterRec->classification = ri1->classification & 0x07;
			clusterRec->finalLinkerType = ri1->finalLinkerType & 0x0F;
			clusterRec->qual = (c1->qual>0x7F) ? 0x7F : c1->qual;
			clusterRec->pairId = ri1->pairId;
			
			clusterRec->tid = c1->tid; clusterRec->pos = c1->pos;
			clusterRec->strand = (BAM_FREVERSE==(c1->flag&BAM_FREVERSE)) ? 1 : 0;
			clusterRec->linkerType = ri1->linker_type & 0x0F;
			clusterRec->uniq = (CPU_MAPSTATE_UNIQUE==anchor1MapState) ? 1 : 0;
			clusterRec->rtid = 2*ri1->readId+ri1->tagId+1;
			clusterRec->lineno = w->n_processed + anchor1ReadIndex; // TODO: 32 bits or 64 bits?
			
			// there is no need to set the value if we want default zeros
			clusterRec->mtid = -1;
			/*
			clusterRec->mtid = 0; clusterRec->pos = 0;
			clusterRec->mstrand = 0;
			clusterRec->mlinkerType = 0;
			clusterRec->muniq = 0;
			clusterRec->mrtid = 0;
			 */
			
			w->tagGroups[i].bamRecs[0] = &(w->bamRecs[anchor1ReadIndex]);
			
		} else {
			// anchor1 avail, anchor2 avail
			c1 = &(w->bamRecs[anchor1ReadIndex].core);
			ri1 = &(w->readids[anchor1ReadIndex]);
			c2 = &(w->bamRecs[anchor2ReadIndex].core);
			ri2 = &(w->readids[anchor2ReadIndex]);
			
			if (c2->tid<c1->tid) clusterRec->swap = 1;
			else if (c2->tid==c1->tid) {
				if (c2->pos<c1->pos) clusterRec->swap = 1;
				else if (c2->pos==c1->pos) {
					if (clusterRec->mendpos<clusterRec->endpos) clusterRec->swap = 1;
				}
			}

			if (clusterRec->swap) {
				const bam1_core_t *c = c1; c1 = c2; c2 = c;
				const cpu_readid_t *ri = ri1; ri1 = ri2; ri2 = ri;
				int32_t endpos = clusterRec->endpos; clusterRec->endpos = clusterRec->mendpos; clusterRec->mendpos = endpos;
				clusterRec->lineno = w->n_processed + anchor2ReadIndex; // TODO: 32 bits or 64 bits?
				clusterRec->mlineno = w->n_processed + anchor1ReadIndex; // TODO: 32 bits or 64 bits?
				
				uint8_t anchorMapState = anchor1MapState; anchor1MapState = anchor2MapState; anchor2MapState = anchorMapState;
			} else {
				clusterRec->lineno = w->n_processed + anchor1ReadIndex; // TODO: 32 bits or 64 bits?
				clusterRec->mlineno = w->n_processed + anchor2ReadIndex; // TODO: 32 bits or 64 bits?
			}
			
			// adjust the clusterlication record
			clusterRec->classification = ri1->classification & 0x07;
			clusterRec->finalLinkerType = ri1->finalLinkerType & 0x0F;
			// we keep the lower of the two scores so that we do not over-reported in the case of repeats
			uint8_t qual = (c1->qual <= c2->qual) ? c1->qual : c2->qual;
			clusterRec->qual = (qual>0x7F) ? 0x7F : qual;
			clusterRec->pairId = ri1->pairId;
			
			clusterRec->tid = c1->tid; clusterRec->pos = c1->pos;
			clusterRec->strand = (BAM_FREVERSE==(c1->flag&BAM_FREVERSE)) ? 1 : 0;
			clusterRec->linkerType = ri1->linker_type & 0x0F;
			clusterRec->uniq = (CPU_MAPSTATE_UNIQUE==anchor1MapState) ? 1 : 0;
			clusterRec->rtid = 2*ri1->readId+ri1->tagId+1;
			
			clusterRec->mtid = c2->tid; clusterRec->mpos = c2->pos;
			clusterRec->mstrand = (BAM_FREVERSE==(c2->flag&BAM_FREVERSE)) ? 1 : 0;
			clusterRec->mlinkerType = ri2->linker_type & 0x0F;
			clusterRec->muniq = (CPU_MAPSTATE_UNIQUE==anchor2MapState) ? 1 : 0;
			clusterRec->mrtid = 2*ri2->readId+ri2->tagId+1;

			w->tagGroups[i].bamRecs[0] = &(w->bamRecs[anchor1ReadIndex]);
			w->tagGroups[i].bamRecs[1] = &(w->bamRecs[anchor2ReadIndex]);
			
			// check if this is not a self-ligation
			if (clusterRec->tid==clusterRec->mtid) {
#if 1
				// for debugging purpose, guoliang version
				{
					int32_t center = (clusterRec->pos + clusterRec->endpos) / 2;
					clusterRec->pos = center;
					clusterRec->endpos = center;
					center = (clusterRec->mpos + clusterRec->mendpos) / 2;
					clusterRec->mpos = center;
					clusterRec->mendpos = center;
				}
#endif
				int32_t tlen = getIntraSpan(clusterRec->pos,clusterRec->endpos,clusterRec->mpos,clusterRec->mendpos);
				if (tlen>w->opt->selfLigation) { clusterRec->usable = 1; }
			} else {
				// trans-interaction
				clusterRec->usable = 1;
			}
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

void cluster_process_alns(const cluster_opt_t *opt, uint32_t n_processed, int n, bam1_t *bamRecs, cpu_readid_t *readids, uint32_t *pNumTagGroups, tag_group_t *tagGroups, cpu_cluster_record_t *clusterRecords, CPUBamBuffer_t* bamsBuffer, uint32_t nUpdateBam)
{
	int i;
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	if (n<=0) return;
	
	w.opt = opt;
	w.n_processed = n_processed;
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
	w.clusterRecords = clusterRecords;
	kt_for(opt->n_threads, cluster_pair_worker, &w, numTagGroups);
	
	// TODO: accumulate all the span information
	for(i=0; i<opt->n_threads; ++i) {
	}
}

int updateClusterTotalRecords(const char *clusterKeyFilename, uint32_t n_pairProcessed) {
	const char *mode = "r+b";
	FILE *fd = fopen(clusterKeyFilename, mode);
	if (fd == NULL) {
		fprintf(stderr, "[E::%s] fail to open '%s' for update: %s\n", __func__, clusterKeyFilename, strerror (errno));
		return errno;
	}
	
	// TODO: header processing synchronization
	fseek(fd, 3*sizeof(uint32_t), SEEK_SET);
	err_fwrite(&n_pairProcessed, sizeof(n_pairProcessed), 1, fd);// filler for total number of records
	fclose(fd);
	
	return 0;
}

// TODO:
cpu_cluster_record_t *readClusterKeys(const char *clusterKeyFilename, uint32_t *nClusterKeys) {
	*nClusterKeys = 0;
	cpu_cluster_record_t *clusterKeys = NULL;
	
	// open file and decide on amount of storage to allocate
	const char *mode = "rb";
	FILE *fd = fopen(clusterKeyFilename, mode);
	if (fd == NULL) {
		fprintf(stderr, "[E::%s] fail to open '%s' for update: %s\n", __func__, clusterKeyFilename, strerror (errno));
		return NULL;
	}
	
	uint32_t nValue = 0;
	err_fread_noeof(&nValue, sizeof(nValue), 1, fd); // cpu version
	err_fread_noeof(&nValue, sizeof(nValue), 1, fd); // sort state: unsorted=0, chunk sorted=1, fully sorted=2
	err_fread_noeof(&nValue, sizeof(nValue), 1, fd); // record size
	err_fread_noeof(&nValue, sizeof(nValue), 1, fd); // filler for total number of records
	*nClusterKeys = nValue;
	if (*nClusterKeys) {
		clusterKeys = (cpu_cluster_record_t *) calloc(*nClusterKeys, sizeof(cpu_cluster_record_t));
		if (clusterKeys) {
			// chunk read the data
			
			// TODO: more error checking needed
			//       what if it is truncated? what if there are more records than allocated?
			uint32_t nRead = 0;
			while (nRead<*nClusterKeys) {
				uint32_t nChunkRecs = 0;
				err_fread_noeof(&nChunkRecs, sizeof(nChunkRecs), 1, fd); // number of records in chunk
				uint32_t nSorted = 0;
				err_fread_noeof(&nSorted, sizeof(nSorted), 1, fd); // chunk sorted?
				if ((nRead+nChunkRecs)>*nClusterKeys) {
					fprintf(stderr, "[E::%s] There are more records than recorded. Specified: %u, Reading: %u.\n", __func__, *nClusterKeys, nRead+nChunkRecs);
					*nClusterKeys = nRead;
					break;
				}
				err_fread_noeof(&(clusterKeys[nRead]), sizeof(cpu_cluster_record_t), nChunkRecs, fd);
				nRead += nChunkRecs;
			}
		} else {
			fprintf(stderr, "[E::%s] out of memory to read %u records from '%s'e: %s\n", __func__, *nClusterKeys, clusterKeyFilename, strerror (errno));
		}
	}
	
	// close file and return prep'd data
	fclose(fd);
	return clusterKeys;
}

void report_iPETs_List (const cpu_cluster_record_t *clusterKeys, uint32_t nClusterKeys, const bam_header_t *header, int printSwap) {
	int32_t i;
	fprintf(stdout, "#pairId\tqual\tusable\tcid\tmcid\toid\tmoid\tiPETs\tlinker\tmlinker\tflinker\tclass\ttag1\ttag2\ttid\tpos\tendpos\tstrand\tmtid\tmpos\tmendpos\tmstrand\tisize");
#if 0
	fprintf(stdout, "\tfield");
#endif
	if (0!=printSwap) fprintf(stdout, "\tswap");
#if 0
	fprintf(stdout, "\tlineno\tmlineno");
#endif
	fprintf(stdout, "\n");
	
	for(i=0; i<nClusterKeys; ++i) {
		
		const cpu_cluster_record_t *clusterRec = &(clusterKeys[i]);
		
		kstring_t str;
		str.l = str.m = 0; str.s = 0;
		
		// pair id
		kputul(clusterRec->pairId, &str);
		
		// mapping quality
		kputc('\t', &str); kputw(clusterRec->qual, &str);
		
		kputc('\t', &str); kputl(clusterRec->usable, &str);
		kputc('\t', &str); kputl(clusterRec->cid, &str);
		kputc('\t', &str); kputl(clusterRec->mcid, &str);
		kputc('\t', &str); kputl(clusterRec->oid, &str);
		kputc('\t', &str); kputl(clusterRec->moid, &str);
		kputc('\t', &str); kputl(clusterRec->iPET, &str);
		
		// linker 1
		kputc('\t', &str);
		const char *text = shortLinkerTypeToStr(clusterRec->linkerType); kputs(text, &str);
		// linker 2
		kputc('\t', &str); text = shortLinkerTypeToStr(clusterRec->mlinkerType); kputs(text, &str);
		
		// final linker
		kputc('\t', &str); text = shortLinkerTypeToStr(clusterRec->finalLinkerType); kputs(text, &str);
		
		// classification
		kputc('\t', &str); text = shortOutputClassToStr(clusterRec->classification); kputs(text, &str);
		
		// read 1/tag 1
		kputc('\t', &str); kputw(clusterRec->rtid, &str);
		
		// read 2/tag 2
		kputc('\t', &str); kputw(clusterRec->mrtid, &str);
		
		kputc('\t', &str);
		if (clusterRec->tid < 0) kputsn("*\t0\t0\t+", 7, &str);
		else {
			if (header) kputs(header->target_name[clusterRec->tid] , &str);
			else kputw(clusterRec->tid, &str);
			kputc('\t', &str); kputw(clusterRec->pos, &str);
			kputc('\t', &str); kputw(clusterRec->endpos, &str);
			kputc('\t', &str); kputc(clusterRec->strand ? '-' : '+', &str);
		}
		
		kputc('\t', &str);
		if (clusterRec->mtid < 0) kputsn("*\t0\t0\t+", 7, &str);
		else {
			if (header) kputs(header->target_name[clusterRec->mtid] , &str);
			else kputw(clusterRec->mtid, &str);
			kputc('\t', &str); kputw(clusterRec->mpos, &str);
			kputc('\t', &str); kputw(clusterRec->mendpos, &str);
			kputc('\t', &str); kputc(clusterRec->mstrand ? '-' : '+', &str);
		}
		
		// we recompute the insert size just in case it is unavailable
		int32_t iSize = 0;
		if (clusterRec->tid == clusterRec->mtid && -1!=clusterRec->tid) iSize = clusterRec->mendpos - clusterRec->pos;
		kputc('\t', &str);kputw(iSize, &str);
		
#if 0
		kputc('\t', &str);kputc('*', &str);
#endif
		if (0!=printSwap) {
			kputc('\t', &str); kputw(clusterRec->swap, &str);
		}
		
#if 0
		kputc('\t', &str); kputul(clusterRec->lineno, &str);
		kputc('\t', &str); kputul(clusterRec->mlineno, &str);
#endif
		kputc('\n', &str);
		
		fprintf(stdout, "%s", str.s);
		free(str.s);
	}
}

int cmpPETTag1(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int32_t i32comp = ia->usable - ib->usable; if (0!=i32comp) return i32comp;
#if 0
	// sort by chrom1, start1, end1,
	//int32_t
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
#else
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
#endif
	// strand1, strand2, final linker, linker1, linker2, score
	i32comp = ia->strand - ib->strand; if (0!=i32comp) return i32comp;
	i32comp = ia->classification - ib->classification; if (0!=i32comp) return i32comp;
	i32comp = ia->finalLinkerType - ib->finalLinkerType; if (0!=i32comp) return i32comp;
	i32comp = ia->linkerType - ib->linkerType; if (0!=i32comp) return i32comp;
	// TODO: qual is higher first!
	//i32comp = ia->qual - ib->qual; if (0!=i32comp) return i32comp;
	i32comp = ib->qual - ia->qual; if (0!=i32comp) return i32comp;
	i32comp = ia->lineno - ib->lineno; //if (!i32comp) return i32comp;
	return i32comp;
}

int cmpPETTag2(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int32_t i32comp = ia->usable - ib->usable; if (0!=i32comp) return i32comp;
#if 0
	// follow by chrom2, start2, end2
	//int32_t
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	// IMPT: additional; experimental
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
#else
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
#endif
	
	// strand1, strand2, final linker, linker1, linker2, score
	// TODO: there are 8x2 comparison here which can be simplified to 4 or 2
	i32comp = ia->mstrand - ib->mstrand; if (0!=i32comp) return i32comp;
	i32comp = ia->classification - ib->classification; if (0!=i32comp) return i32comp;
	i32comp = ia->finalLinkerType - ib->finalLinkerType; if (0!=i32comp) return i32comp;
	i32comp = ia->mlinkerType - ib->mlinkerType; if (0!=i32comp) return i32comp;
	// TODO: qual is higher first!
	//i32comp = ia->qual - ib->qual; if (0!=i32comp) return i32comp;
	i32comp = ib->qual - ia->qual; if (0!=i32comp) return i32comp;
	i32comp = ia->lineno - ib->lineno; //if (!i32comp) return i32comp;
	return i32comp;
}

int cmpPETBothTags(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int32_t i32comp = ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	// sort by chrom1, chrom2
	//int32_t
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	
	// sort by start1, end1,
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	// follow by start2, end2
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;

	// TODO: do we need all these comparisons? maybe only the tie-breaker?
	
	// strand1, strand2, final linker, linker1, linker2, score
	// TODO: there are 8x2 comparison here which can be simplified to 4 or 2
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
}

int cmpPETTag1Order(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
#if 0
	int32_t i32comp = ia->oid - ib->oid;
#else
	int32_t i32comp = ia->usable - ib->usable; if (0!=i32comp) return i32comp;
	// sort by chrom1, start1, end1,
#if 0
	//int32_t
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
#else
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
#endif
	// strand1, strand2, final linker, linker1, linker2, score
	i32comp = ia->strand - ib->strand; if (0!=i32comp) return i32comp;
	i32comp = ia->classification - ib->classification; if (0!=i32comp) return i32comp;
	i32comp = ia->finalLinkerType - ib->finalLinkerType; if (0!=i32comp) return i32comp;
	i32comp = ia->linkerType - ib->linkerType; if (0!=i32comp) return i32comp;
	// TODO: qual is higher first!
	//i32comp = ia->qual - ib->qual; if (0!=i32comp) return i32comp;
	i32comp = ib->qual - ia->qual; if (0!=i32comp) return i32comp;
	i32comp = ia->lineno - ib->lineno; //if (!i32comp) return i32comp;
#endif
	return i32comp;
}

int cmpPETTag2Order(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
#if 0
	int32_t i32comp = ia->moid - ib->moid;
#else
	int32_t i32comp = ia->usable - ib->usable; if (0!=i32comp) return i32comp;
	// follow by chrom2, start2, end2
#if 0
	//int32_t
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	// IMPT: additional; experimental
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
#else
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
#endif
	
	// strand1, strand2, final linker, linker1, linker2, score
	// TODO: there are 8x2 comparison here which can be simplified to 4 or 2
	i32comp = ia->mstrand - ib->mstrand; if (0!=i32comp) return i32comp;
	i32comp = ia->classification - ib->classification; if (0!=i32comp) return i32comp;
	i32comp = ia->finalLinkerType - ib->finalLinkerType; if (0!=i32comp) return i32comp;
	i32comp = ia->mlinkerType - ib->mlinkerType; if (0!=i32comp) return i32comp;
	// TODO: qual is higher first!
	//i32comp = ia->qual - ib->qual; if (0!=i32comp) return i32comp;
	i32comp = ib->qual - ia->qual; if (0!=i32comp) return i32comp;
	i32comp = ia->lineno - ib->lineno; //if (!i32comp) return i32comp;
#endif
	return i32comp;
}

int cmpPETBothTagsOrder(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int32_t i32comp;
	//= ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	// sort by chrom1, chrom2
	//int32_t
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	
	// TODO: any proxy for optimization?
	
	// sort by start1, end1,
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	// follow by start2, end2
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	// TODO: do we need all these comparisons? maybe only the tie-breaker?
	
	// strand1, strand2, final linker, linker1, linker2, score
	// TODO: there are 8x2 comparison here which can be simplified to 4 or 2
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
}

// TODO: this routine can be further optimized
static void tags_extension_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	cpu_cluster_record_t *clusterRecord = &(w->clusterRecords[i]);

	if (-1!=clusterRecord->tid) {
		int32_t pos5p;
		int32_t pos3p;
		
#if 1
		// TODO: this for checking against guoliang's program output
		{
			int32_t center = (clusterRecord->pos + clusterRecord->endpos) / 2;
			clusterRecord->pos = center;
			clusterRecord->endpos = center;
		}
#endif
		
		if (5==w->opt->extension5p.source) {
			// use 5' position as reference
			if (clusterRecord->strand) {
				pos5p = clusterRecord->endpos - w->opt->extension5p.size;
			} else {
				pos5p = clusterRecord->pos + w->opt->extension5p.size;
			}
		} else {
			// use 3' position as reference
			if (clusterRecord->strand) {
				pos5p = clusterRecord->pos - w->opt->extension5p.size;
			} else {
				pos5p = clusterRecord->endpos + w->opt->extension5p.size;
			}
		}
		
		if (5==w->opt->extension3p.source) {
			// use 5' position as reference
			if (clusterRecord->strand) {
				pos3p = clusterRecord->endpos - w->opt->extension3p.size;
			} else {
				pos3p = clusterRecord->pos + w->opt->extension3p.size;
			}
		} else {
			// use 3' position as reference
			if (clusterRecord->strand) {
				pos3p = clusterRecord->pos - w->opt->extension3p.size;
			} else {
				pos3p = clusterRecord->endpos + w->opt->extension3p.size;
			}
		}
		
		// check that we are within the chromatin length!!!
		if (pos5p<0) { pos5p = 0; }
		else if (pos5p>w->targetLens[clusterRecord->tid]) { pos5p = w->targetLens[clusterRecord->tid]; }
		if (pos3p<0) { pos3p = 0; }
		else if (pos3p>w->targetLens[clusterRecord->tid]) { pos3p = w->targetLens[clusterRecord->tid]; }

		if (clusterRecord->strand) {
			w->clusterRecords[i].pos = pos3p;
			w->clusterRecords[i].endpos = pos5p;
		} else {
			w->clusterRecords[i].pos = pos5p;
			w->clusterRecords[i].endpos = pos3p;
		}
	}

	if (-1!=clusterRecord->mtid) {
		int32_t pos5p;
		int32_t pos3p;
		
#if 1
		// TODO: this for checking against guoliang's program output
		{
			int32_t center = (clusterRecord->mpos + clusterRecord->mendpos) / 2;
			clusterRecord->mpos = center;
			clusterRecord->mendpos = center;
		}
#endif
		
		if (5==w->opt->extension5p.source) {
			// use 5' position as reference
			if (clusterRecord->mstrand) {
				pos5p = clusterRecord->mendpos - w->opt->extension5p.size;
			} else {
				pos5p = clusterRecord->mpos + w->opt->extension5p.size;
			}
		} else {
			// use 3' position as reference
			if (clusterRecord->mstrand) {
				pos5p = clusterRecord->mpos - w->opt->extension5p.size;
			} else {
				pos5p = clusterRecord->mendpos + w->opt->extension5p.size;
			}
		}
		
		if (5==w->opt->extension3p.source) {
			// use 5' position as reference
			if (clusterRecord->mstrand) {
				pos3p = clusterRecord->mendpos - w->opt->extension3p.size;
			} else {
				pos3p = clusterRecord->mpos + w->opt->extension3p.size;
			}
		} else {
			// use 3' position as reference
			if (clusterRecord->mstrand) {
				pos3p = clusterRecord->mpos - w->opt->extension3p.size;
			} else {
				pos3p = clusterRecord->mendpos + w->opt->extension3p.size;
			}
		}
		
		// check that we are within the chromatin length!!!
		if (pos5p<0) { pos5p = 0; }
		else if (pos5p>w->targetLens[clusterRecord->mtid]) { pos5p = w->targetLens[clusterRecord->mtid]; }
		if (pos3p<0) { pos3p = 0; }
		else if (pos3p>w->targetLens[clusterRecord->mtid]) { pos3p = w->targetLens[clusterRecord->mtid]; }
		
		if (clusterRecord->mstrand) {
			w->clusterRecords[i].mpos = pos3p;
			w->clusterRecords[i].mendpos = pos5p;
		} else {
			w->clusterRecords[i].mpos = pos5p;
			w->clusterRecords[i].mendpos = pos3p;
		}
	}
}

void extend_PETs(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys, uint32_t *targetLens, int32_t nTargets) {
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	w.opt = opt;
	w.clusterRecords = clusterRecords;
	w.nTargets = nTargets;
	w.targetLens = targetLens;
	kt_for(opt->n_threads, tags_extension_worker, &w, nClusterKeys);
}

static void PETag1_cluster_step1_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	cpu_cluster_record_t *clusterRecord = &(w->clusterRecords[i]);
	clusterRecord->oid = i + 1;
	if (-1==clusterRecord->tid || 0==clusterRecord->usable) {
		clusterRecord->cid = -1;
	} else {
		
		cpu_cluster_record_t *pClusterRecord = clusterRecord - 1;
		if (pClusterRecord->tid==clusterRecord->tid && is_overlap(pClusterRecord->pos, pClusterRecord->endpos, clusterRecord->pos, clusterRecord->endpos)) {
			clusterRecord->cid = 1;
		} /* else { // default 0 }*/
	}
}

void cluster_PETags1(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	int32_t i;
	worker_t w;
	int32_t cid = 0; // cluster id
	
	w.opt = opt;
	w.clusterRecords = clusterRecords+1;
	clusterRecords[0].oid = 0;
	kt_for(opt->n_threads, PETag1_cluster_step1_worker, &w, nClusterKeys-1);
	
	// TODO: use serial approach to set up the correct cid
	for(i=0; i<nClusterKeys; ++i) {
		if (0==clusterRecords[i].cid) {
#if 0
			cid = i;
#else
			cid = clusterRecords[i].oid;
#endif
			clusterRecords[i].cid = cid;
		} else if (1==clusterRecords[i].cid) {
			clusterRecords[i].cid = cid;
		}
	}
}

static void PETag2_cluster_step1_worker(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	cpu_cluster_record_t *clusterRecord = &(w->clusterRecords[i]);
	clusterRecord->moid = i + 1;
	if (-1==clusterRecord->mtid || 0==clusterRecord->usable) {
		clusterRecord->mcid = -1;
	} else {
		cpu_cluster_record_t *pClusterRecord = clusterRecord - 1;
		if (pClusterRecord->mtid==clusterRecord->mtid && is_overlap(pClusterRecord->mpos, pClusterRecord->mendpos, clusterRecord->mpos, clusterRecord->mendpos)) {
			clusterRecord->mcid = 1;
		} /* else { // default 0 }*/
	}
}

void cluster_PETags2(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	int32_t i;
	worker_t w;
	int32_t mcid = 0; // cluster id
	
	w.opt = opt;
	w.clusterRecords = clusterRecords+1;
	clusterRecords[0].moid = 0;
	kt_for(opt->n_threads, PETag2_cluster_step1_worker, &w, nClusterKeys-1);
	
	// TODO: use serial approach to set up the correct cid
	for(i=0; i<nClusterKeys; ++i) {
		if (0==clusterRecords[i].mcid) {
#if 0
			mcid = i;
#else
			mcid = clusterRecords[i].moid;
#endif
			clusterRecords[i].mcid = mcid;
		} else if (1==clusterRecords[i].mcid) {
			clusterRecords[i].mcid = mcid;
		}
	}
}

void cluster_count_iPET(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	int32_t i;
	//worker_t w;
	
	//w.opt = opt;
	//w.clusterRecords = clusterRecords+1;
	//kt_for(opt->n_threads, PETag2_cluster_step1_worker, &w, nClusterKeys-1);
	
	// TODO: use serial approach to set up the correct cid
	// clusterRec->iPET = 1;
	if (nClusterKeys>1) {
		int32_t clusterIndex=0;
		int32_t cid=-1; // cluster id
		int32_t mcid=-1; // cluster id
		// TODO : optimization: we know where the usable end!!!
		for(i=0; i<nClusterKeys; ++i) {
			if (0==clusterRecords[i].usable) {
				break;
			}
			if (cid==clusterRecords[i].cid && mcid==clusterRecords[i].mcid) {
				clusterRecords[clusterIndex].iPET++;
				
				// get new two dimensional ids
				clusterRecords[i].cid = clusterRecords[clusterIndex].oid;
				clusterRecords[i].mcid = clusterRecords[clusterIndex].moid;
				
				// TODO: overloading
				clusterRecords[clusterIndex].swap = 1;
				
			} else {
				clusterRecords[i].iPET = 1;
				cid=clusterRecords[i].cid;
				mcid=clusterRecords[i].mcid;
				clusterIndex = i;

				// get new two dimensional ids
				clusterRecords[i].cid = clusterRecords[i].oid;
				clusterRecords[i].mcid = clusterRecords[i].moid;
				
				// TODO: overloading, if there is only a single record, nothing to check
				clusterRecords[i].swap = 0;
			}
		}
	}
}

void report_Total_iPET(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	int32_t i;
	int32_t totaliPET = 0;
	int32_t usableRows = 0;
	
	for(i=0; i<nClusterKeys; ++i) {
		if (0==clusterRecords[i].usable) {
			break;
		}
		totaliPET += clusterRecords[i].iPET;
		usableRows++;
	}
	
	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Total %d iPETs, %d usable iPETs, %d iPET colsum, %d self-ligated iPETs..\n", __func__,
				nClusterKeys, usableRows, totaliPET, nClusterKeys-usableRows);
	}
}

int32_t cluster_bin(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	int32_t i,j;
	int32_t nDiscordantPETs = 0;
	
	if (nClusterKeys>1) {
		// TODO : optimization: we know where the usable end!!!
		for(i=0; i<nClusterKeys; ++i) {
			if (0==clusterRecords[i].usable) {
				break;
			}
			
			// iPET will indicate the number of memebers in the bin
			//if (clusterRecords[i].iPET>1) {
			if (clusterRecords[i].swap) {
				// if >=2 members, we will need to check if the overlap is proper
				int32_t niPET = clusterRecords[i].iPET;
				
				// TODO: FOR DEBUGGING ONLY
				int nDebugDump = 0;
				int nBinDiscordantPETs = 0;
				
#if 0
				for(j=0; j<niPET; ++j) {
					if (
#if 1
						78979327
#endif
#if 0
						23110849
#endif
#if 0
						73199568
#endif
						==clusterRecords[i+j].pairId) {
						fprintf(stderr, "[WARNING:DEBUG] check record\n");
						nDebugDump = 1;
						break;
					}
				}
#endif
				
				// let's sort based on tag right and process them
				if (0!=nDebugDump) {
					fprintf(stdout, "[WARNING:DEBUG] initial \n");
					report_iPETs_List (&(clusterRecords[i]), niPET, NULL, -1);
				}
				qsort(&(clusterRecords[i]), niPET, sizeof(cpu_cluster_record_t), cmpPETTag2Order);
				// let's process sorted right tags
				//cluster_PETags2(opt, clusterKeys, nClusterKeys);
				clusterRecords[i].mcid = 0;
				for(j=1; j<niPET; ++j) {
					cpu_cluster_record_t *pClusterRecord = &(clusterRecords[i+j-1]);
					cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
					if (is_overlap(pClusterRecord->mpos, pClusterRecord->mendpos, clusterRecord->mpos, clusterRecord->mendpos)) {
						clusterRecord->mcid = 1;
					} else {
						clusterRecord->mcid = 0;
					}
				}
				int32_t mcid = 0; // cluster id
				for(j=0; j<niPET; ++j) {
					cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
					if (0==clusterRecord->mcid) {
						mcid = clusterRecord->moid;
						clusterRecord->mcid = mcid;
					} else if (1==clusterRecord->mcid) {
						clusterRecord->mcid = mcid;
					}
				}
				if (0!=nDebugDump) {
					fprintf(stdout, "[WARNING:DEBUG] Tag2 ordered \n");
					report_iPETs_List (&(clusterRecords[i]), niPET, NULL, -1);
				}
				
				// let's sort based on tag left and process them
				qsort(&(clusterRecords[i]), niPET, sizeof(cpu_cluster_record_t), cmpPETTag1Order);
				// let's process sorted left tags
				//cluster_PETags1(opt, clusterKeys, nClusterKeys);
				clusterRecords[i].cid = 0;
				for(j=1; j<niPET; ++j) {
					cpu_cluster_record_t *pClusterRecord = &(clusterRecords[i+j-1]);
					cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
					if (is_overlap(pClusterRecord->pos, pClusterRecord->endpos, clusterRecord->pos, clusterRecord->endpos)) {
						clusterRecord->cid = 1;
					} else {
						clusterRecord->cid = 0;
					}
				}
				int32_t cid = 0; // cluster id
				for(j=0; j<niPET; ++j) {
					cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
					if (0==clusterRecord->cid) {
						cid = clusterRecord->oid;
						clusterRecord->cid = cid;
					} else if (1==clusterRecord->cid) {
						clusterRecord->cid = cid;
					}
				}
				if (0!=nDebugDump) {
					fprintf(stdout, "[WARNING:DEBUG] Tag1 ordered \n");
					report_iPETs_List (&(clusterRecords[i]), niPET, NULL, -1);
				}
				
				// let's sort based on combine tag
				qsort(&(clusterRecords[i]), niPET, sizeof(cpu_cluster_record_t), cmpPETBothTagsOrder);
				// let's process them
				// TODO: cluster_count_iPET(opt, clusterKeys, nClusterKeys);
				int32_t clusterIndex=i;
				cid=-1; // cluster id
				mcid=-1; // cluster id
				for(j=0; j<niPET; ++j) {
					cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
					if (cid==clusterRecord->cid && mcid==clusterRecord->mcid) {
						clusterRecords[clusterIndex].iPET++;
						clusterRecords[clusterIndex].swap = 1;
						clusterRecord->swap = 0;
						
						clusterRecord->cid = clusterRecords[clusterIndex].oid;
						clusterRecord->mcid = clusterRecords[clusterIndex].moid;
						
						// folded reset into current loop
						clusterRecord->iPET = 0;
					} else {
						clusterRecord->iPET = 1;
						cid=clusterRecord->cid;
						mcid=clusterRecord->mcid;
						clusterIndex = i+j;
						
						clusterRecord->swap = 0;
						
						// update the 2D ids
						clusterRecord->cid = clusterRecord->oid;
						clusterRecord->mcid = clusterRecord->moid;
					}
				}
				if (0!=nDebugDump) {
					fprintf(stdout, "[WARNING:DEBUG] Combined Tag1 and Tag2 ordered \n");
					report_iPETs_List (&(clusterRecords[i]), niPET, NULL, -1);
				}
				
				{
					// for TESTING: check that we do not need further processing
					clusterIndex=i;
					cid=clusterRecords[i].cid; // cluster id
					mcid=clusterRecords[i].mcid; // cluster id
					int32_t s1=clusterRecords[i].pos;
					int32_t e1=clusterRecords[i].endpos;
					int32_t s2=clusterRecords[i].mpos;
					int32_t e2=clusterRecords[i].mendpos;
					int32_t discordantIndex = -1; // nope to start off with
					for(j=0; j<niPET; ++j) {
						cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
						if (cid==clusterRecord->cid && mcid==clusterRecord->mcid) {
							// same bin grouping, so must be overlapping
							clusterRecord->swap = 0;
							if (0==is_overlap(s1, e1, clusterRecord->pos, clusterRecord->endpos)) {
								// does not overlap! this is a problem, report it
								FILE *logFH = (0!=nDebugDump) ? stdout : stderr;
								fprintf(logFH, "[WARNING:tag1] i=%d, j=%d, iPET=%d, leader=%d, s1=%d(%d), e1=%d(%d), s2=%d, e2=%d, oid=%d(%d), moid=%d(%d)\n",
										i, j, niPET, clusterIndex, s1, clusterRecord->pos, e1, clusterRecord->endpos, s2, e2,
										clusterRecords[clusterIndex].oid, clusterRecord->oid,
										clusterRecords[clusterIndex].moid, clusterRecord->moid);
								nBinDiscordantPETs++;
								
								clusterRecords[clusterIndex].iPET--;
								if (-1==discordantIndex) {
									discordantIndex = i+j;
									clusterRecords[discordantIndex].cid = clusterRecords[discordantIndex].oid;
									clusterRecords[discordantIndex].mcid = clusterRecords[discordantIndex].moid;
									clusterRecords[discordantIndex].iPET = 1;
									clusterRecords[discordantIndex].swap = 0;
								} else {
									clusterRecord->cid = clusterRecords[discordantIndex].oid;
									clusterRecord->mcid = clusterRecords[discordantIndex].moid;
									clusterRecords[discordantIndex].iPET++;
									clusterRecords[discordantIndex].swap = 1;
								}
								
							} else {
								if (0==is_overlap(s2, e2, clusterRecord->mpos, clusterRecord->mendpos)) {
									// does not overlap! this is a problem, report it
									FILE *logFH = (0!=nDebugDump) ? stdout : stderr;
									fprintf(logFH, "[WARNING:tag2] i=%d, j=%d, iPET=%d, leader=%d, s1=%d, e1=%d, s2=%d(%d), e2=%d(%d), oid=%d(%d), moid=%d(%d)\n",
											i, j, niPET, clusterIndex, s1, e1, s2, clusterRecord->mpos, e2, clusterRecord->mendpos,
											clusterRecords[clusterIndex].oid, clusterRecord->oid,
											clusterRecords[clusterIndex].moid, clusterRecord->moid);
									nBinDiscordantPETs++;
									
									clusterRecords[clusterIndex].iPET--;
									if (-1==discordantIndex) {
										discordantIndex = i+j;
										clusterRecords[discordantIndex].cid = clusterRecords[discordantIndex].oid;
										clusterRecords[discordantIndex].mcid = clusterRecords[discordantIndex].moid;
										clusterRecords[discordantIndex].iPET = 1;
										clusterRecords[discordantIndex].swap = 0;
									} else {
										clusterRecord->cid = clusterRecords[discordantIndex].oid;
										clusterRecord->mcid = clusterRecords[discordantIndex].moid;
										clusterRecords[discordantIndex].iPET++;
										clusterRecords[discordantIndex].swap = 1;
									}
									
								} else {
									// both anchors overlapped, let's expand our anchor regions
									if (clusterRecord->pos<s1) {
										s1 = clusterRecord->pos;
									}
									if (e1<clusterRecord->endpos) {
										e1 = clusterRecord->endpos;
									}
									if (clusterRecord->mpos<s2) {
										s2 = clusterRecord->mpos;
									}
									if (e2<clusterRecord->mendpos) {
										e2 = clusterRecord->mendpos;
									}
								}
							}
							
						} else {
							// we are starting a new bin grouping, no checking
							cid=clusterRecord->cid;
							mcid=clusterRecord->mcid;
							clusterIndex = i+j;
							
							s1=clusterRecords[clusterIndex].pos;
							e1=clusterRecords[clusterIndex].endpos;
							s2=clusterRecords[clusterIndex].mpos;
							e2=clusterRecords[clusterIndex].mendpos;
							
							// at the end of this check, the group is confirmed
							// outliers will be expelled to discordantIndex set
							clusterRecords[clusterIndex].swap = 0;
							
							discordantIndex = -1;
						}
					}
					
					if (0!=nDebugDump) {
						fprintf(stdout, "[WARNING:DEBUG] After checking correction\n");
						report_iPETs_List (&(clusterRecords[i]), niPET, NULL, -1);
					}
					
					// we need to re-order the array
					if (nBinDiscordantPETs>0) {
						qsort(&(clusterRecords[i]), niPET, sizeof(cpu_cluster_record_t), cmpPETBothTagsOrder);
						
						if (0!=nDebugDump) {
							fprintf(stdout, "[WARNING:DEBUG] Final re-ordering\n");
							report_iPETs_List (&(clusterRecords[i]), niPET, NULL, -1);
						}
					}
					
					// TODO: check if recursion is a good way here as we might be unable multi-thread as a result
					if (nBinDiscordantPETs>0) {
#if 0
						fprintf(stdout, "[WARNING:DEBUG] Calling nested subclustering with discordantIndex=%d, %d discordant iPETs\n", discordantIndex, nBinDiscordantPETs);
						report_iPETs_List (&(clusterRecords[i]), niPET, NULL, -1);
#endif
						int32_t nRetryDiscordantPETs = cluster_bin(opt, &(clusterRecords[i]), niPET);
						if (0==nRetryDiscordantPETs) {
							// there is no discordant set and thus we attempt to merge if we can
						}
#if 0
						fprintf(stdout, "[WARNING:DEBUG] After nested subclustering with with discordantIndex=%d, %d discordant iPETs, %d discordard after trying\n", discordantIndex, nBinDiscordantPETs, nRetryDiscordantPETs);
						report_iPETs_List (&(clusterRecords[i]), niPET, NULL, -1);
#endif
					}
				}
				
				{
					// TODO: for DEBUGGING; check that we have not lost the iPET count in each of the original bin group
					int32_t niPETinBin = 0;
					for(j=0; j<niPET; ++j) {
						cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
						niPETinBin += clusterRecord->iPET;
					}
					if (niPET!=niPETinBin) {
						fprintf(stderr, "[WARNING:iPET] i=%d, iPET=%d, iPETinBin=%d\n", i, niPET, niPETinBin);
						fprintf(stderr, "[WARNING:iPET] i=%d, iPETs=", i);
						for(j=0; j<niPET; ++j) {
							cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
							fprintf(stderr, "%d,", clusterRecord->iPET);
						}
						fprintf(stderr, "\n");
					}
				}
				
				
				i += (niPET-1);
				
				nDiscordantPETs += nBinDiscordantPETs;
			}
		}
		
		
		// TODO: should move this elsewhere or as a routine to call
		//       update the leader of the bin to encompass the anchor region
		if (0==nDiscordantPETs) {
			for(i=0; i<nClusterKeys; ++i) {
				if (0==clusterRecords[i].usable) {
					break;
				}
				
				if (clusterRecords[i].iPET>1) {
					cpu_cluster_record_t *clusterRecord = &(clusterRecords[i]);
					for(j=1; j<clusterRecords[i].iPET; ++j) {
						clusterRecord++;
						
						if (clusterRecord->pos<clusterRecords[i].pos) {
							clusterRecords[i].pos = clusterRecord->pos;
						}
						if (clusterRecords[i].endpos<clusterRecord->endpos) {
							clusterRecords[i].endpos = clusterRecord->endpos;
						}
						if (clusterRecord->mpos<clusterRecords[i].mpos) {
							clusterRecords[i].mpos = clusterRecord->mpos;
						}
						if (clusterRecords[i].mendpos<clusterRecord->mendpos) {
							clusterRecords[i].mendpos = clusterRecord->mendpos;
						}
					}
				}
			}
		}
	}
	
	return nDiscordantPETs;
}

int cmpSubclusterTag1(const void *a, const void *b) {
	const cpu_subcluster_record_t *ia = (const cpu_subcluster_record_t *)a;
	const cpu_subcluster_record_t *ib = (const cpu_subcluster_record_t *)b;
	
	// sort by start1, end1,
	int i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->lineno - ia->lineno; //if (!i32comp) return i32comp;
	
	return i32comp;
}

int cmpSubclusterTag2(const void *a, const void *b) {
	const cpu_subcluster_record_t *ia = (const cpu_subcluster_record_t *)a;
	const cpu_subcluster_record_t *ib = (const cpu_subcluster_record_t *)b;
	
	int i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->lineno - ia->lineno; //if (!i32comp) return i32comp;

	return i32comp;
}

int cmpSubclusterBothTags(const void *a, const void *b) {
	const cpu_subcluster_record_t *ia = (const cpu_subcluster_record_t *)a;
	const cpu_subcluster_record_t *ib = (const cpu_subcluster_record_t *)b;
	
	int32_t i32comp;
	i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	
	// TODO: any proxy for optimization?
	
	// sort by start1, end1,
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	// follow by start2, end2
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->lineno - ia->lineno; //if (!i32comp) return i32comp;
	
	return i32comp;
}

void report_Subclusters_List (FILE *clustersfd, const cpu_subcluster_record_t *subclusters, uint32_t nSubclusters, const bam_header_t *header) {
	int32_t i;
	
	fprintf(stdout, "#pairId\toffset\tcid\tmcid\toid\tmoid\tiPETs\tpos\tendpos\tmpos\tmendpos\n");
	
	for(i=0; i<nSubclusters; ++i) {
		
		const cpu_subcluster_record_t *subcluster = &(subclusters[i]);
#if 1
		fprintf(stderr, "%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				subcluster->pairId, subcluster->offset,
				subcluster->cid, subcluster->mcid,
				subcluster->oid, subcluster->moid,
				subcluster->iPET,
				subcluster->pos, subcluster->endpos,
				subcluster->mpos, subcluster->mendpos);
#else
		kstring_t str;
		str.l = str.m = 0; str.s = 0;
		
		kputul(subcluster->pairId, &str);
		kputc('\t', &str); kputl(subcluster->offset, &str);
		
		kputc('\t', &str); kputl(subcluster->cid, &str);
		kputc('\t', &str); kputl(subcluster->mcid, &str);
		kputc('\t', &str); kputl(subcluster->oid, &str);
		kputc('\t', &str); kputl(subcluster->moid, &str);
		
		kputc('\t', &str); kputl(subcluster->iPET, &str);
		
		kputc('\t', &str); kputw(subcluster->pos, &str);
		kputc('\t', &str); kputw(subcluster->endpos, &str);
		
		kputc('\t', &str); kputw(subcluster->mpos, &str);
		kputc('\t', &str); kputw(subcluster->mendpos, &str);
		
		kputc('\n', &str);
		
		fprintf(stdout, "%s", str.s);
		free(str.s);
#endif
	}
}

int32_t merge_subclusters (const cluster_opt_t *opt, cpu_subcluster_record_t *subclusters, uint32_t nSubclusters, uint32_t nFinalSubclusters) {
	int32_t numSubclusters = nSubclusters;
	int32_t j;
	
	FILE *logFH	= stderr;
	//FILE *logFH	= stdout;
	
	if (nSubclusters>1) {
		// TODO: FOR DEBUGGING ONLY
		//int nDebugDump = 0;
		int nDebugDump = 1;
		
#if 0
		for(j=0; j<niPET; ++j) {
			if (
#if 1
				78979327
#endif
				==subclusters[j].pairId) {
				fprintf(logFH, "[WARNING:DEBUG:SubCluster] check record\n");
				nDebugDump = 1;
				break;
			}
		}
#endif
		
		// let's sort based on tag right and process them
		if (0!=nDebugDump) {
			fprintf(logFH, "[WARNING:DEBUG:SubCluster] initial \n");
			report_Subclusters_List (logFH, subclusters, nSubclusters, NULL);
		}
		
		qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterTag2);
		// let's process sorted right tags
		//cluster_PETags2(opt, clusterKeys, nClusterKeys);
		subclusters[0].mcid = 0;
		for(j=1; j<nSubclusters; ++j) {
			cpu_subcluster_record_t *pSubcluster = &(subclusters[j-1]);
			cpu_subcluster_record_t *subcluster = &(subclusters[j]);
			if (is_overlap(pSubcluster->mpos, pSubcluster->mendpos, subcluster->mpos, subcluster->mendpos)) {
				subcluster->mcid = 1;
			} else {
				subcluster->mcid = 0;
			}
		}
		int32_t mcid = 0; // cluster id
		for(j=0; j<nSubclusters; ++j) {
			cpu_subcluster_record_t *subcluster = &(subclusters[j]);
			if (0==subcluster->mcid) {
				mcid = subcluster->moid;
				subcluster->mcid = mcid;
			} else if (1==subcluster->mcid) {
				subcluster->mcid = mcid;
			}
		}
		
		if (0!=nDebugDump) {
			fprintf(logFH, "[WARNING:DEBUG:SubCluster] Tag2 ordered \n");
			report_Subclusters_List (logFH, subclusters, nSubclusters, NULL);
		}
		
		// let's sort based on tag left and process them
		qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterTag1);
		// let's process sorted left tags
		//cluster_PETags1(opt, clusterKeys, nClusterKeys);
		subclusters[0].cid = 0;
		for(j=1; j<nSubclusters; ++j) {
			cpu_subcluster_record_t *pSubcluster = &(subclusters[j-1]);
			cpu_subcluster_record_t *subcluster = &(subclusters[j]);
			if (is_overlap(pSubcluster->pos, pSubcluster->endpos, subcluster->pos, subcluster->endpos)) {
				subcluster->cid = 1;
			} else {
				subcluster->cid = 0;
			}
		}
		int32_t cid = 0; // cluster id
		for(j=0; j<nSubclusters; ++j) {
			cpu_subcluster_record_t *subcluster = &(subclusters[j]);
			if (0==subcluster->cid) {
				cid = subcluster->oid;
				subcluster->cid = cid;
			} else if (1==subcluster->cid) {
				subcluster->cid = cid;
			}
		}
		
		if (0!=nDebugDump) {
			fprintf(logFH, "[WARNING:DEBUG:SubCluster] Tag1 ordered \n");
			report_Subclusters_List (logFH, subclusters, nSubclusters, NULL);
		}
		
		// let's sort based on combine tag
		qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterBothTags);
		// let's process them
		int32_t subclusterIndex=0;
		cid=-1; // cluster id
		mcid=-1; // cluster id
		numSubclusters = 0;
		for(j=0; j<nSubclusters; ++j) {
			cpu_subcluster_record_t *subcluster = &(subclusters[j]);
			if (cid==subcluster->cid && mcid==subcluster->mcid) {
				//subclusters[subclusterIndex].iPET+=subcluster->iPET;
				
				subcluster->cid = subclusters[subclusterIndex].oid;
				subcluster->mcid = subclusters[subclusterIndex].moid;
				
				// NOTE: if we fold here but find that it is not overlapping, we have to unfold! Delay folding!
				// folded reset into current loop
				// subcluster->iPET = 0;
			} else {
				//subcluster->iPET = 1;
				cid=subcluster->cid;
				mcid=subcluster->mcid;
				subclusterIndex = j;
				
				// update the 2D ids
				subcluster->cid = subcluster->oid;
				subcluster->mcid = subcluster->moid;
				
				numSubclusters++;
			}
		}
		
		if (0!=nDebugDump) {
			fprintf(logFH, "[WARNING:DEBUG:SubCluster] Combined Tag1 and Tag2 ordered, estimated %d subclusters\n", numSubclusters);
			report_Subclusters_List (logFH, subclusters, nSubclusters, NULL);
		}
		
		// TODO: check condition of merging and update the ranges
		{
			int32_t nMergedSubclusters = 0;
			// for TESTING: check that we do not need further processing
			int32_t subclusterIndex=0;
			cid=-1; // cluster id
			mcid=-1; // cluster id
			int32_t s1=-1;
			int32_t e1=-1;
			int32_t s2=-1;
			int32_t e2=-1;
			int32_t discordantIndex = -1; // nope to start off with
			for(j=0; j<nSubclusters; ++j) {
				cpu_subcluster_record_t *subcluster = &(subclusters[j]);
				if (cid==subcluster->cid && mcid==subcluster->mcid) {
					// same bin grouping, so must be overlapping
					if (0==is_overlap(s1, e1, subcluster->pos, subcluster->endpos)) {
						// does not overlap! this is a problem, report it
						FILE *logFH = (0!=nDebugDump) ? stdout : stderr;
						fprintf(logFH, "[WARNING:SubCluster:tag1] j=%d, iPET=%d, leader=%d, s1=%d(%d), e1=%d(%d), s2=%d, e2=%d, oid=%d(%d), moid=%d(%d)\n",
								j, subcluster->iPET, subclusterIndex, s1, subcluster->pos, e1, subcluster->endpos, s2, e2,
								subclusters[subclusterIndex].oid, subcluster->oid,
								subclusters[subclusterIndex].moid, subcluster->moid);
						//nMergedSubclusters++;
						
						//subcluster[subclusterIndex].iPET--;
						if (-1==discordantIndex) {
							discordantIndex = j;
							subclusters[discordantIndex].cid = subclusters[discordantIndex].oid;
							subclusters[discordantIndex].mcid = subclusters[discordantIndex].moid;
							//subclusters[discordantIndex].iPET = 1;
							//subclusters[discordantIndex].swap = 0;
						} else {
							subcluster->cid = subclusters[discordantIndex].oid;
							subcluster->mcid = subclusters[discordantIndex].moid;
							//subclusters[discordantIndex].iPET++;
							//subclusters[discordantIndex].swap = 1;
						}
						
					} else {
						if (0==is_overlap(s2, e2, subcluster->mpos, subcluster->mendpos)) {
							// does not overlap! this is a problem, report it
							//FILE *logFH = (0!=nDebugDump) ? stdout : stderr;
							fprintf(logFH, "[WARNING:SubCluster:tag2] j=%d, iPET=%d, leader=%d, s1=%d, e1=%d, s2=%d(%d), e2=%d(%d), oid=%d(%d), moid=%d(%d)\n",
									j, subcluster->iPET, subclusterIndex, s1, e1, s2, subcluster->mpos, e2, subcluster->mendpos,
									subclusters[subclusterIndex].oid, subcluster->oid,
									subclusters[subclusterIndex].moid, subcluster->moid);
							//nMergedSubclusters++;
							
							//subclusters[clusterIndex].iPET--;
							if (-1==discordantIndex) {
								discordantIndex = j;
								subclusters[discordantIndex].cid = subclusters[discordantIndex].oid;
								subclusters[discordantIndex].mcid = subclusters[discordantIndex].moid;
								//subclusters[discordantIndex].iPET = 1;
								//subclusters[discordantIndex].swap = 0;
							} else {
								subcluster->cid = subclusters[discordantIndex].oid;
								subcluster->mcid = subclusters[discordantIndex].moid;
								//subclusters[discordantIndex].iPET++;
								//subclusters[discordantIndex].swap = 1;
							}
							
						} else {
							// both anchors overlapped, let's expand our anchor regions
							if (subcluster->pos<s1) {
								s1 = subcluster->pos;
							}
							if (e1<subcluster->endpos) {
								e1 = subcluster->endpos;
							}
							if (subcluster->mpos<s2) {
								s2 = subcluster->mpos;
							}
							if (e2<subcluster->mendpos) {
								e2 = subcluster->mendpos;
							}
							nMergedSubclusters++;
						}
					}
					
				} else {
					// we are starting a new bin grouping, no checking
					cid=subcluster->cid;
					mcid=subcluster->mcid;
					subclusterIndex = j;
					
					s1=subclusters[subclusterIndex].pos;
					e1=subclusters[subclusterIndex].endpos;
					s2=subclusters[subclusterIndex].mpos;
					e2=subclusters[subclusterIndex].mendpos;
					
					// at the end of this check, the group is confirmed
					// outliers will be expelled to discordantIndex set
					
					discordantIndex = -1;
				}
			}
			
			if (0!=nDebugDump) {
				fprintf(logFH, "[WARNING:DEBUG:SubCluster] After checking correction, with %d merger(s)\n", nMergedSubclusters);
				report_Subclusters_List (logFH, subclusters, nSubclusters, NULL);
			}
			
			// we need to re-order the array
			if (nMergedSubclusters>0) {
				qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterBothTags);
				
				if (0!=nDebugDump) {
					fprintf(logFH, "[WARNING:DEBUG:SubCluster] Final re-ordering\n");
					report_Subclusters_List (logFH, subclusters, nSubclusters, NULL);
				}
				
				// TODO: check condition of invocation & termination
				//if (numSubclusters!=nSubclusters) {
				//	int32_t final_numSubclusters = merge_subclusters (opt, &(subclusters[i]), numSubclusters, numSubclusters);
				//}
			}
		}
	}
	
	return numSubclusters;
}

int32_t merge_subclusters_in_bin(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t niPET) {
	int32_t numSubclusters = 0;
	int32_t i,j;
	{
		// report the number of unique subcluster within a bin (>1 members)
		for(j=0; j<niPET; ++j) {
			cpu_cluster_record_t *clusterRecord = &(clusterRecords[j]);
			if (clusterRecord->iPET>0) {
				numSubclusters++;
			}
		}
		
		if (numSubclusters>1) {
			kstring_t str;
			str.l = str.m = 0; str.s = 0;
			kputs("[WARNING:SubCluster] pairId=", &str); kputul(clusterRecords[0].pairId, &str);
			kputs(" iPET=", &str); kputl(niPET, &str);
			kputs(" #subcluster=", &str); kputl(numSubclusters, &str);
			kputc('\n', &str);
			fprintf(stderr, "%s", str.s);
			free(str.s);
		}
	}

	if (numSubclusters>1) {
		// prepare the summarized range
		cpu_subcluster_record_t *subclusters = (cpu_subcluster_record_t *) calloc(numSubclusters, sizeof(cpu_subcluster_record_t));
		i = 0;
		for(j=0; j<niPET; ++j) {
			cpu_cluster_record_t *clusterRecord = &(clusterRecords[j]);
			if (clusterRecord->iPET>0) {
				cpu_subcluster_record_t *subcluster = &(subclusters[i]);
				subcluster->pos = clusterRecord->pos; subcluster->endpos = clusterRecord->endpos;
				subcluster->mpos = clusterRecord->mpos; subcluster->mendpos = clusterRecord->mendpos;
				subcluster->cid = clusterRecord->cid; subcluster->mcid = clusterRecord->mcid;
				subcluster->iPET = clusterRecord->iPET;
				subcluster->oid = clusterRecord->oid; subcluster->moid = clusterRecord->moid;
				subcluster->offset = i;
				subcluster->lineno = clusterRecord->lineno;
				subcluster->pairId = clusterRecord->pairId;
				i++;
			}
		}
		
		// attempt to merge the ranges
		int32_t final_numSubclusters = merge_subclusters (opt, subclusters, numSubclusters, numSubclusters);
		
		// TODO: propagate the new merged ranges back to the cluster record!
		
		free(subclusters);
	}
	
	return numSubclusters;
}

int32_t cluster_iPETs(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	int32_t i,j;
	if (nClusterKeys>1) {
		// TODO : optimization: we know where the usable end!!!
		for(i=0; i<nClusterKeys; ++i) {
			if (0==clusterRecords[i].usable) {
				break;
			}
			
			// iPET will indicate the number of memebers in the bin
			//if (clusterRecords[i].iPET>1) {
			if (clusterRecords[i].swap) {
				// if >=2 members, we will need to check if the overlap is proper
				int32_t niPET = clusterRecords[i].iPET;
				int32_t nRetryDiscordantPETs = cluster_bin(opt, &(clusterRecords[i]), niPET);
				int32_t numSubclusters = merge_subclusters_in_bin(opt, &(clusterRecords[i]), niPET);
				
				{
					// TODO: for DEBUGGING; check that we have not lost the iPET count in each of the original bin group
					int32_t niPETinBin = 0;
					for(j=0; j<niPET; ++j) {
						cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
						niPETinBin += clusterRecord->iPET;
					}
					if (niPET!=niPETinBin) {
						fprintf(stderr, "[WARNING:iPET] i=%d, iPET=%d, iPETinBin=%d\n", i, niPET, niPETinBin);
						fprintf(stderr, "[WARNING:iPET] i=%d, iPETs=", i);
						for(j=0; j<niPET; ++j) {
							cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
							fprintf(stderr, "%d,", clusterRecord->iPET);
						}
						fprintf(stderr, "\n");
					}
				}
				
				// TODO: we will now try to merge the subcluster(s) in this bin if possible
				
				i += (niPET-1);
			}
		}
	}
	
	return 0;
}

void report_Clusters_List (FILE *clustersfd, const cpu_cluster_record_t *clusterKeys, uint32_t nClusterKeys, const bam_header_t *header) {
	int32_t i;
	int32_t cid=0; // cluster id
	int32_t mcid=0; // cluster id
	for(i=0; i<nClusterKeys; ++i) {
		const cpu_cluster_record_t *clusterRec = &(clusterKeys[i]);
		int32_t toPrint = 0;
		
		if (0==clusterRec->usable) {
			// just dump the record
			toPrint = 0;
		} else {
			// make sure that it is different before printing
			if (cid!=clusterRec->cid || mcid!=clusterRec->mcid) {
				toPrint = 1;
				
				cid = clusterRec->cid;
				mcid = clusterRec->mcid;
			}
		}
		
		if (0!=toPrint) {
			kstring_t str;
			str.l = str.m = 0; str.s = 0;
			
			if (clusterRec->tid < 0) kputsn("*\t0\t0", 5, &str);
			else {
				if (header) kputs(header->target_name[clusterRec->tid] , &str);
				else kputw(clusterRec->tid, &str);
				kputc('\t', &str);kputw(clusterRec->pos, &str);
				kputc('\t', &str);kputw(clusterRec->endpos, &str);
			}
			kputsn("\t0\t0\t0", 6, &str);
			
			if (clusterRec->mtid < 0) kputsn("\t*\t0\t0", 6, &str);
			else {
				kputc('\t', &str);
				if (header) kputs(header->target_name[clusterRec->mtid] , &str);
				else kputw(clusterRec->mtid, &str);
				kputc('\t', &str); kputw(clusterRec->mpos, &str);
				kputc('\t', &str); kputw(clusterRec->mendpos, &str);
			}
			kputsn("\t0\t0\t0", 6, &str);
			
			kputc('\t', &str); kputl(clusterRec->iPET, &str);
			kputc('\t', &str); kputw((clusterRec->tid == clusterRec->mtid) ? 1 : 0, &str);
			
			int32_t dist;
			if (clusterRec->tid != clusterRec->mtid) {
				dist = 0x7FFFFFFF;
			} else if (-1==clusterRec->tid) {
				dist = 0x7FFFFFFF;
			} else {
				dist = ((clusterRec->mpos + clusterRec->mendpos) - (clusterRec->pos + clusterRec->endpos)) / 2 ;
			}
			kputc('\t', &str); kputl(dist, &str);
			
			kputsn("\t0.0\t---", 8, &str);
			kputc('\n', &str);
			
			fprintf(clustersfd, "%s", str.s);
			free(str.s);
		}
	}
}

int main_cluster(int argc, char *argv[])
{
	int c, n, i;
	char *p = 0;
	cluster_opt_t *opt;
	
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
	
	FILE *outfd = NULL;
	kstring_t clusterKeyFilename = {0,0,0};

	FILE *clustersfd = NULL;
	kstring_t clustersFilename = {0,0,0};
	
	init_CPUBamBuffer(&bamsBuffer);
	
	opt = cluster_opt_init();
	strcpy(in_mode, "r");
	while ((c = getopt(argc, argv, "disSt:O:l:u:b:5:3:m:c:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'l') opt->outputSortedList = atoi(optarg), opt->outputSortedList = opt->outputSortedList > 0 ? opt->outputSortedList : 0;
		else if (c == 'c') opt->cycle = atoi(optarg);
		else if (c == 'd') opt->disponoParallel = 1;
		else if (c == 'S') is_bamin = 0;
		else if (c == 'O') { ks_set(optarg, &(opt->outputPrefix)); }
		else if (c == 's') opt->selfLigation = atoi(optarg), opt->selfLigation = opt->selfLigation > 0 ? opt->selfLigation : G_SELF_LIGATION;
		else if (c == '5') {
			opt->extension5p.source = (int32_t) strtol(optarg, &p, 10);
			// TODO: check that it is 5' or 3'
			if (*p != 0 && ispunct(*p) && (isdigit(p[1])||'+'==p[1]||'-'==p[1]))
				opt->extension5p.size = (int32_t) strtol(p+1, &p, 10);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] 5' extension by %d bp from %d'\n",
						__func__, opt->extension5p.size, opt->extension5p.source);
		}
		else if (c == '3') {
			opt->extension3p.source = (int32_t) strtol(optarg, &p, 10);
			// TODO: check that it is 5' or 3'
			if (*p != 0 && ispunct(*p) && (isdigit(p[1])||'+'==p[1]||'-'==p[1]))
				opt->extension3p.size = (int32_t) strtol(p+1, &p, 10);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] 3' extension by %d bp from %d'\n",
						__func__, opt->extension3p.size, opt->extension3p.source);
		} else {
			cluster_opt_terminate(opt);
			free(opt);
			return 1;
		}
	}
	
	if (opt->n_threads < 1) opt->n_threads = 1;
	
	if (is_bamin) strcat(in_mode, "b");

	if (optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cpu cluster [options] <in.sam/.bam>\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -d         use parallel sort\n");
		fprintf(stderr, "       -5 INT,INT extension source %d' and size [%d bp]\n", opt->extension5p.source, opt->extension5p.size);
		fprintf(stderr, "       -3 INT,INT extension source %d' and size [%d bp]\n", opt->extension3p.source, opt->extension3p.size);
		fprintf(stderr, "       -s INT     self-ligation distance [<=%d bp]\n", opt->selfLigation);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -S         input sam format\n");
		fprintf(stderr, "       -l         output sorted PETs to stdout\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		cluster_opt_terminate(opt);
		free(opt);
		return 1;
	}
	
	// open file handlers
	if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open \"%s\" for reading.\n", __func__, argv[optind]);
		
		free(fn_list);
		cluster_opt_terminate(opt);
		free(opt);
		
		return 1;
	}
	if (in->header == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to read the header from \"%s\".\n", __func__, argv[optind]);
		
		free(fn_list);
		// close files, free and return
		samclose(in);
		
		cluster_opt_terminate(opt);
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
	
	// we keep the data that we are going to sort in a binary files for single read processing
	// in future, if we have decided on merge sort, we can perform other form of combining the results
	// also, we can already sort the entries before writing them into the binary files
	{
		ksprintf(&clusterKeyFilename, "%s.cpu.cluster", opt->outputPrefix.s);
		const char *mode = "wb";
		outfd = fopen(clusterKeyFilename.s, mode);
		if (outfd == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing: %s\n", __func__, clusterKeyFilename.s, strerror (errno));
			cluster_opt_terminate(opt);
			free(opt);
			return 1;
		}
		
		ksprintf(&clustersFilename, "%s.clusters.txt", opt->outputPrefix.s);
		clustersfd = fopen(clustersFilename.s, "w");
		if (clustersfd == NULL) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing: %s\n", __func__, clustersFilename.s, strerror (errno));
			cluster_opt_terminate(opt);
			free(opt);
			fclose(outfd);
			return 1;
		}
	}
	
	// we write header information for future usage
	// TODO: can we safely read text lines from a binary files?
	//       if so, we will like to record the software version and command line
	// TODO: get a unsigned int32 version of the record
	uint32_t nValue = (0x00<<16 | 0x01 <<8 | 0x00); err_fwrite(&nValue, sizeof(nValue), 1, outfd); // cpu version
	nValue = 0; err_fwrite(&nValue, sizeof(nValue), 1, outfd); // sort state: unsorted=0, chunk sorted=1, fully sorted=2
	nValue = sizeof(cpu_cluster_record_t); err_fwrite(&nValue, sizeof(nValue), 1, outfd); // record size
	err_fwrite(&n_pairProcessed, sizeof(n_pairProcessed), 1, outfd);// filler for total number of records
	
	// PROCESSING
	if (argc == optind + 1) { // convert/print the entire file
		if (CPU_DEDUP_DEBUG_CONVERT_LIST==(CPU_DEDUP_DEBUG_CONVERT_LIST&opt->outputSortedList)) {
			fprintf(stdout, "\n\n===\tCONVERTED LIST\n");
		}
		
		int nNumRequested = opt->chunk_size * opt->n_threads;
		t_timepointProcessing = realtime();
		while ((bamRecs = readCPUBam(&n, nNumRequested, in, &bamsBuffer)) != 0) {
			t_timepointIO = realtime();
			t_diffIO = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] read %d records..\n", __func__, n-bamsBuffer.unprocessed);
			
			// processing the pairing
			cpu_readid_t *readids = calloc(n, sizeof(cpu_readid_t));
			tag_group_t *tagGroups = calloc(n, sizeof(tag_group_t));
			cpu_cluster_record_t *clusterRecords = calloc(n, sizeof(cpu_cluster_record_t));
			
			uint32_t numTagGroups = 0;
			cluster_process_alns(opt, n_processed, n, bamRecs, readids, &numTagGroups, tagGroups, clusterRecords, &bamsBuffer, 0);
			t_timepointProcessing = realtime();
			t_diffProcessing = t_timepointProcessing - t_timepointIO;
			
			// TODO: write the bam results
			// TODO: for testing purposes
			// TODO: need to make sure that we have either 1 record (SE) or 2 records (PE)
			
			if (CPU_DEDUP_DEBUG_CONVERT_LIST==(CPU_DEDUP_DEBUG_CONVERT_LIST&opt->outputSortedList)) {
				report_iPETs_List (clusterRecords, numTagGroups, in->header, -1);
			}
			
			// TODO: write the binary data file
			err_fwrite(&numTagGroups, sizeof(numTagGroups), 1, outfd); // write the number of record filler for total number of records
			uint32_t nSorted = 0; err_fwrite(&nSorted, sizeof(nSorted), 1, outfd); // sort state: unsorted=0, chunk sorted=1, fully sorted=2
			err_fwrite(clusterRecords, sizeof(cpu_cluster_record_t), numTagGroups, outfd);
			
			t_timepointIO = realtime();
			t_diffIO += (t_timepointIO - t_timepointProcessing);
			// clean up
			for(i=0; i<n; ++i) free(bamRecs[i].data);
			free(bamRecs);
			free(readids);
			free(clusterRecords);
			
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
		updateClusterTotalRecords(clusterKeyFilename.s, n_pairProcessed);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] update record count, i/o %.2f sec (%.2f min)..\n", __func__, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// we next read all the binary records for sorting
		uint32_t nClusterKeys = 0;
		t_timepointProcessing = realtime();
		cpu_cluster_record_t *clusterKeys = readClusterKeys(clusterKeyFilename.s, &nClusterKeys);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] read %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		if (!clusterKeys) {
			fprintf(stderr, "[E::%s] No duplicated key reads, terminating..\n", __func__);
			// TODO: clean up
		}
		
		if (CPU_DEDUP_DEBUG_READ_LIST==(CPU_DEDUP_DEBUG_READ_LIST&opt->outputSortedList)) {
			// for debugging purposes only
			t_timepointProcessing = realtime();
			fprintf(stdout, "\n\n===\tREAD LIST\n");
			report_iPETs_List (clusterKeys, nClusterKeys, in->header, -1);
			fprintf(stdout, "===\tEND : READ LIST\n");
			fflush(stdout);
			t_timepointIO = realtime();
			t_diffProcessing = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3) {
				fprintf(stderr, "[M::%s] Writing read %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
			}
			// END - for debugging purposes only, remove from production code
		}
		
		// let's extend the tags accordingly
		t_timepointProcessing = realtime();
		extend_PETs(opt, clusterKeys, nClusterKeys, in->header->target_len, in->header->n_targets);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Extend %u PET records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		
		// let's sort based on tag right and process them
		t_timepointProcessing = realtime();
		qsort(clusterKeys, nClusterKeys, sizeof(cpu_cluster_record_t), cmpPETTag2);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Sorted %u key right records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// let's process sorted right tags
		t_timepointProcessing = realtime();
		cluster_PETags2(opt, clusterKeys, nClusterKeys);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Clustered %u key right records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// let's sort based on tag left and process them
		t_timepointProcessing = realtime();
		qsort(clusterKeys, nClusterKeys, sizeof(cpu_cluster_record_t), cmpPETTag1);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Sorted %u left key records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// let's process sorted left tags
		t_timepointProcessing = realtime();
		cluster_PETags1(opt, clusterKeys, nClusterKeys);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Clustered %u left key records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// let's sort based on combine tag
		t_timepointProcessing = realtime();
		qsort(clusterKeys, nClusterKeys, sizeof(cpu_cluster_record_t), cmpPETBothTags);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Sorted %u combined key records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// let's process them
		t_timepointProcessing = realtime();
		cluster_count_iPET(opt, clusterKeys, nClusterKeys);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Clustered %u key right records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		// for DEBUGGING
		report_Total_iPET(opt, clusterKeys, nClusterKeys);
		if (CPU_DEDUP_DEBUG_SORTED_LIST==(CPU_DEDUP_DEBUG_SORTED_LIST&opt->outputSortedList) && (0==opt->cycle)) {
			// for debugging purposes only, remove from production code
			t_timepointProcessing = realtime();
			fprintf(stdout, "\n\n===\tSORTED LIST - Check Cycle #%d\n", 0);
			report_iPETs_List (clusterKeys, nClusterKeys, in->header, -1);
			fprintf(stdout, "===\tEND : SORTED LIST - Check Cycle #%d\n", 0);
			fflush(stdout);
			
			t_timepointIO = realtime();
			t_diffProcessing = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3) {
				fprintf(stderr, "[M::%s] Writing sorted %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
			}
			// TODO: END - for debugging purposes only, remove from production code
		}
		
		// let's sub-cluster the approximated cluster bin
		int32_t nCheckDepth = 0;
		int32_t nDiscordantPETs = 0;
		do {
			t_timepointProcessing = realtime();
			nDiscordantPETs = cluster_iPETs(opt, clusterKeys, nClusterKeys);
			t_timepointIO = realtime();
			t_diffProcessing = t_timepointIO - t_timepointProcessing;
			nCheckDepth++;
			if (bwa_verbose >= 3) {
				fprintf(stderr, "[M::%s] Check#%d, Sub-clustered %u key right records, %d remaining discordant PETs, processing %.2f sec ( %.2f min)..\n", __func__, nCheckDepth, nClusterKeys, nDiscordantPETs, t_diffProcessing, t_diffProcessing/60.0);
			}
			
			// for DEBUGGING
			report_Total_iPET(opt, clusterKeys, nClusterKeys);
			
			if (CPU_DEDUP_DEBUG_SORTED_LIST==(CPU_DEDUP_DEBUG_SORTED_LIST&opt->outputSortedList) && (opt->cycle==nCheckDepth)) {
				// for debugging purposes only, remove from production code
				t_timepointProcessing = realtime();
				fprintf(stdout, "\n\n===\tSORTED LIST - Check Cycle #%d\n", nCheckDepth);
				report_iPETs_List (clusterKeys, nClusterKeys, in->header, -1);
				fprintf(stdout, "===\tEND : SORTED LIST - Check Cycle #%d\n", nCheckDepth);
				fflush(stdout);
				
				t_timepointIO = realtime();
				t_diffProcessing = t_timepointIO - t_timepointProcessing;
				if (bwa_verbose >= 3) {
					fprintf(stderr, "[M::%s] Writing sorted %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
				}
				// TODO: END - for debugging purposes only, remove from production code
			}
		} while (nDiscordantPETs>0);
		
		// TOOD:future
		// writeDedupKeys(clusterKeyFilename.s, dedupKeys, nDedupKeys);

		if (CPU_DEDUP_DEBUG_SORTED_LIST==(CPU_DEDUP_DEBUG_SORTED_LIST&opt->outputSortedList) && opt->cycle<0) {
			// for debugging purposes only, remove from production code
			t_timepointProcessing = realtime();
			fprintf(stdout, "\n\n===\tSORTED LIST\n");
			report_iPETs_List (clusterKeys, nClusterKeys, in->header, -1);
			fprintf(stdout, "===\tEND : SORTED LIST\n");
			fflush(stdout);
			
			t_timepointIO = realtime();
			t_diffProcessing = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3) {
				fprintf(stderr, "[M::%s] Writing sorted %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
			}
			// TODO: END - for debugging purposes only, remove from production code
		}
		
		// print results out in same format for comparison
		report_Clusters_List (clustersfd, clusterKeys, nClusterKeys, in->header);
		
		// release the sorted records as we no longer need them
		free(clusterKeys);
		
		//-------------------------------------------------------------------------------------------------------------------
		
		samclose(in);
		
		// we reset the counters
		n_processed = 0;
		n_pairProcessed = 0;
		// TODO: we loop to compute and output the cluster information
	}
	// END - PROCESSING
	
	free(fn_list);
	
	// close files, free and return
	// TODO: decide if we need additional closing due to re-opening
	//samclose(in);
	
	terminate_CPUBamBuffer(&bamsBuffer);
	
	//TODO: fclose(outfd);
	free(clusterKeyFilename.s);
	
	//TODO: fclose(clustersfd);
	free(clustersFilename.s);
	
	cluster_opt_terminate(opt);
	free(opt);
	
	return ret;
}


