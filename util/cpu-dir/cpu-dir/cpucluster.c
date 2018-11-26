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
// [X] 7. -j with gzip output
// [X] 8. -x will exclude self-ligation in -j output
// [X] 9. -g to write inter-chromosomal and intra-chromosomal separately
// [X] 10. -f potentially writing gzip if the downstream can take it

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

#define SUBCLUSTER_ITERACTION_BLOCK 200
#if 1
#define BULK_IPET_TRIGGER			1000
#define BULK_IPET_BLOCK				 500
#else
#define BULK_IPET_TRIGGER			5000
#define BULK_IPET_BLOCK				2000
#endif

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

#define CPU_CLUSTER_JUICER               0x01
#define CPU_CLUSTER_JUICER_SELF_LIGATION 0x02
#define CPU_CLUSTER_FORMAT_CHIASIG_STR   "chiasig"
#define CPU_CLUSTER_FORMAT_BEDPE_STR     "bedpe"
#define CPU_CLUSTER_FORMAT_INTRA         0x01
#define CPU_CLUSTER_FORMAT_INTER         0x02
#define CPU_CLUSTER_FORMAT_SEPARATE      0x04
#define CPU_CLUSTER_FORMAT_CHIASIG       0x10
#define CPU_CLUSTER_FORMAT_BEDPE         0x20
#define CPU_CLUSTER_FORMAT_MASK          0xF0

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

#define FUNC_CALL_COUNTING
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
#ifdef FUNC_CALL_COUNTING
	uint32_t clusterCalls;
	uint32_t mergeCalls;
#endif
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
	int32_t merged;  // -1 if there is no merging, otherwise 0..n = cluster id to merged into
	int32_t level;  // the depth level in the dependency graph
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
	int midPoint;
	
	int disponoParallel; // TODO: we should have bit based flags!
	int outputSortedList;
	int cycle; // TODO: we should have bit based flags!
	int iterationProgressBlock;
	
	int bulkiPETsTrigger;
	int bulkiPETsBlock;
	
	kstring_t outputPrefix;
	
	int juicerOutput;
	int outputFormat;
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
	o->midPoint = 0;
	
	o->disponoParallel = 0;
	o->outputSortedList = 0;
	o->cycle = -1;
	o->iterationProgressBlock = SUBCLUSTER_ITERACTION_BLOCK;
	
	o->bulkiPETsTrigger = BULK_IPET_TRIGGER;
	o->bulkiPETsBlock = BULK_IPET_BLOCK;
	
	memset(&(o->outputPrefix), 0, sizeof(kstring_t));
	
	o->juicerOutput = CPU_CLUSTER_JUICER_SELF_LIGATION;
	o->outputFormat = (CPU_CLUSTER_FORMAT_CHIASIG|CPU_CLUSTER_FORMAT_INTRA|CPU_CLUSTER_FORMAT_INTER);

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
	
	int (*computeIntraSpan)(int32_t,int32_t,int32_t,int32_t);
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
	// set up the result in w->clusterRecords[i]
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
	
	// consideration
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
			// TODO: FUTURE: we ignore sceondary alignment for now
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
							// this is considered fatal, cannot continue
							//if (bwa_verbose >= 1)
							fprintf(stderr, "[E::%s] Unknown XT tag \"%c\" for \"%s\". Make sure that mapper produce XT tag like bwa aln..\n", __func__, a_xt, bam1_qname(bam));
						}
					} else {
						// this is considered fatal, cannot continue
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
			
			clusterRec->tid = c1->tid;
			clusterRec->pos = c1->pos;
			clusterRec->strand = (BAM_FREVERSE==(c1->flag&BAM_FREVERSE)) ? 1 : 0;
			clusterRec->linkerType = ri1->linker_type & 0x0F;
			clusterRec->uniq = (CPU_MAPSTATE_UNIQUE==anchor1MapState) ? 1 : 0;
			clusterRec->rtid = 2*ri1->readId+ri1->tagId+1;
			
			clusterRec->mtid = c2->tid;
			clusterRec->mpos = c2->pos;
			clusterRec->mstrand = (BAM_FREVERSE==(c2->flag&BAM_FREVERSE)) ? 1 : 0;
			clusterRec->mlinkerType = ri2->linker_type & 0x0F;
			clusterRec->muniq = (CPU_MAPSTATE_UNIQUE==anchor2MapState) ? 1 : 0;
			clusterRec->mrtid = 2*ri2->readId+ri2->tagId+1;

			w->tagGroups[i].bamRecs[0] = &(w->bamRecs[anchor1ReadIndex]);
			w->tagGroups[i].bamRecs[1] = &(w->bamRecs[anchor2ReadIndex]);
	
			// check if this is not a self-ligation
			if (clusterRec->tid==clusterRec->mtid) {
				int32_t tlen = w->computeIntraSpan(clusterRec->pos,clusterRec->endpos,clusterRec->mpos,clusterRec->mendpos);
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
				anchor1Bam->core.flag |= (BAM_FPAIRED|BAM_FPROPER_PAIR); //FLAG, read paired
				if (BAM_FREVERSE & anchor2Bam->core.flag) anchor1Bam->core.flag |= BAM_FMREVERSE;
				anchor1Bam->core.flag |= BAM_FREAD1;
				
				anchor2Bam->core.flag |= (BAM_FPAIRED|BAM_FPROPER_PAIR); //FLAG, read paired
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

void cluster_process_alns(const cluster_opt_t *opt, uint32_t n_processed, int n, bam1_t *bamRecs, cpu_readid_t *readids, uint32_t *pNumTagGroups, tag_group_t *tagGroups, cpu_cluster_record_t *clusterRecords, CPUBamBuffer_t* bamsBuffer, uint32_t nUpdateBam, int (*computeIntraSpan)(int32_t,int32_t,int32_t,int32_t))
{
	int i;
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	if (n<=0) return;
	
	w.opt = opt;
	w.n_processed = n_processed; // TODO:
	// TODO: w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	
	w.bamRecs = bamRecs;
	w.readids = readids;
	w.nUpdateBam = nUpdateBam;
	w.computeIntraSpan = computeIntraSpan;
	
	// parse the read id
	kt_for(opt->n_threads, readid_worker, &w, n);
	
	// tag grouping
	// TODO: FUTURE: parallelize later
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

void report_iPETs_List (FILE *logFH, const cpu_cluster_record_t *clusterKeys, uint32_t nClusterKeys, const bam_header_t *header, int printSwap) {
	int32_t i=0;
	fprintf(logFH, "#pairId\tqual\tusable\tcid\tmcid\toid\tmoid\tiPETs\tlinker\tmlinker\tflinker\tclass\ttag1\ttag2\ttid\tpos\tendpos\tstrand\tmtid\tmpos\tmendpos\tmstrand\tisize");
	if (0!=printSwap) fprintf(logFH, "\tswap");
#if 0
	fprintf(logFH, "\tlineno\tmlineno");
#endif
	fprintf(logFH, "\n");

	if (NULL != clusterKeys) {
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
			
			if (0!=printSwap) {
				kputc('\t', &str); kputw(clusterRec->swap, &str);
			}
			
#if 0
			kputc('\t', &str); kputul(clusterRec->lineno, &str);
			kputc('\t', &str); kputul(clusterRec->mlineno, &str);
#endif
	
			// DEBUG: check that we are okie in term of the coordinates system!
#if 1
			if (clusterRec->tid==clusterRec->mtid && clusterRec->pos>clusterRec->mpos) {
				kputs("\t*coord*", &str);
			}
#endif
			kputc('\n', &str);
			
			fprintf(logFH, "%s", str.s);
			free(str.s);
		}
	}
}

int cmpPETTag1(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int32_t i32comp = ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	//int32_t i32comp = ia->usable - ib->usable; if (0!=i32comp) return i32comp;
	
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	
	// strand1, strand2, final linker, linker1, linker2, score
	i32comp = ia->strand - ib->strand; if (0!=i32comp) return i32comp;
	i32comp = ia->classification - ib->classification; if (0!=i32comp) return i32comp;
	i32comp = ia->finalLinkerType - ib->finalLinkerType; if (0!=i32comp) return i32comp;
	i32comp = ia->linkerType - ib->linkerType; if (0!=i32comp) return i32comp;
	// TODO: FUTURE: qual is higher first!
	//i32comp = ia->qual - ib->qual; if (0!=i32comp) return i32comp;
	i32comp = ib->qual - ia->qual; if (0!=i32comp) return i32comp;
	i32comp = ia->lineno - ib->lineno; //if (!i32comp) return i32comp;
	return i32comp;
}

int cmpPETTag2(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int32_t i32comp = ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	//int32_t i32comp = ia->usable - ib->usable; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->mendpos - ia->mendpos; if (0!=i32comp) return i32comp;
	i32comp = ib->pos - ia->pos; if (0!=i32comp) return i32comp;
	
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

	//new
	i32comp = (ib->mendpos-ib->pos) - (ia->mendpos-ia->pos); if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	
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
	
	// int32_t i32comp = ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	//int32_t i32comp = ia->usable - ib->usable; if (0!=i32comp) return i32comp;
	// sort by chrom1, start1, end1,
	int i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	
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

int cmpPETTag2Order(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	//int32_t i32comp = ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	//int32_t i32comp = ia->usable - ib->usable; if (0!=i32comp) return i32comp;
	// follow by chrom2, start2, end2
	int i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->mendpos - ia->mendpos; if (0!=i32comp) return i32comp;
	i32comp = ib->pos - ia->pos; if (0!=i32comp) return i32comp;
	
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

int cmpPETBothTagsOrder(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	//int32_t i32comp;
	//= ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	// sort by chrom1, chrom2
	//int32_t
	int i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	
	// TODO: any proxy for optimization?
	
	// new!!!
	i32comp = (ib->mendpos-ib->pos) - (ia->mendpos-ia->pos); if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	
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

int cmpPETBothTagsAfterMerger(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	
	// TODO: any proxy for optimization?
	i32comp = ib->iPET - ia->iPET; if (0!=i32comp) return i32comp;
	
	//TODO: QUESTION: is any ordering violated?
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;

	// NOTE: we do NOT need the rest of the fields comparison on
	//       strand1, strand2, final linker, linker1, linker2, score
	i32comp = ia->lineno - ib->lineno; //if (!i32comp) return i32comp;
	return i32comp;
}

// TODO: this routine can be further optimized
// TOOD: to implement code
static void tags_midpoint_then_extension_worker(void *data, int i, int tid)
{
	// -M: mid-point first, then extension
	int32_t pos5p;
	int32_t pos3p;
	int32_t center;
	
	worker_t *w = (worker_t*)data;
	cpu_cluster_record_t *clusterRecord = &(w->clusterRecords[i]);
	
	center = (clusterRecord->pos + clusterRecord->endpos) / 2;
	clusterRecord->pos = center;
	clusterRecord->endpos = center;
	
	if (5==w->opt->extension5p.source) {
		// use 5' position as reference
		if (clusterRecord->strand) pos5p = clusterRecord->endpos - w->opt->extension5p.size;
		else pos5p = clusterRecord->pos + w->opt->extension5p.size;
	} else {
		// use 3' position as reference
		if (clusterRecord->strand) pos5p = clusterRecord->pos - w->opt->extension5p.size;
		else pos5p = clusterRecord->endpos + w->opt->extension5p.size;
	}
	
	if (5==w->opt->extension3p.source) {
		// use 5' position as reference
		if (clusterRecord->strand) pos3p = clusterRecord->endpos - w->opt->extension3p.size;
		else pos3p = clusterRecord->pos + w->opt->extension3p.size;
	} else {
		// use 3' position as reference
		if (clusterRecord->strand) pos3p = clusterRecord->pos - w->opt->extension3p.size;
		else pos3p = clusterRecord->endpos + w->opt->extension3p.size;
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
	
	center = (clusterRecord->mpos + clusterRecord->mendpos) / 2;
	clusterRecord->mpos = center;
	clusterRecord->mendpos = center;
	
	if (5==w->opt->extension5p.source) {
		// use 5' position as reference
		if (clusterRecord->mstrand) pos5p = clusterRecord->mendpos - w->opt->extension5p.size;
		else pos5p = clusterRecord->mpos + w->opt->extension5p.size;
	} else {
		// use 3' position as reference
		if (clusterRecord->mstrand) pos5p = clusterRecord->mpos - w->opt->extension5p.size;
		else pos5p = clusterRecord->mendpos + w->opt->extension5p.size;
	}
	
	if (5==w->opt->extension3p.source) {
		// use 5' position as reference
		if (clusterRecord->mstrand) pos3p = clusterRecord->mendpos - w->opt->extension3p.size;
		else pos3p = clusterRecord->mpos + w->opt->extension3p.size;
	} else {
		// use 3' position as reference
		if (clusterRecord->mstrand) pos3p = clusterRecord->mpos - w->opt->extension3p.size;
		else pos3p = clusterRecord->mendpos + w->opt->extension3p.size;
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
	
	// direction of extension can violate invariant of tag1<tag2
	if (clusterRecord->mtid==clusterRecord->tid) {
		// TODO: we either allow extension to go beyond beyond the 5' of tag1 and 3' of tag2 or we do NOT
#if 0
		int swap = 0;
		if (clusterRecord->mpos<clusterRecord->pos) {
			swap = 1;
		}
		else if (clusterRecord->mpos==clusterRecord->pos && clusterRecord->mendpos<clusterRecord->endpos) {
			swap = 1;
		}
		
		if (swap) {
			int32_t value;
			value = clusterRecord->pos; clusterRecord->pos = clusterRecord->mpos; clusterRecord->mpos = value;
			value = clusterRecord->endpos; clusterRecord->endpos = clusterRecord->mendpos; clusterRecord->mendpos = value;
			value = clusterRecord->lineno; clusterRecord->lineno = clusterRecord->mlineno; clusterRecord->mlineno = value;
			
			// TODO: optimization; we should be able to use bits-mask to get the value
			value = clusterRecord->strand; clusterRecord->strand = clusterRecord->mstrand; clusterRecord->mstrand = value;
			value = clusterRecord->linkerType; clusterRecord->linkerType = clusterRecord->mlinkerType; clusterRecord->mlinkerType = value;
			value = clusterRecord->uniq; clusterRecord->uniq = clusterRecord->muniq; clusterRecord->muniq = value;
			value = clusterRecord->rtid; clusterRecord->rtid = clusterRecord->mrtid; clusterRecord->mrtid = value;
		}
#else
		if (clusterRecord->mpos<clusterRecord->pos) {
			clusterRecord->mpos = clusterRecord->pos;
		}
		if (clusterRecord->endpos>clusterRecord->mendpos) {
			clusterRecord->endpos = clusterRecord->mendpos;
		}
#endif
	}
	
#if 0
	if (bwa_verbose>=5) {
		// DEBUG: check that coordinate system is correct
		if (clusterRecord->mtid<clusterRecord->tid) {
			fprintf(stderr, "coordinates out of order!\n");
		}
		else if (clusterRecord->mtid==clusterRecord->tid) {
			if (clusterRecord->mpos<clusterRecord->pos) {
				fprintf(stderr, "coordinates out of order!\n");
			}
			else if (clusterRecord->mpos==clusterRecord->pos) {
				if (clusterRecord->mendpos<clusterRecord->endpos) {
					fprintf(stderr, "coordinates out of order!\n");
				}
			}
		}
	}
#endif
}

static void tags_extension_then_midpoint_worker(void *data, int i, int tid)
{
	// -m: extend first, the compute mid-point
	int32_t pos5p;
	int32_t pos3p;
	
	worker_t *w = (worker_t*)data;
	cpu_cluster_record_t *clusterRecord = &(w->clusterRecords[i]);

	if (5==w->opt->extension5p.source) {
		// use 5' position as reference
		if (clusterRecord->strand) pos5p = clusterRecord->endpos - w->opt->extension5p.size;
		else pos5p = clusterRecord->pos + w->opt->extension5p.size;
	} else {
		// use 3' position as reference
		if (clusterRecord->strand) pos5p = clusterRecord->pos - w->opt->extension5p.size;
		else pos5p = clusterRecord->endpos + w->opt->extension5p.size;
	}
	
	if (5==w->opt->extension3p.source) {
		// use 5' position as reference
		if (clusterRecord->strand) pos3p = clusterRecord->endpos - w->opt->extension3p.size;
		else pos3p = clusterRecord->pos + w->opt->extension3p.size;
	} else {
		// use 3' position as reference
		if (clusterRecord->strand) pos3p = clusterRecord->pos - w->opt->extension3p.size;
		else pos3p = clusterRecord->endpos + w->opt->extension3p.size;
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

	if (5==w->opt->extension5p.source) {
		// use 5' position as reference
		if (clusterRecord->mstrand) pos5p = clusterRecord->mendpos - w->opt->extension5p.size;
		else pos5p = clusterRecord->mpos + w->opt->extension5p.size;
	} else {
		// use 3' position as reference
		if (clusterRecord->mstrand) pos5p = clusterRecord->mpos - w->opt->extension5p.size;
		else pos5p = clusterRecord->mendpos + w->opt->extension5p.size;
	}
	
	if (5==w->opt->extension3p.source) {
		// use 5' position as reference
		if (clusterRecord->mstrand) pos3p = clusterRecord->mendpos - w->opt->extension3p.size;
		else pos3p = clusterRecord->mpos + w->opt->extension3p.size;
	} else {
		// use 3' position as reference
		if (clusterRecord->mstrand) pos3p = clusterRecord->mpos - w->opt->extension3p.size;
		else pos3p = clusterRecord->mendpos + w->opt->extension3p.size;
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
	
	// direction of extension can violate invariant of tag1<tag2
	if (clusterRecord->mtid==clusterRecord->tid) {
		// TODO: we either allow extension to go beyond beyond the 5' of tag1 and 3' of tag2 or we do NOT
#if 0
		int swap = 0;
		if (clusterRecord->mpos<clusterRecord->pos) {
			swap = 1;
		}
		else if (clusterRecord->mpos==clusterRecord->pos && clusterRecord->mendpos<clusterRecord->endpos) {
			swap = 1;
		}
		
		if (swap) {
			int32_t value;
			value = clusterRecord->pos; clusterRecord->pos = clusterRecord->mpos; clusterRecord->mpos = value;
			value = clusterRecord->endpos; clusterRecord->endpos = clusterRecord->mendpos; clusterRecord->mendpos = value;
			value = clusterRecord->lineno; clusterRecord->lineno = clusterRecord->mlineno; clusterRecord->mlineno = value;
			
			// TODO: optimization; we should be able to use bits-mask to get the value
			value = clusterRecord->strand; clusterRecord->strand = clusterRecord->mstrand; clusterRecord->mstrand = value;
			value = clusterRecord->linkerType; clusterRecord->linkerType = clusterRecord->mlinkerType; clusterRecord->mlinkerType = value;
			value = clusterRecord->uniq; clusterRecord->uniq = clusterRecord->muniq; clusterRecord->muniq = value;
			value = clusterRecord->rtid; clusterRecord->rtid = clusterRecord->mrtid; clusterRecord->mrtid = value;
		}
#else
		if (clusterRecord->mpos<clusterRecord->pos) {
			// tag 2 extended beyond 5' of tag 1
			clusterRecord->mpos=clusterRecord->pos;
		}
		if (clusterRecord->endpos>clusterRecord->mendpos) {
			// tag 1 extended beyond 3' of tag2
			clusterRecord->endpos=clusterRecord->mendpos;
		}
#endif
	}

#if 0
	if (bwa_verbose>=5) {
		// DEBUG: check that coordinate system is correct
		if (clusterRecord->mtid<clusterRecord->tid) {
			fprintf(stderr, "coordinates out of order!\n");
		}
		else if (clusterRecord->mtid==clusterRecord->tid) {
			if (clusterRecord->mpos<clusterRecord->pos) {
				fprintf(stderr, "coordinates out of order!\n");
			}
			else if (clusterRecord->mpos==clusterRecord->pos) {
				if (clusterRecord->mendpos<clusterRecord->endpos) {
					fprintf(stderr, "coordinates out of order!\n");
				}
			}
		}
	}
#endif
}

void extend_PETs(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys, uint32_t *targetLens, int32_t nTargets, void (*extension_worker)(void*,int,int)) {
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	w.opt = opt;
	w.clusterRecords = clusterRecords;
	w.nTargets = nTargets;
	w.targetLens = targetLens;
	kt_for(opt->n_threads, extension_worker, &w, nClusterKeys);
}

void cluster_PETags1(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	if (nClusterKeys>0) {
		// TODO: FUTURE: use serial approach to set up the correct cid
		int32_t i=0;
		int32_t cid = 0; // cluster id
		int32_t s=clusterRecords[i].pos;
		int32_t e=clusterRecords[i].endpos;
		clusterRecords[i].oid = i;
		clusterRecords[i].cid = 0;
		for(i=1; i<nClusterKeys; ++i) {
			clusterRecords[i].oid = i;
			if (-1==clusterRecords[i].tid || 0==clusterRecords[i].usable) {
				// not a tag-pair that we can use!
				clusterRecords[i].cid = -1;
				s=-1;
				e=-1;
			} else {
				if (clusterRecords[i-1].tid==clusterRecords[i].tid) {
					if (is_overlap(s, e, clusterRecords[i].pos, clusterRecords[i].endpos)) {
						clusterRecords[i].cid = 1;
						if (clusterRecords[i].pos<s) s=clusterRecords[i].pos;
						if (clusterRecords[i].endpos>e) e=clusterRecords[i].endpos;
					} else {
						// new cluster grouping
						// clusterRecords[i].cid = 0;
						s=clusterRecords[i].pos;
						e=clusterRecords[i].endpos;
					}
				} else {
					// new cluster grouping
					// clusterRecords[i].cid = 0;
					s=clusterRecords[i].pos;
					e=clusterRecords[i].endpos;
				}
			}
		}
		cid = 0; // cluster id
		for(i=0; i<nClusterKeys; ++i) {
			if (0==clusterRecords[i].cid) {
				cid = clusterRecords[i].oid;
				clusterRecords[i].cid = cid;
			} else if (1==clusterRecords[i].cid) {
				clusterRecords[i].cid = cid;
			}
		}
	}
}

void cluster_PETags2(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	if (nClusterKeys>0) {
		// TODO: FUTURE: use serial approach to set up the correct cid
		int32_t i=0;
		int32_t mcid = 0; // cluster id
		int32_t s=clusterRecords[i].mpos;
		int32_t e=clusterRecords[i].mendpos;
		clusterRecords[i].moid = i;
		clusterRecords[i].mcid = 0;
		for(i=1; i<nClusterKeys; ++i) {
			clusterRecords[i].moid = i;
			if (-1==clusterRecords[i].mtid || 0==clusterRecords[i].usable) {
				// not a tag-pair that we can use!
				clusterRecords[i].mcid = -1;
				s=-1;
				e=-1;
			} else {
				if (clusterRecords[i-1].mtid==clusterRecords[i].mtid) {
					if (is_overlap(s, e, clusterRecords[i].mpos, clusterRecords[i].mendpos)) {
						clusterRecords[i].mcid = 1;
						if (clusterRecords[i].mpos<s) s=clusterRecords[i].mpos;
						if (clusterRecords[i].mendpos>e) e=clusterRecords[i].mendpos;
					} else {
						// new cluster grouping
						// clusterRecords[i].mcid = 0;
						s=clusterRecords[i].mpos;
						e=clusterRecords[i].mendpos;
					}
				} else {
					// new cluster grouping
					// clusterRecords[i].mcid = 0;
					s=clusterRecords[i].mpos;
					e=clusterRecords[i].mendpos;
				}
			}
		}
		mcid = 0; // cluster id
		for(i=0; i<nClusterKeys; ++i) {
			if (0==clusterRecords[i].mcid) {
				mcid = clusterRecords[i].moid;
				clusterRecords[i].mcid = mcid;
			} else if (1==clusterRecords[i].mcid) {
				clusterRecords[i].mcid = mcid;
			}
		}
	}
}

int cmp_int32_t_desc (const void * a, const void * b)
{
	return ( *(int*)b - *(int*)a );
}

void report_iPETs_Distribution(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys, uint32_t numClusterBins) {
	if (bwa_verbose>=3) {
		if (nClusterKeys>1) {
			// allocate the memory for the count, which we will sort and iterate sequentially
			double t_start = realtime();
			uint32_t nClusterBins = 0;
			int32_t *iPETs = (int32_t *) calloc(numClusterBins, sizeof(int32_t));
			if (!iPETs) {
				fprintf(stderr, "[E::%s] out of memory for %u iPETs bin counts e: %s\n", __func__, numClusterBins, strerror (errno));
				return;
			}
			
			int32_t i;
			for(i=0; i<nClusterKeys; ++i) {
				if (clusterRecords[i].iPET>0) {
					iPETs[nClusterBins] = clusterRecords[i].iPET;
					nClusterBins++;
				}
			}
			qsort(iPETs, numClusterBins, sizeof(int32_t), cmp_int32_t_desc);
			
			double t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s] Computed %u iPETs distribution across %u cluster bins, processing %.2f sec (%.2f min)..\n", __func__, nClusterKeys, numClusterBins, t_diff, t_diff/60.0);

			// report log scale when 3==bwa_verbose and detail scale when 4=<bwa_verbose
			if (bwa_verbose>=4) {
				fprintf(stderr, "[M::%s] #iPETinBin\t#clusterBin\n", __func__);
				int32_t count = 0;
				int32_t iPET = iPETs[0];
				for(i=0; i<numClusterBins; ++i) {
					if (iPETs[i]==iPET) {
						count++;
					} else {
						if (0!=iPET) fprintf(stderr, "[M::%s] %d\t%d\n", __func__, iPET, count);
						count++;
						iPET = iPETs[i];
					}
				}
				free(iPETs);
				if (0!=iPET) fprintf(stderr, "[M::%s] %d\t%d\n", __func__, iPET, count);
			} else {
				// range version
				// 1,2,3,4,5,6,7,8,9,10
				// 11-50, 51-100
				// 101-500, 501-1000
				// 1001-5000, 5001-10000
				// 10001-50000, 50001-100000
				// 100001-500000, 500001-1000000
				int32_t ranges[20][2] = {
					{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,9},{10,10},
					{11,50},{51,100},{101,500},{501,1000},{1001,5000},{5001,10000},
					{10001,50000},{50001,100000},{100001,500000},{500001,1000000}
				};
				fprintf(stderr, "[M::%s] #iPETinBin\t#clusterBin\n", __func__);
				int32_t count = 0;
				int32_t lastRangeIndex = 19;
				int32_t rangeIndex = lastRangeIndex;
				if (iPETs[0]<=ranges[rangeIndex][1]) {
					for(rangeIndex = lastRangeIndex; rangeIndex>=0; --rangeIndex) {
						if (ranges[rangeIndex][0]<=iPETs[0] && iPETs[0]<=ranges[rangeIndex][1]) {
							break;
						}
					}
				} else {
					ranges[rangeIndex][1] = iPETs[0];
				}
				for(i=0; i<numClusterBins; ++i) {
					if (ranges[rangeIndex][0]<=iPETs[i] && iPETs[i]<=ranges[rangeIndex][1]) {
						// same range
						count++;
					} else {
						// different range
						if (rangeIndex>=0) fprintf(stderr, "[M::%s] [%d,%d]\t%d\n", __func__, ranges[rangeIndex][0], ranges[rangeIndex][1], count);
						count++;
						// set up the new range
						for(; rangeIndex>=0; --rangeIndex) {
							if (ranges[rangeIndex][0]<=iPETs[i] && iPETs[i]<=ranges[rangeIndex][1]) {
								break;
							}
						}
					}
				}
				free(iPETs);
				if (rangeIndex>=0) fprintf(stderr, "[M::%s] [%d,%d]\t%d\n", __func__, ranges[rangeIndex][0], ranges[rangeIndex][1], count);
			}
		}
	}
}

void setup_cluster_bin(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys, uint32_t *numClusterBins) {
	int32_t i;

	// use serial approach to set up the correct cid
	// we also overloaded ".swap" to indicate if the cluster need to be checked
	// cluster of size 1 need NO checking, i.e. 0==.swap, 1 otherwise
	*numClusterBins = 0;
	if (nClusterKeys>1) {
		uint32_t nClusterBins = 0;
		int32_t clusterIndex=0;
		int32_t cid=-1; // cluster id
		int32_t mcid=-1; // cluster id
		for(i=0; i<nClusterKeys; ++i) {
			if (cid==clusterRecords[i].cid && mcid==clusterRecords[i].mcid) {
				clusterRecords[clusterIndex].iPET++;
				
				// get new two dimensional ids
				clusterRecords[i].cid = clusterRecords[clusterIndex].oid;
				clusterRecords[i].mcid = clusterRecords[clusterIndex].moid;
				
				clusterRecords[clusterIndex].swap = 1;
				
			} else {
				clusterRecords[i].iPET = 1;
				cid=clusterRecords[i].cid;
				mcid=clusterRecords[i].mcid;
				clusterIndex = i;

				// get new two dimensional ids
				clusterRecords[i].cid = clusterRecords[i].oid;
				clusterRecords[i].mcid = clusterRecords[i].moid;
				
				// we assume it is only a single record, nothing to check
				clusterRecords[i].swap = 0;
				
				nClusterBins++;
			}
		}
		*numClusterBins = nClusterBins;
	}
}

void report_Total_iPET(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
	if (bwa_verbose >= 3) {
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
	
		fprintf(stderr, "[M::%s] Total %d iPETs, %d usable iPETs, %d iPET colsum, %d self-ligated iPETs..\n", __func__,
				nClusterKeys, usableRows, totaliPET, nClusterKeys-usableRows);
	}
}

int32_t cluster_bin(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys, int iteration, long pairId, uint32_t *pNumFuncCall) {
	int32_t i,j;
	int32_t nDiscordantPETs = 0;
	
	*pNumFuncCall+=1;
	if (nClusterKeys>1) {
		// TODO : optimization: we know where the usable end!!!
		for(i=0; i<nClusterKeys; ++i) {
			
			// iPET will indicate the number of memebers in the bin
			//if (clusterRecords[i].iPET>1) {
			if (clusterRecords[i].swap) {
				// if >=2 members, we will need to check if the overlap is proper
				int32_t niPET = clusterRecords[i].iPET;
				
				// TODO: FOR DEBUGGING ONLY
				int nBinDiscordantPETs = 0;
				
				// let's sort based on tag right and process them
				/*if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] initial\n", __func__);
					report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
				}*/
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
				/*if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] Tag2 ordered\n", __func__);
					report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
				}*/
				
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
				/*if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] Tag1 ordered\n", __func__);
					report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
				}*/
				
				// let's sort based on combine tag
				qsort(&(clusterRecords[i]), niPET, sizeof(cpu_cluster_record_t), cmpPETBothTagsOrder);
				// let's process them
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
				if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] Combined Tag1 and Tag2 ordered\n", __func__);
					report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
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
								if (bwa_verbose>=5) {
									fprintf(stderr, "[D:%s:non_overlap_tag1] i=%d, j=%d, iPET=%d, leader=%d, s1=%d(%d), e1=%d(%d), s2=%d, e2=%d, oid=%d(%d), moid=%d(%d)\n",
											__func__,
											i, j, niPET, clusterIndex, s1, clusterRecord->pos, e1, clusterRecord->endpos, s2, e2,
											clusterRecords[clusterIndex].oid, clusterRecord->oid,
											clusterRecords[clusterIndex].moid, clusterRecord->moid);
								}
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
									if (bwa_verbose>=5) {
										fprintf(stderr, "[W:%s:non_overlap_tag2] i=%d, j=%d, iPET=%d, leader=%d, s1=%d, e1=%d, s2=%d(%d), e2=%d(%d), oid=%d(%d), moid=%d(%d)\n",
												__func__,
												i, j, niPET, clusterIndex, s1, e1, s2, clusterRecord->mpos, e2, clusterRecord->mendpos,
												clusterRecords[clusterIndex].oid, clusterRecord->oid,
												clusterRecords[clusterIndex].moid, clusterRecord->moid);
									}
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
					
					if (bwa_verbose>=5) {
						fprintf(stderr, "[D:%s] After checking correction\n", __func__);
						report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
					}
					
					// we need to re-order the array
					if (nBinDiscordantPETs>0) {
						qsort(&(clusterRecords[i]), niPET, sizeof(cpu_cluster_record_t), cmpPETBothTagsOrder);
						
						if (bwa_verbose>=5) {
							fprintf(stderr, "[D:%s] Final re-ordering\n", __func__);
							report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
						}
					}
					
					// TODO: check if recursion is a good way here as we might be unable multi-thread as a result
					if (nBinDiscordantPETs>0) {
						if (bwa_verbose>=4) {
							if (0==(iteration%opt->iterationProgressBlock))
								fprintf(stderr, "[D:%s] Calling nested subclustering on %ld iteration %d with discordantIndex=%d, %d discordant iPETs\n",
										__func__, pairId, iteration, discordantIndex, nBinDiscordantPETs);
						}
							
						if (bwa_verbose>=5) {
							fprintf(stderr, "[D:%s] Calling nested subclustering with discordantIndex=%d, %d discordant iPETs\n",
									__func__, discordantIndex, nBinDiscordantPETs);
							report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
						}
						int32_t nRetryDiscordantPETs = cluster_bin(opt, &(clusterRecords[i]), niPET, iteration+1, pairId, pNumFuncCall);
						if (0==nRetryDiscordantPETs) {
							// there is no discordant set and thus we attempt to merge if we can
						}
						if (bwa_verbose>=4) {
							if (0==(iteration%opt->iterationProgressBlock))
								fprintf(stderr, "[D:%s] After nested subclustering on %ld iteration %d with with discordantIndex=%d, %d discordant iPETs, %d discordard after trying\n", __func__, pairId, iteration, discordantIndex, nBinDiscordantPETs, nRetryDiscordantPETs);
						}
						if (bwa_verbose>=5) {
							fprintf(stderr, "[D:%s] After nested subclustering with with discordantIndex=%d, %d discordant iPETs, %d discordard after trying\n", __func__, discordantIndex, nBinDiscordantPETs, nRetryDiscordantPETs);
							report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
						}
					}
				}
				
				if (bwa_verbose>=5) {
					// TODO: for DEBUGGING; check that we have not lost the iPET count in each of the original bin group
					int32_t niPETinBin = 0;
					for(j=0; j<niPET; ++j) {
						cpu_cluster_record_t *clusterRecord = &(clusterRecords[i+j]);
						niPETinBin += clusterRecord->iPET;
					}
					if (niPET!=niPETinBin) {
						fprintf(stderr, "[D:%s:iPET] i=%d, iPET=%d, iPETinBin=%d\n", __func__, i, niPET, niPETinBin);
						fprintf(stderr, "[D:%s:iPET] i=%d, iPETs=", __func__, i);
						for(j=0; j<niPET; ++j) {
							fprintf(stderr, "%d,", clusterRecords[i+j].iPET);
						}
						fprintf(stderr, "\n");
					}
				}
				
				
				i += (niPET-1);
				
				nDiscordantPETs += nBinDiscordantPETs;
			}
		}
	}
	
	return nDiscordantPETs;
}

int cmpSubclusterTag1(const void *a, const void *b) {
	const cpu_subcluster_record_t *ia = (const cpu_subcluster_record_t *)a;
	const cpu_subcluster_record_t *ib = (const cpu_subcluster_record_t *)b;
	
	// sort by start1, end1,
	int i32comp;
	i32comp = ia->merged - ib->merged; if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->endpos - ib->endpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->lineno - ia->lineno; //if (!i32comp) return i32comp;
	
	return i32comp;
}

int cmpSubclusterTag2(const void *a, const void *b) {
	const cpu_subcluster_record_t *ia = (const cpu_subcluster_record_t *)a;
	const cpu_subcluster_record_t *ib = (const cpu_subcluster_record_t *)b;
	
	int i32comp;
	i32comp = ia->merged - ib->merged; if (0!=i32comp) return i32comp;

	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
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
	i32comp = ia->merged - ib->merged; if (0!=i32comp) return i32comp;
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

int cmpSubclusterFlattenDependency(const void *a, const void *b) {
	const cpu_subcluster_record_t *ia = (const cpu_subcluster_record_t *)a;
	const cpu_subcluster_record_t *ib = (const cpu_subcluster_record_t *)b;
	
	int32_t i32comp;
	i32comp = ib->level - ia->level; if (0!=i32comp) return i32comp;
	i32comp = ia->merged - ib->merged; if (0!=i32comp) return i32comp;
	i32comp = ia->offset - ib->offset; //if (0!=i32comp) return i32comp;
	return i32comp;
}

int cmpSubclusterUpdateOrderByLevel(const void *a, const void *b) {
	const cpu_subcluster_record_t *ia = (const cpu_subcluster_record_t *)a;
	const cpu_subcluster_record_t *ib = (const cpu_subcluster_record_t *)b;
	
	int32_t i32comp;
	i32comp = ia->level - ib->level; if (0!=i32comp) return i32comp;
	i32comp = ib->merged - ia->merged; if (0!=i32comp) return i32comp;
	i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	i32comp = ia->lineno - ia->lineno; //if (!i32comp) return i32comp;
	return i32comp;
}

int cmpSubclusterUpdateOrder(const void *a, const void *b) {
	const cpu_subcluster_record_t *ia = (const cpu_subcluster_record_t *)a;
	const cpu_subcluster_record_t *ib = (const cpu_subcluster_record_t *)b;
	
	int32_t i32comp;
	i32comp = ib->merged - ia->merged; if (0!=i32comp) return i32comp;
	i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	
	// TODO: any proxy for optimization?
	
	//new
	i32comp = (ib->mendpos-ib->pos) - (ia->mendpos-ia->pos); if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	
	i32comp = ia->lineno - ia->lineno; //if (!i32comp) return i32comp;
	
	return i32comp;
}

void report_Reads_Juicer (gzFile *juicerfd, const cpu_cluster_record_t *clusterKeys, uint32_t nClusterKeys, const bam_header_t *header, uint32_t excludeSelfLigation) {
	// Medium format
	// <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>
	// For read, we do not need the <score> 'cos it is always 1
	// NOTE: do not report clusterRec->iPET as it is zero to start with; NOT clustered yet.
	// str = strand (0 for forward, anything else for reverse)
	// chr = chromosome (must be a chromosome in the genome)
	// pos = position
	// frag = restriction site fragment
	// score = the score imputed to this read

	// short chr format!!!
	// check self-ligation
	int32_t nOffset = 0;
	if (header) {
		// TODO: should have a better way to determine this!!!
		if ('c'==header->target_name[0][0] || 'C'==header->target_name[0][0]) {
			if ('h'==header->target_name[0][1] || 'H'==header->target_name[0][1]) {
				if ('r'==header->target_name[0][2] || 'R'==header->target_name[0][2]) {
					nOffset = 3;
				}
			}
		}
	}
	if (0==excludeSelfLigation) {
		double t_diffProcessing;
		double t_timepointIO, t_timepointProcessing;
		t_timepointProcessing = realtime();
		uint32_t i;
		for(i=0; i<nClusterKeys; ++i) {
			const cpu_cluster_record_t *clusterRec = &(clusterKeys[i]);
			gzprintf(*juicerfd, "%d %s %d 0 %d %s %d 1\n",
					 clusterRec->strand, header->target_name[clusterRec->tid]+nOffset,
					 (0==clusterRec->strand) ? clusterRec->pos : clusterRec->endpos,
					 clusterRec->mstrand, header->target_name[clusterRec->mtid]+nOffset,
					 (0==clusterRec->mstrand) ? clusterRec->mpos : clusterRec->mendpos);
		}
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Written %u iPET records for juicer, with self-ligation iPETs. %.2f sec ( %.2f min)..\n", __func__, i, t_diffProcessing, t_diffProcessing/60.0);
		}
	} else {
		double t_diffProcessing;
		double t_timepointIO, t_timepointProcessing;
		t_timepointProcessing = realtime();
		uint32_t i;
		for(i=0; i<nClusterKeys; ++i) {
			const cpu_cluster_record_t *clusterRec = &(clusterKeys[i]);
			if (1==clusterRec->usable) {
				gzprintf(*juicerfd, "%d %s %d 0 %d %s %d 1\n",
						 clusterRec->strand, header->target_name[clusterRec->tid]+nOffset,
						 (0==clusterRec->strand) ? clusterRec->pos : clusterRec->endpos,
						 clusterRec->mstrand, header->target_name[clusterRec->mtid]+nOffset,
						 (0==clusterRec->mstrand) ? clusterRec->mpos : clusterRec->mendpos);
			} else {
				// self-ligation start only after non-self-ligation
				break;
			}
		}
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] Written %u iPET records for juicer, where self-ligation iPETs are excluded. %.2f sec ( %.2f min)..\n", __func__, i, t_diffProcessing, t_diffProcessing/60.0);
		}
	}
}

void report_Clusters_List_CHIASIG (gzFile *clustersfd, gzFile *clustersTransfd, const cpu_cluster_record_t *clusterKeys, uint32_t nClusterKeys, const bam_header_t *header) {
	// <chrom1> <start1> <end1> <chrom2> <start2> <end2> <score>
	if (NULL==*clustersTransfd) {
		int32_t i;
		int32_t cid=-1; // cluster id
		int32_t mcid=-1; // cluster id
		for(i=0; i<nClusterKeys; ++i) {
			const cpu_cluster_record_t *clusterRec = &(clusterKeys[i]);
			if (1==clusterRec->usable) {
				gzprintf(*clustersfd, "%s\t%d\t%d\t%s\t%d\t%d\t%d\n",
						 header->target_name[clusterRec->tid], clusterRec->pos, clusterRec->endpos,
						 header->target_name[clusterRec->mtid], clusterRec->mpos, clusterRec->mendpos,
						 clusterRec->iPET);
				if (cid==clusterRec->cid && mcid==clusterRec->mcid) {
					fprintf(stderr, "[E::report_Clusters_List_CHIASIG] Still the same cluster after advancement, i=%d, cid=%d, mcid=%d, oid=%d, moid=%d, pairId=%lu\n",
							i, clusterRec->cid, clusterRec->mcid, clusterRec->oid, clusterRec->moid, clusterRec->pairId);
				}
				cid = clusterRec->cid;
				mcid = clusterRec->mcid;
				i += ((clusterRec->iPET>0?clusterRec->iPET : 1)-1);
			}
		}
	} else {
		// we have to write trans and cis seprately!!!
		int32_t i;
		int32_t cid=-1; // cluster id
		int32_t mcid=-1; // cluster id
		for(i=0; i<nClusterKeys; ++i) {
			const cpu_cluster_record_t *clusterRec = &(clusterKeys[i]);
			if (1==clusterRec->usable) {
				if (clusterRec->tid == clusterRec->mtid) {
					gzprintf(*clustersfd, "%s\t%d\t%d\t%s\t%d\t%d\t%d\n",
							 header->target_name[clusterRec->tid], clusterRec->pos, clusterRec->endpos,
							 header->target_name[clusterRec->mtid], clusterRec->mpos, clusterRec->mendpos,
							 clusterRec->iPET);
				} else {
					gzprintf(*clustersTransfd, "%s\t%d\t%d\t%s\t%d\t%d\t%d\n",
							 header->target_name[clusterRec->tid], clusterRec->pos, clusterRec->endpos,
							 header->target_name[clusterRec->mtid], clusterRec->mpos, clusterRec->mendpos,
							 clusterRec->iPET);
				}
				
				if (cid==clusterRec->cid && mcid==clusterRec->mcid) {
					fprintf(stderr, "[E::report_Clusters_List_CHIASIG] Still the same cluster after advancement, i=%d, cid=%d, mcid=%d, oid=%d, moid=%d, pairId=%lu\n",
							i, clusterRec->cid, clusterRec->mcid, clusterRec->oid, clusterRec->moid, clusterRec->pairId);
				}
				cid = clusterRec->cid;
				mcid = clusterRec->mcid;
				i += ((clusterRec->iPET>0?clusterRec->iPET : 1)-1);
			}
		}
	}
}

void report_Clusters_List_BEDPE (gzFile *clustersfd, gzFile *clustersTransfd, const cpu_cluster_record_t *clusterKeys, uint32_t nClusterKeys, const bam_header_t *header) {
	// <chrom1> <start1> <end1> <chrom2> <start2> <end2> <name> <score> <strand1> <strand2> <user-defined11> <user-defined12> ...
	if (NULL==*clustersTransfd) {
		int32_t i;
		int32_t cid=-1; // cluster id
		int32_t mcid=-1; // cluster id
		for(i=0; i<nClusterKeys; ++i) {
			const cpu_cluster_record_t *clusterRec = &(clusterKeys[i]);
			if (1==clusterRec->usable) {
				gzprintf(*clustersfd, "%s\t%d\t%d\t%s\t%d\t%d\t%c_%d_%d\t%d\n",
						 header->target_name[clusterRec->tid], clusterRec->pos, clusterRec->endpos,
						 header->target_name[clusterRec->mtid], clusterRec->mpos, clusterRec->mendpos,
						 (clusterRec->tid == clusterRec->mtid) ? 'C' : 'T', clusterRec->cid, clusterRec->mcid,
						 clusterRec->iPET);
				if (cid==clusterRec->cid && mcid==clusterRec->mcid) {
					fprintf(stderr, "[E::report_Clusters_List_BEDPE] Still the same cluster after advancement, i=%d, cid=%d, mcid=%d, oid=%d, moid=%d, pairId=%lu\n",
							i, clusterRec->cid, clusterRec->mcid, clusterRec->oid, clusterRec->moid, clusterRec->pairId);
				}
				cid = clusterRec->cid;
				mcid = clusterRec->mcid;
				i += ((clusterRec->iPET>0?clusterRec->iPET : 1)-1);
			}
		}
	} else {
		// we have to write trans and cis seprately!!!
		int32_t i;
		int32_t cid=-1; // cluster id
		int32_t mcid=-1; // cluster id
		for(i=0; i<nClusterKeys; ++i) {
			const cpu_cluster_record_t *clusterRec = &(clusterKeys[i]);
			if (1==clusterRec->usable) {
				if (clusterRec->tid == clusterRec->mtid) {
					gzprintf(*clustersfd, "%s\t%d\t%d\t%s\t%d\t%d\tC_%d_%d\t%d\n",
							 header->target_name[clusterRec->tid], clusterRec->pos, clusterRec->endpos,
							 header->target_name[clusterRec->mtid], clusterRec->mpos, clusterRec->mendpos,
							 clusterRec->cid, clusterRec->mcid,
							 clusterRec->iPET);
				} else {
					gzprintf(*clustersfd, "%s\t%d\t%d\t%s\t%d\t%d\tT_%d_%d\t%d\n",
							 header->target_name[clusterRec->tid], clusterRec->pos, clusterRec->endpos,
							 header->target_name[clusterRec->mtid], clusterRec->mpos, clusterRec->mendpos,
							 clusterRec->cid, clusterRec->mcid,
							 clusterRec->iPET);
				}
				if (cid==clusterRec->cid && mcid==clusterRec->mcid) {
					fprintf(stderr, "[E::report_Clusters_List_BEDPE] Still the same cluster after advancement, i=%d, cid=%d, mcid=%d, oid=%d, moid=%d, pairId=%lu\n",
							i, clusterRec->cid, clusterRec->mcid, clusterRec->oid, clusterRec->moid, clusterRec->pairId);
				}
				cid = clusterRec->cid;
				mcid = clusterRec->mcid;
				i += ((clusterRec->iPET>0?clusterRec->iPET : 1)-1);
			}
		}
	}
}

void report_Subclusters_List (FILE *clustersfd, const cpu_subcluster_record_t *subclusters, uint32_t nSubclusters, int32_t *subclustersLevel, const bam_header_t *header) {
	int32_t i;
	
	fprintf(stderr, "#pairId\toffset\tmerged\tlevel\tcid\tmcid\toid\tmoid\tiPETs\tpos\tendpos\tmpos\tmendpos\n");
	
	for(i=0; i<nSubclusters; ++i) {
		
		const cpu_subcluster_record_t *subcluster = &(subclusters[i]);
#if 1
		// TODO: we need to report two different level information; subcluster->level, subclustersLevel[]
		fprintf(stderr, "%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				subcluster->pairId, subcluster->offset, subcluster->merged, subclustersLevel ? subclustersLevel[subcluster->offset] : subcluster->level,
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
		kputc('\t', &str); kputl(subcluster->merged, &str);
		kputc('\t', &str); kputl(subcluster->level, &str);
		
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

int32_t merge_subclusters (const cluster_opt_t *opt, cpu_subcluster_record_t *subclusters, int32_t *subclustersLevel, int32_t *pLevel, uint32_t nSubclusters, uint32_t nFinalSubclusters, uint32_t *pNumFuncCall) {
	int32_t numSubclusters = nSubclusters;
	int32_t j;
	
	FILE *logFH	= stderr;
	//FILE *logFH	= stdout;
	
	*pNumFuncCall+=1;
	
	if (nSubclusters>1) {
		// TODO: FOR DEBUGGING ONLY
		
		// let's sort based on tag right and process them
		/*if (bwa_verbose >= 5) {
			fprintf(logFH, "[D:%s] initial\n", __func__);
			report_Subclusters_List (logFH, subclusters, nSubclusters, subclustersLevel, NULL);
		}*/
		
		qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterTag2);
		// let's process sorted right tags
		//cluster_PETags2(opt, clusterKeys, nClusterKeys);
		int32_t s=subclusters[0].mpos;
		int32_t e=subclusters[0].mendpos;
		subclusters[0].mcid = 0;
		for(j=1; j<nSubclusters; ++j) {
			//cpu_subcluster_record_t *pSubcluster = &(subclusters[j-1]);
			cpu_subcluster_record_t *subcluster = &(subclusters[j]);
			//if (is_overlap(pSubcluster->mpos, pSubcluster->mendpos, subcluster->mpos, subcluster->mendpos)) {
			if (is_overlap(s, e, subcluster->mpos, subcluster->mendpos)) {
				subcluster->mcid = 1;
				if (subcluster->mpos<s) s=subcluster->mpos;
				if (subcluster->mendpos>e) e=subcluster->mendpos;
			} else {
				subcluster->mcid = 0;
				s=subcluster->mpos;
				e=subcluster->mendpos;
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
		
		/*if (bwa_verbose >= 5) {
			fprintf(logFH, "[D:%s] Tag2 ordered\n", __func__);
			report_Subclusters_List (logFH, subclusters, nSubclusters, subclustersLevel, NULL);
		}*/
		
		// let's sort based on tag left and process them
		qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterTag1);
		// let's process sorted left tags
		//cluster_PETags1(opt, clusterKeys, nClusterKeys);
		s=subclusters[0].pos;
		e=subclusters[0].endpos;
		subclusters[0].cid = 0;
		for(j=1; j<nSubclusters; ++j) {
			//cpu_subcluster_record_t *pSubcluster = &(subclusters[j-1]);
			cpu_subcluster_record_t *subcluster = &(subclusters[j]);
			//if (is_overlap(pSubcluster->pos, pSubcluster->endpos, subcluster->pos, subcluster->endpos)) {
			if (is_overlap(s, e, subcluster->pos, subcluster->endpos)) {
				subcluster->cid = 1;
				if (subcluster->pos<s) s=subcluster->pos;
				if (subcluster->endpos>e) e=subcluster->endpos;
			} else {
				subcluster->cid = 0;
				s=subcluster->pos;
				e=subcluster->endpos;
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
		
		/*if (bwa_verbose >= 5) {
			fprintf(logFH, "[D:%s] Tag1 ordered\n", __func__);
			report_Subclusters_List (logFH, subclusters, nSubclusters, subclustersLevel, NULL);
		}*/
		
		// let's sort based on combine tag
		qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterBothTags);
		// let's process them
		//s=subclusters[0].mpos;
		//e=subclusters[0].mendpos;
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
		
		if (bwa_verbose >= 5) {
			fprintf(logFH, "[D:%s] Combined Tag1 and Tag2 ordered, estimated %d subclusters\n", __func__, numSubclusters);
			report_Subclusters_List (logFH, subclusters, nSubclusters, subclustersLevel, NULL);
		}
		
		// TODO: check condition of merging and update the ranges
		*pLevel = *pLevel + 1;
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
						if (bwa_verbose >= 5) {
							fprintf(logFH, "[D:%s:tag1] j=%d, iPET=%d, leader=%d, s1=%d(%d), e1=%d(%d), s2=%d, e2=%d, oid=%d(%d), moid=%d(%d)\n",
									__func__, j, subcluster->iPET, subclusterIndex, s1, subcluster->pos, e1, subcluster->endpos, s2, e2,
									subclusters[subclusterIndex].oid, subcluster->oid,
									subclusters[subclusterIndex].moid, subcluster->moid);
						}
						//nMergedSubclusters++;
						
						//subcluster[subclusterIndex].iPET--;
						if (-1==discordantIndex) {
							discordantIndex = j;
							subcluster->cid = subcluster->oid;
							subcluster->mcid = subcluster->moid;
							//subclusters[discordantIndex].iPET = 1;
							//subclusters[discordantIndex].swap = 0;
							
							/*WCH*///
							subclustersLevel[subcluster->offset] = *pLevel;
							
							subcluster->merged = -1;
						} else {
							subcluster->cid = subclusters[discordantIndex].oid;
							subcluster->mcid = subclusters[discordantIndex].moid;
							//subclusters[discordantIndex].iPET++;
							//subclusters[discordantIndex].swap = 1;
							
							/*WCH*///
							subclustersLevel[subcluster->offset] = *pLevel;
							
							subclusters[discordantIndex].merged = -2;
							subcluster->merged = -2;
						}
						
					} else {
						if (0==is_overlap(s2, e2, subcluster->mpos, subcluster->mendpos)) {
							// does not overlap! this is a problem, report it
							if (bwa_verbose >= 5) {
								fprintf(logFH, "[D:%s:tag2] j=%d, iPET=%d, leader=%d, s1=%d, e1=%d, s2=%d(%d), e2=%d(%d), oid=%d(%d), moid=%d(%d)\n",
										__func__, j, subcluster->iPET, subclusterIndex, s1, e1, s2, subcluster->mpos, e2, subcluster->mendpos,
										subclusters[subclusterIndex].oid, subcluster->oid,
										subclusters[subclusterIndex].moid, subcluster->moid);
							}
							//nMergedSubclusters++;
							
							//subclusters[clusterIndex].iPET--;
							if (-1==discordantIndex) {
								discordantIndex = j;
								subcluster->cid = subcluster->oid;
								subcluster->mcid = subcluster->moid;
								//subclusters[discordantIndex].iPET = 1;
								//subclusters[discordantIndex].swap = 0;
								
								/*WCH*///
								subclustersLevel[subcluster->offset] = *pLevel;
								
								subcluster->merged = -1;
							} else {
								subcluster->cid = subclusters[discordantIndex].oid;
								subcluster->mcid = subclusters[discordantIndex].moid;
								//subclusters[discordantIndex].iPET++;
								//subclusters[discordantIndex].swap = 1;
								
								/*WCH*///
								subclustersLevel[subcluster->offset] = *pLevel;
								
								subclusters[discordantIndex].merged = -2;
								subcluster->merged = -2;
							}
							
						} else {
							// both anchors overlapped, let's expand our anchor regions
							if (subcluster->pos<s1) {
								s1 = subcluster->pos;
								subclusters[subclusterIndex].pos = s1;
							}
							if (e1<subcluster->endpos) {
								e1 = subcluster->endpos;
								subclusters[subclusterIndex].endpos = e1;
							}
							if (subcluster->mpos<s2) {
								s2 = subcluster->mpos;
								subclusters[subclusterIndex].mpos = s2;
							}
							if (e2<subcluster->mendpos) {
								e2 = subcluster->mendpos;
								subclusters[subclusterIndex].mendpos = e2;
							}
							
							subcluster->merged = subclusters[subclusterIndex].offset;
							subclusters[subclusterIndex].iPET += subcluster->iPET;
							//subcluster->iPET = 0;
							
							/*WCH*///
							subclustersLevel[subcluster->offset] = *pLevel;
							
							nMergedSubclusters++;
						}
					}
					
				} else {
					// we are starting a new bin grouping, no checking
					cid=subcluster->cid;
					mcid=subcluster->mcid;
					subclusterIndex = j;
					subcluster->merged = -1; // TODO: check if this gives us problem
					
					s1=subclusters[subclusterIndex].pos;
					e1=subclusters[subclusterIndex].endpos;
					s2=subclusters[subclusterIndex].mpos;
					e2=subclusters[subclusterIndex].mendpos;
					
					/*WCH*///
					subclustersLevel[subcluster->offset] = *pLevel;
					// at the end of this check, the group is confirmed
					// outliers will be expelled to discordantIndex set
					
					discordantIndex = -1;
				}
			}
			
			if (bwa_verbose >= 5) {
				fprintf(logFH, "[D:%s] After checking correction, with %d merger(s)\n", __func__, nMergedSubclusters);
				report_Subclusters_List (logFH, subclusters, nSubclusters, subclustersLevel, NULL);
			}
			
			int32_t pendingSubclusters = 0;
			for(j=0; j<nSubclusters; ++j) {
				if (-2==subclusters[j].merged) pendingSubclusters++;
			}
			if (pendingSubclusters>0) {
				qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterBothTags);
				
				if (bwa_verbose >= 5) {
					fprintf(logFH, "[D:%s] Pending subcluster combination\n", __func__);
					report_Subclusters_List (logFH, subclusters, nSubclusters, subclustersLevel, NULL);
				}
				int32_t final_numSubclusters = merge_subclusters (opt, subclusters, subclustersLevel, pLevel, pendingSubclusters, pendingSubclusters, pNumFuncCall);
				if (bwa_verbose >= 5) {
					fprintf(logFH, "[D:%s] Pending subcluster combination merging, triggered by %d, subclusters: %d \n", __func__, pendingSubclusters, final_numSubclusters);
				}
				// if there is no changes, we do no need to trigger anything later
				if (final_numSubclusters==pendingSubclusters) {
					pendingSubclusters = 0;
					
					// TODO: how do we know that there has been mergers and thus, should be rolled up to this level?!
					
					// TODO: check if we need to re-order!
					//qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterBothTags);
				}
			}
			
			
			// we need to re-order the array
			int32_t newNumSubclusters = 0;
			for(j=0; j<nSubclusters; ++j) {
				if (-1==subclusters[j].merged) newNumSubclusters++;
			}
			if (nMergedSubclusters>0 || pendingSubclusters>0) {
				qsort(subclusters, nSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterBothTags);
				
				if (bwa_verbose >= 5) {
					fprintf(logFH, "[D:%s] Final re-ordering\n", __func__);
					report_Subclusters_List (logFH, subclusters, nSubclusters, subclustersLevel, NULL);
				}
				
				// TODO: check condition of invocation & termination
				// TODO: decide if we need to tune the number of clusters!!
				//if (numSubclusters!=nSubclusters) {
				int32_t final_numSubclusters = merge_subclusters (opt, subclusters, subclustersLevel, pLevel, newNumSubclusters, newNumSubclusters, pNumFuncCall);
				if (bwa_verbose >= 5) {
					fprintf(logFH, "[D:%s] Re-cur merging, triggered by %d, subclusters: %d -> %d \n", __func__, nMergedSubclusters, nSubclusters, final_numSubclusters);
				}
				return final_numSubclusters;
				//}
			} else {
				return newNumSubclusters;
			}
		}
	}
	
	return numSubclusters;
}

int32_t merge_subclusters_in_bin(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t niPET, uint32_t *pNumFuncCall) {
	int32_t numSubclusters = 0;
	int32_t i,j;
	
	*pNumFuncCall+=1;
	{
		// update the leader of the bin to encompass the anchor region
		i = 0;
		for(j=0; j<niPET; ++j) {
			cpu_cluster_record_t *clusterRecord = &(clusterRecords[j]);
			if (clusterRecord->iPET>0) {
				numSubclusters++;
				i = j;
			} else {
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
		
		if (bwa_verbose>=5) {
			if (numSubclusters>1) {
				kstring_t str;
				str.l = str.m = 0; str.s = 0;
				kputs("[D:SubCluster] pairId=", &str); kputul(clusterRecords[0].pairId, &str);
				kputs(" iPET=", &str); kputl(niPET, &str);
				kputs(" #subcluster=", &str); kputl(numSubclusters, &str);
				kputc('\n', &str);
				fprintf(stderr, "%s", str.s);
				free(str.s);
			}
		}
	}

	if (numSubclusters>1) {
		// prepare the summarized range
		i = 0;
		cpu_subcluster_record_t *subclusters = (cpu_subcluster_record_t *) calloc(numSubclusters, sizeof(cpu_subcluster_record_t));
		int32_t *subclustersLevel = (int32_t *) calloc(niPET, sizeof(int32_t));
		int32_t level = 0;
		if (subclusters && subclustersLevel) {
			for(j=0; j<niPET; ++j) {
				cpu_cluster_record_t *clusterRecord = &(clusterRecords[j]);
				if (clusterRecord->iPET>0) {
					cpu_subcluster_record_t *subcluster = &(subclusters[i]);
					subcluster->pos = clusterRecord->pos; subcluster->endpos = clusterRecord->endpos;
					subcluster->mpos = clusterRecord->mpos; subcluster->mendpos = clusterRecord->mendpos;
					subcluster->cid = clusterRecord->cid; subcluster->mcid = clusterRecord->mcid;
					subcluster->iPET = clusterRecord->iPET;
					subcluster->oid = clusterRecord->oid; subcluster->moid = clusterRecord->moid;
					subcluster->offset = j;
					subcluster->lineno = clusterRecord->lineno;
					subcluster->pairId = clusterRecord->pairId;
					subcluster->merged = -1;
					i++;
				}
			}
			
			// attempt to make new use of oid and moid
			qsort(subclusters, numSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterTag2);
			for(i=0; i<numSubclusters; ++i) { subclusters[i].moid = i; }
			qsort(subclusters, numSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterTag1);
			for(i=0; i<numSubclusters; ++i) { subclusters[i].oid = i; }
			
			if (bwa_verbose>=5) {
				fprintf(stderr, "[D:%s] Before merge_subclusters\n", __func__);
				report_Subclusters_List (stderr, subclusters, numSubclusters, subclustersLevel, NULL);
			}
			
			// attempt to merge the ranges
			int32_t final_numSubclusters = merge_subclusters (opt, subclusters, subclustersLevel, &level, numSubclusters, numSubclusters, pNumFuncCall);
			/*WCH*///final_numSubclusters = numSubclusters;
			
			// TODO: propagate the new merged ranges back to the cluster record!
			if (final_numSubclusters!=numSubclusters) {
				// TODO: update the subclusters level for proper update sequence
				for(i=0; i<numSubclusters; ++i) {
					subclusters[i].level = subclustersLevel[subclusters[i].offset];
				}
				
				if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] After merge_subclusters, #subclusters %d -> %d\n", __func__, numSubclusters, final_numSubclusters);
					report_Subclusters_List (stderr, subclusters, numSubclusters, subclustersLevel, NULL);
				}
				
				// TODO: need to flatten the update as a single level rather than a hierarchy order !!!
				qsort(subclusters, numSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterFlattenDependency);
				for(i=0; i<numSubclusters; ++i) {
					if (-1==subclusters[i].merged) {
						subclustersLevel[subclusters[i].offset] = subclusters[i].offset;
					} else {
						subclustersLevel[subclusters[i].offset] = subclustersLevel[subclusters[i].merged];
						subclusters[i].merged = subclustersLevel[subclusters[i].offset];
					}
				}
				
				if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] Merge_subclusters direct merge update, #subclusters %d -> %d\n", __func__, numSubclusters, final_numSubclusters);
					report_Subclusters_List (stderr, subclusters, numSubclusters, NULL, NULL);
				}
				
				if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] Before adjustment\n", __func__);
					report_iPETs_List (stderr, clusterRecords, niPET, NULL, -1);
				}
				
				// order largest to smaller so that we fold up rather than random
				//qsort(subclusters, numSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterUpdateOrder);
				qsort(subclusters, numSubclusters, sizeof(cpu_subcluster_record_t), cmpSubclusterUpdateOrderByLevel);
				if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] merge_subclusters update order\n", __func__);
					report_Subclusters_List (stderr, subclusters, numSubclusters, subclustersLevel, NULL);
				}
				
				// update the records with the new (cid, mcid) and the iPET counts
				int32_t numMergedSubclusters = 0;
				for(i=0; i<numSubclusters; ++i) {
					//if (-1!=subclusters[i].merged) {
					if (-1==subclusters[i].merged) {
						//these entries are not merged, retain
						break;
					}
					cpu_cluster_record_t *repCluster = &(clusterRecords[subclusters[i].merged]);
					cpu_cluster_record_t *toMergeCluster = &(clusterRecords[subclusters[i].offset]);
					
					// update the coordinates system
					if (repCluster->pos>toMergeCluster->pos) repCluster->pos = toMergeCluster->pos;
					if (repCluster->endpos<toMergeCluster->endpos) repCluster->endpos = toMergeCluster->endpos;
					if (repCluster->mpos>toMergeCluster->mpos) repCluster->mpos = toMergeCluster->mpos;
					if (repCluster->mendpos<toMergeCluster->mendpos) repCluster->mendpos = toMergeCluster->mendpos;
					
					// update the tags grouping
					cpu_cluster_record_t *cluster = &(clusterRecords[subclusters[i].offset]);
					for(j=0; j<toMergeCluster->iPET; ++j, ++cluster) {
						cluster->cid = repCluster->cid;
						cluster->mcid = repCluster->mcid;
					}
					
					// update the counts
					repCluster->iPET += toMergeCluster->iPET; toMergeCluster->iPET = 0;
					
					numMergedSubclusters++;
				}
				
				if (numMergedSubclusters>0) {
					//qsort(clusterRecords, niPET, sizeof(cpu_cluster_record_t), cmpPETBothTagsOrder);
					qsort(clusterRecords, niPET, sizeof(cpu_cluster_record_t), cmpPETBothTagsAfterMerger);
				}
				
				if (bwa_verbose>=5) {
					fprintf(stderr, "[D:%s] After adjustment\n", __func__);
					report_iPETs_List (stderr, clusterRecords, niPET, NULL, -1);
				}
				
			}
			
			free(subclusters);
			free(subclustersLevel);
		} else {
			if (subclusters) free(subclusters);
			else {
				fprintf(stderr, "[E::%s] out of memory for %d summary range e: %s\n", __func__, numSubclusters, strerror (errno));
			}
			if (subclustersLevel) free(subclustersLevel);
			else {
				fprintf(stderr, "[E::%s] out of memory for %d summary range levels e: %s\n", __func__, numSubclusters, strerror (errno));
			}
			//return -1;
		}
	}
	
	return numSubclusters;
}

int32_t cluster_iPETs(const cluster_opt_t *opt, cpu_cluster_record_t *clusterRecords, uint32_t nClusterKeys) {
//#define TRACE_CASE_DETAILS
//#define TRACE_BULK_CASES
	int32_t i;
	//bwa_verbose = 5;
	if (nClusterKeys>1) {
		for(i=0; i<nClusterKeys; ++i) {
			
#ifdef FUNC_CALL_COUNTING
			double t_started = realtime();
			clusterRecords[i].clusterCalls = 0;
			clusterRecords[i].mergeCalls = 0;
#endif
			// iPET will indicate the number of memebers in the bin
			// only those which needs further checking has .swap!=0
			if (clusterRecords[i].swap) {
				// if >=2 members, we will need to check if the overlap is proper
				int32_t niPET = clusterRecords[i].iPET;
				
#if 0
				if (2066!=niPET && 10664!=niPET && 19060!=niPET && 36807!=niPET && 22307!=niPET) continue;
#endif
				
#ifdef TRACE_CASE_DETAILS
				int bwa_verbose_orig = bwa_verbose;
				int32_t j;
				int nDump = 0;
				//if (915==niPET) {
					for(j=0; j<niPET; ++j) {
						if (84429932==clusterRecords[i+j].pairId) {
							nDump = 1;
							bwa_verbose = 5;
							break;
						}
					}
				//}
				
				if (0!=nDump) {
					fprintf(stderr, "[D:%s] DEBUGGING\n", __func__);
				}
				/*if (bwa_verbose>=3) {
					fprintf(stderr, "[D:%s] Processing i=%d pairId=%lu (%d iPETs)\n",
							__func__, i, clusterRecords[i].pairId, niPET);
				}*/
#endif
				if (niPET>=opt->bulkiPETsTrigger) {
#ifdef TRACE_BULK_CASES
					int bwa_verbose_orig = bwa_verbose;
					bwa_verbose = 5;
#endif
					int32_t nRetryDiscordantPETs = 0;
					long pairId = clusterRecords[i].pairId;
					int32_t j,k;
					
					uint32_t nClusterFuncCalls = 0;
					uint32_t nMergeFuncCalls = 0;
					
					if (bwa_verbose>=4) {
						fprintf(stderr, "[D:%s] Initiate chunk processing on large (%d iPETs) cluster %ld\n",
								__func__, niPET, pairId);
					}
					
					// update the system to ensure chunk processing
					int32_t cid=clusterRecords[i].cid; int32_t mcid=clusterRecords[i].mcid;
					int32_t bulkiPETsBlock = opt->bulkiPETsBlock;
					// TODO: heuristics to get the optimal run time
					//       need to get a better parameter function
					if (niPET>=30000) {
						bulkiPETsBlock = 500;
					} else if (niPET>=10000) {
						bulkiPETsBlock = 250;
					} else {
						bulkiPETsBlock = 100;
					}
					for(j=0; j<niPET; j+=bulkiPETsBlock) {
						int32_t nChunkSize = niPET - j;
						if (nChunkSize>bulkiPETsBlock) nChunkSize = bulkiPETsBlock;

						if (clusterRecords[i+j].cid!=cid || clusterRecords[i+j].mcid!=mcid) {
							fprintf(stderr, "[D:%s] Different cluster group! pairId=%ld cid=%d (expected %d), mcid=%d (expected %d)\n", __func__, clusterRecords[i+j].pairId, clusterRecords[i+j].cid, cid, clusterRecords[i+j].mcid, mcid);
						}
						
						for(k=0; k<nChunkSize; ++k) {
							clusterRecords[i+j+k].cid = clusterRecords[i+j].oid;
							clusterRecords[i+j+k].mcid = clusterRecords[i+j].moid;
						}
						clusterRecords[i+j].iPET = nChunkSize;
						clusterRecords[i+j].swap = 1;
					}
					
					// iterate thru' the new chunk set to get results
					int32_t nProcessedSize = 0;
					for(j=0; j<niPET; j+=bulkiPETsBlock) {
						int32_t nChunkSize = niPET - j;
						if (nChunkSize>bulkiPETsBlock) nChunkSize = bulkiPETsBlock;
						nProcessedSize += nChunkSize;
						
						if (bwa_verbose>=4) {
							fprintf(stderr, "[D:%s] Calling chunk (%d-%d,size=%d) clustering on %ld (%d iPETs)\n",
									__func__, j, j+nChunkSize, nChunkSize, pairId, niPET);
						}
						
						int32_t nChunkRetryDiscordantPETs = cluster_bin(opt, &(clusterRecords[i+j]), nChunkSize, 1, pairId, &nClusterFuncCalls);
						nRetryDiscordantPETs += nChunkRetryDiscordantPETs;
						
						if (bwa_verbose>=4) {
							fprintf(stderr, "[D:%s] Calling chunk (%d-%d,size=%d) merging on %ld (%d iPETs) with %d discordant iPETs\n",
									__func__, j, j+nChunkSize, nChunkSize, pairId, niPET, nChunkRetryDiscordantPETs);
						}
						
#if 1
						int32_t numChunkSubclusters = merge_subclusters_in_bin(opt, &(clusterRecords[i+j]), nChunkSize, &nMergeFuncCalls);
#else
						// we attempt to merge progressively, later test merging only some fraction?
						int32_t numChunkSubclusters = merge_subclusters_in_bin(opt, &(clusterRecords[i]), nProcessedSize);
#endif
						
						if (bwa_verbose>=4) {
							fprintf(stderr, "[D:%s] Chunk (%d-%d,size=%d) processed on %ld (%d iPETs) with %d discordant iPETs and %d subclusters\n",
									__func__, j, j+nChunkSize, nChunkSize, pairId, niPET, nChunkRetryDiscordantPETs, numChunkSubclusters);
						}
						
					}
					
					// DEBUG:
					if (bwa_verbose>=5) {
						fprintf(stderr, "[D:%s] Before merging all chunks from niPET==%d\n", __func__, niPET);
						report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
					}
					// DEBUG:END
					// perform the final merger
					int32_t numSubclusters = merge_subclusters_in_bin(opt, &(clusterRecords[i]), niPET, &nMergeFuncCalls);
					
					// DEBUG:
					if (bwa_verbose>=5) {
						fprintf(stderr, "[D:%s] After merging all chunks from niPET==%d\n", __func__, niPET);
						report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
					}
					// DEBUG:END

					if (bwa_verbose>=4) {
						fprintf(stderr, "[D:%s] Large (%d iPETs) cluster %ld processed with %d discordant iPETs and %d subclusters\n",
								__func__, niPET, pairId, nRetryDiscordantPETs, numSubclusters);
					}
					
					//if (bwa_verbose >= 5)
					{
						// for DEBUGGING; check that we have not lost the iPET count in each of the original bin group
						int32_t j;
						int32_t niPETinBin = 0;
						for(j=0; j<niPET; ++j) niPETinBin += clusterRecords[i+j].iPET;
						if (niPET!=niPETinBin) {
							fprintf(stderr, "[D:%s:iPET] i=%d, iPET=%d, iPETinBin=%d, #discordant=%d, #subclusters=%d\n", __func__, i, niPET, niPETinBin, nRetryDiscordantPETs, numSubclusters);
							fprintf(stderr, "[D:%s:iPET] i=%d, iPETs=", __func__, i);
							for(j=0; j<niPET; ++j)  fprintf(stderr, "%d,", clusterRecords[i+j].iPET);
							fprintf(stderr, "\n");
						}
					}
					
#ifdef FUNC_CALL_COUNTING
					clusterRecords[i].clusterCalls = nClusterFuncCalls;
					clusterRecords[i].mergeCalls = nMergeFuncCalls;
#endif
					
#ifdef TRACE_BULK_CASES
					bwa_verbose = bwa_verbose_orig;
#endif
					
				} else {
#ifdef TRACE_CASE_DETAILS
					// DEBUG:
					if (0!=nDump) {
						fprintf(stderr, "[D:%s] Chunk niPET==%d before cluster bin\n", __func__, niPET);
						report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
					}
					// DEBUG:END
#endif
					uint32_t nClusterFuncCalls = 0;
					uint32_t nMergeFuncCalls = 0;
					
					int32_t nRetryDiscordantPETs = cluster_bin(opt, &(clusterRecords[i]), niPET, 1, clusterRecords[i].pairId, &nClusterFuncCalls);
#ifdef TRACE_CASE_DETAILS
					// DEBUG:
					if (0!=nDump) {
						fprintf(stderr, "[D:%s] Chunk niPET==%d before cluster bin merger\n", __func__, niPET);
						report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
					}
					// DEBUG:END
#endif
					// TODO: we will now try to merge the subcluster(s) in this bin if possible
					int32_t numSubclusters = merge_subclusters_in_bin(opt, &(clusterRecords[i]), niPET, &nMergeFuncCalls);
					
#ifdef TRACE_CASE_DETAILS
					if (0!=nDump) {
						bwa_verbose = bwa_verbose_orig;
						fprintf(stderr, "[D:%s] Chunk niPET==%d after final merger\n", __func__, niPET);
						report_iPETs_List (stderr, &(clusterRecords[i]), niPET, NULL, -1);
					}
#endif
					//if (bwa_verbose >= 5)
					{
						// for DEBUGGING; check that we have not lost the iPET count in each of the original bin group
						int32_t j;
						int32_t niPETinBin = 0;
						for(j=0; j<niPET; ++j) niPETinBin += clusterRecords[i+j].iPET;
						if (niPET!=niPETinBin) {
							fprintf(stderr, "[D:%s:iPET] i=%d, iPET=%d, iPETinBin=%d, #discordant=%d, #subclusters=%d\n", __func__, i, niPET, niPETinBin, nRetryDiscordantPETs, numSubclusters);
							fprintf(stderr, "[D:%s:iPET] i=%d, iPETs=", __func__, i);
							for(j=0; j<niPET; ++j)  fprintf(stderr, "%d,", clusterRecords[i+j].iPET);
							fprintf(stderr, "\n");
						}
					}
					
#ifdef FUNC_CALL_COUNTING
					clusterRecords[i].clusterCalls = nClusterFuncCalls;
					clusterRecords[i].mergeCalls = nMergeFuncCalls;
#endif
					
				}
#ifdef TRACE_CASE_DETAILS
				bwa_verbose = bwa_verbose_orig;
#endif
				
#ifdef FUNC_CALL_COUNTING
				double t_diffProcessing = realtime() - t_started;
				if (bwa_verbose>=4) {
					fprintf(stderr, "[D:%s] Function calls; pairid %lu iPETs %d clusterfn %u mergefn %u duration %.2f sec\n", __func__, clusterRecords[i].pairId, niPET, clusterRecords[i].clusterCalls, clusterRecords[i].mergeCalls, t_diffProcessing);
				}
#endif
				i += (niPET-1);
			}
		}
	}
	
	return 0;
}

void firstPassClusterParition (const cluster_opt_t *opt, cpu_cluster_record_t *clusterKeys, uint32_t nClusterKeys, uint32_t *pNumUsables, uint32_t *pNumClusterBins)
{
	uint32_t i;
	double t_diffProcessing;
	double t_timepointIO, t_timepointProcessing;
	
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
	
	// Let's get the number of usable rows
	*pNumUsables = nClusterKeys;
	for(i=0; i<nClusterKeys; ++i) {
		if (0==clusterKeys[i].usable) {
			*pNumUsables = i;
			break;
		}
	}
	
	// let's sort based on tag left and process them
	t_timepointProcessing = realtime();
	//qsort(clusterKeys, nClusterKeys, sizeof(cpu_cluster_record_t), cmpPETTag1);
	qsort(clusterKeys, *pNumUsables, sizeof(cpu_cluster_record_t), cmpPETTag1);
	t_timepointIO = realtime();
	t_diffProcessing = t_timepointIO - t_timepointProcessing;
	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Sorted %u left key records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
	}
	
	// let's process sorted left tags
	t_timepointProcessing = realtime();
	//cluster_PETags1(opt, clusterKeys, nClusterKeys);
	cluster_PETags1(opt, clusterKeys, *pNumUsables);
	t_timepointIO = realtime();
	t_diffProcessing = t_timepointIO - t_timepointProcessing;
	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Clustered %u left key records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
	}
	
	// let's sort based on combine tag
	t_timepointProcessing = realtime();
	qsort(clusterKeys, *pNumUsables, sizeof(cpu_cluster_record_t), cmpPETBothTags);
	t_timepointIO = realtime();
	t_diffProcessing = t_timepointIO - t_timepointProcessing;
	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Sorted %u combined key records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
	}
	
	// let's process them
	t_timepointProcessing = realtime();
	*pNumClusterBins = 0;
	setup_cluster_bin(opt, clusterKeys, *pNumUsables, pNumClusterBins);
	t_timepointIO = realtime();
	t_diffProcessing = t_timepointIO - t_timepointProcessing;
	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Clustered %u combined key records, %u cluster bins, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, *pNumClusterBins, t_diffProcessing, t_diffProcessing/60.0);
	}
}

int cmpPETForJuicerAll(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int32_t i32comp; //= ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	// sort by chrom1, chrom2
	//int32_t
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	
	//new
	i32comp = (ib->mendpos-ib->pos) - (ia->mendpos-ia->pos); if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	
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

int cmpPETForJuicerDampenSelfLigation(const void *a, const void *b) {
	const cpu_cluster_record_t *ia = (const cpu_cluster_record_t *)a;
	const cpu_cluster_record_t *ib = (const cpu_cluster_record_t *)b;
	
	int32_t i32comp = ib->usable - ia->usable; if (0!=i32comp) return i32comp;
	// sort by chrom1, chrom2
	//int32_t
	i32comp = ia->tid - ib->tid; if (0!=i32comp) return i32comp;
	i32comp = ia->mtid - ib->mtid; if (0!=i32comp) return i32comp;
	
	i32comp = ia->cid - ib->cid; if (0!=i32comp) return i32comp;
	i32comp = ia->mcid - ib->mcid; if (0!=i32comp) return i32comp;
	
	//new
	i32comp = (ib->mendpos-ib->pos) - (ia->mendpos-ia->pos); if (0!=i32comp) return i32comp;
	
	i32comp = ia->pos - ib->pos; if (0!=i32comp) return i32comp;
	i32comp = ia->mendpos - ib->mendpos; if (0!=i32comp) return i32comp;
	
	i32comp = ib->endpos - ia->endpos; if (0!=i32comp) return i32comp;
	i32comp = ia->mpos - ib->mpos; if (0!=i32comp) return i32comp;
	
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

void orderForJuicer (cpu_cluster_record_t *clusterKeys, uint32_t nClusterKeys, uint32_t excludeSelfLigation)
{
	double t_diffProcessing;
	double t_timepointIO, t_timepointProcessing;
	
	// let's sort based on combine tag
	t_timepointProcessing = realtime();
	if (0!=excludeSelfLigation) {
		qsort(clusterKeys, nClusterKeys, sizeof(cpu_cluster_record_t), cmpPETForJuicerDampenSelfLigation);
	} else {
		qsort(clusterKeys, nClusterKeys, sizeof(cpu_cluster_record_t), cmpPETForJuicerAll);
	}
	t_timepointIO = realtime();
	t_diffProcessing = t_timepointIO - t_timepointProcessing;
	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Sorted %u iPETs for juicer, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
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
	
	//void (*extension_worker)(void*,int,int) = tags_midpoint_then_extension_worker;
	void (*extension_worker)(void*,int,int) = tags_extension_then_midpoint_worker;
	int (*computeIntraSpan)(int32_t,int32_t,int32_t,int32_t) = getIntraSpan;

	
	FILE *outfd = NULL;
	kstring_t clusterKeyFilename = {0,0,0};

	gzFile clustersfd = NULL;
	kstring_t clustersFilename = {0,0,0};
	
	gzFile clustersTransfd = NULL;
	kstring_t clustersTransFilename = {0,0,0};
	
	gzFile juicerfd = NULL;
	kstring_t juicerFilename = {0,0,0};
	
	init_CPUBamBuffer(&bamsBuffer);
	
	opt = cluster_opt_init();
	strcpy(in_mode, "r");
	while ((c = getopt(argc, argv, "mMdSjxgf:s:t:O:l:u:b:5:3:i:v:B:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'v') { bwa_verbose = atoi(optarg)<0 ? 3 : atoi(optarg); }
		else if (c == 'l') opt->outputSortedList = atoi(optarg), opt->outputSortedList = opt->outputSortedList > 0 ? opt->outputSortedList : 0;
		else if (c == 'B') opt->iterationProgressBlock = atoi(optarg), opt->iterationProgressBlock = opt->iterationProgressBlock > 0 ? opt->iterationProgressBlock : SUBCLUSTER_ITERACTION_BLOCK;
		else if (c == 'b') {
			opt->bulkiPETsTrigger = (int32_t) strtol(optarg, &p, 10);
			if (opt->bulkiPETsTrigger<=0) opt->bulkiPETsTrigger = BULK_IPET_TRIGGER;
			if (*p != 0 && ispunct(*p)) {
				opt->bulkiPETsBlock = (int32_t) strtol(p+1, &p, 10);
				if (opt->bulkiPETsBlock<=0) opt->bulkiPETsBlock = BULK_IPET_BLOCK;
			}
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] Bulk cluster trigger with %d iPETs in chunk of %d iEPTs\n",
						__func__, opt->bulkiPETsTrigger, opt->bulkiPETsBlock);
		}
		else if (c == 'i') opt->cycle = atoi(optarg);
		else if (c == 'd') opt->disponoParallel = 1;
		else if (c == 'S') is_bamin = 0;
		else if (c == 'M') {
			if (1==opt->midPoint) {
				fprintf(stderr, "[M::%s] Options -M and -m are mutually exclusive\n",__func__);
				free(opt);
				return 1;
			}
			opt->midPoint = -1;
			extension_worker = tags_midpoint_then_extension_worker;
			computeIntraSpan = getIntraSpan_MidPoint;

			if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] Using mid-point of fragment as reference for extension and span computation\n",__func__);
		}
		else if (c == 'm') {
			if (-1==opt->midPoint) {
				fprintf(stderr, "[M::%s] Options -m and -M are mutually exclusive\n",__func__);
				free(opt);
				return 1;
			}
			opt->midPoint = 1;
			extension_worker = tags_extension_then_midpoint_worker;
			computeIntraSpan = getIntraSpan_MidPoint;
			if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] Using mid-point of extended fragment for span computation\n",__func__);
		}
		else if (c == 'O') {
			ks_set(optarg, &(opt->outputPrefix));
		}
		else if (c == 's') opt->selfLigation = atoi(optarg), opt->selfLigation = opt->selfLigation > 0 ? opt->selfLigation : G_SELF_LIGATION;
		else if (c == '5') {
			opt->extension5p.source = (int32_t) strtol(optarg, &p, 10);
			if (5!=opt->extension5p.source && 3!=opt->extension5p.source)
				fprintf(stderr, "[E::%s] 5' extension source must be either 5 or 3\n", __func__);
			if (*p != 0 && ispunct(*p) && (isdigit(p[1])||'+'==p[1]||'-'==p[1]))
				opt->extension5p.size = (int32_t) strtol(p+1, &p, 10);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] 5' extension by %d bp from %d'\n",
						__func__, opt->extension5p.size, opt->extension5p.source);
		}
		else if (c == '3') {
			opt->extension3p.source = (int32_t) strtol(optarg, &p, 10);
			if (5!=opt->extension3p.source && 3!=opt->extension3p.source)
				fprintf(stderr, "[E::%s] 3' extension source must be either 5 or 3\n", __func__);
			if (*p != 0 && ispunct(*p) && (isdigit(p[1])||'+'==p[1]||'-'==p[1]))
				opt->extension3p.size = (int32_t) strtol(p+1, &p, 10);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] 3' extension by %d bp from %d'\n",
						__func__, opt->extension3p.size, opt->extension3p.source);
		} else if (c == 'j') {
			opt->juicerOutput |= CPU_CLUSTER_JUICER;
		} else if (c == 'x') {
			opt->juicerOutput &= ~CPU_CLUSTER_JUICER_SELF_LIGATION;
		} else if (c == 'f') {
			if (0==strcmp(optarg, CPU_CLUSTER_FORMAT_CHIASIG_STR)) {
				opt->outputFormat &= ~CPU_CLUSTER_FORMAT_MASK;
				opt->outputFormat |= CPU_CLUSTER_FORMAT_CHIASIG;
			} else if (0==strcmp(optarg, CPU_CLUSTER_FORMAT_BEDPE_STR)) {
				opt->outputFormat &= ~CPU_CLUSTER_FORMAT_MASK;
				opt->outputFormat |= CPU_CLUSTER_FORMAT_BEDPE;
			} else {
				fprintf(stderr, "[E::%s] Unrecognized output format %s\n", __func__, optarg);
			}
		} else if (c == 'g') {
			opt->outputFormat |= CPU_CLUSTER_FORMAT_SEPARATE;
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
		fprintf(stderr, "       -m         use mid-point after fragment extension\n");
		fprintf(stderr, "       -M         use mid-point before fragment extension\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -S         input sam format\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -O STR     output prefix [%s]\n", opt->outputPrefix.s);
		fprintf(stderr, "       -j         output all reads in juicer short format\n");
		fprintf(stderr, "       -x         exclude self-ligation from juicer output (used along with -j)\n");
		fprintf(stderr, "       -f STR     output format: [%s]\n", CPU_CLUSTER_FORMAT_CHIASIG_STR);
		fprintf(stderr, "                    %s - bedpe with Col 7=count\n", CPU_CLUSTER_FORMAT_CHIASIG_STR);
		fprintf(stderr, "                    %s - as defined by bedtools; Col 7=name, Col 8=score\n", CPU_CLUSTER_FORMAT_BEDPE_STR);
		fprintf(stderr, "       -g         output intra-and-inter-chromosomal interaction separately\n");
		fprintf(stderr, "       -l         output sorted PETs to stdout\n");
		fprintf(stderr, "       -i INT     output iteration#i's sorted PETs to stdout\n");
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
		
		if (CPU_CLUSTER_FORMAT_SEPARATE==(opt->outputFormat & CPU_CLUSTER_FORMAT_SEPARATE)) {
			if (CPU_CLUSTER_FORMAT_CHIASIG == (opt->outputFormat & CPU_CLUSTER_FORMAT_MASK))  {
				ksprintf(&clustersFilename, "%s.clusters.cis.chiasig.gz", opt->outputPrefix.s);
			} else if (CPU_CLUSTER_FORMAT_BEDPE == (opt->outputFormat & CPU_CLUSTER_FORMAT_MASK))  {
				ksprintf(&clustersFilename, "%s.clusters.cis.bedpe.gz", opt->outputPrefix.s);
			}
			const char *mode = "w6";
			clustersfd = gzopen (clustersFilename.s, mode);
			if (clustersfd == NULL) {
				fprintf(stderr, "[E::%s] fail to open '%s' for writing: %s\n", __func__, clustersFilename.s, strerror (errno));
				cluster_opt_terminate(opt);
				free(opt);
				fclose(outfd);
				return 1;
			}
			gzbuffer(clustersfd, G_GZIP_BUFFER_SIZE);

			if (CPU_CLUSTER_FORMAT_CHIASIG == (opt->outputFormat & CPU_CLUSTER_FORMAT_MASK))  {
				ksprintf(&clustersTransFilename, "%s.clusters.trans.chiasig.gz", opt->outputPrefix.s);
			} else if (CPU_CLUSTER_FORMAT_BEDPE == (opt->outputFormat & CPU_CLUSTER_FORMAT_MASK))  {
				ksprintf(&clustersTransFilename, "%s.clusters.trans.bedpe.gz", opt->outputPrefix.s);
			}
			clustersTransfd = gzopen (clustersTransFilename.s, mode);
			if (clustersTransfd == NULL) {
				fprintf(stderr, "[E::%s] fail to open '%s' for writing: %s\n", __func__, clustersTransFilename.s, strerror (errno));
				cluster_opt_terminate(opt);
				free(opt);
				gzflush(clustersfd, Z_FINISH); // there was an error.. no point tracking these
				gzclose(clustersfd); // there was an error.. no point tracking these
				fclose(outfd);
				return 1;
			}
			gzbuffer(clustersTransfd, G_GZIP_BUFFER_SIZE);

		} else {
			if (CPU_CLUSTER_FORMAT_CHIASIG == (opt->outputFormat & CPU_CLUSTER_FORMAT_MASK))  {
				ksprintf(&clustersFilename, "%s.clusters.chiasig.gz", opt->outputPrefix.s);
			} else if (CPU_CLUSTER_FORMAT_BEDPE == (opt->outputFormat & CPU_CLUSTER_FORMAT_MASK))  {
				ksprintf(&clustersFilename, "%s.clusters.bedpe.gz", opt->outputPrefix.s);
			}
			const char *mode = "w6";
			clustersfd = gzopen (clustersFilename.s, mode);
			if (clustersfd == NULL) {
				fprintf(stderr, "[E::%s] fail to open '%s' for writing: %s\n", __func__, clustersFilename.s, strerror (errno));
				cluster_opt_terminate(opt);
				free(opt);
				fclose(outfd);
				return 1;
			}
			gzbuffer(clustersfd, G_GZIP_BUFFER_SIZE);
		}
		
		if (CPU_CLUSTER_JUICER==(opt->juicerOutput & CPU_CLUSTER_JUICER)) {
			// TODO: write as gz file!
			ksprintf(&juicerFilename, "%s.juice.gz", opt->outputPrefix.s);
			const char *mode = "w6";
			juicerfd = gzopen (juicerFilename.s, mode);
			if (juicerfd == NULL) {
				fprintf(stderr, "[E::%s] fail to open '%s' for writing: %s\n", __func__, juicerFilename.s, strerror (errno));
				cluster_opt_terminate(opt);
				free(opt);
				gzflush(clustersfd, Z_FINISH); // there was an error.. no point tracking these
				gzclose(clustersfd); // there was an error.. no point tracking these
				if (CPU_CLUSTER_FORMAT_SEPARATE==(opt->outputFormat & CPU_CLUSTER_FORMAT_SEPARATE)) {
					gzflush(clustersTransfd, Z_FINISH); // there was an error.. no point tracking these
					gzclose(clustersTransfd); // there was an error.. no point tracking these
				}
				fclose(outfd);
				return 1;
			}
			gzbuffer(juicerfd, G_GZIP_BUFFER_SIZE);
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
		
		//int nNumRequested = opt->chunk_size * opt->n_threads;
		int nNumRequested = opt->chunk_size;
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
			cluster_process_alns(opt, n_processed, n, bamRecs, readids, &numTagGroups, tagGroups, clusterRecords, &bamsBuffer, 0, computeIntraSpan);
			t_timepointProcessing = realtime();
			t_diffProcessing = t_timepointProcessing - t_timepointIO;
			
			// TODO: future: write the bam results; for testing purposes
			// TODO: future: need to make sure that we have either 1 record (SE) or 2 records (PE)
			
			if (CPU_DEDUP_DEBUG_CONVERT_LIST==(CPU_DEDUP_DEBUG_CONVERT_LIST&opt->outputSortedList)) {
				report_iPETs_List (stdout, clusterRecords, numTagGroups, in->header, -1);
			}
			
			// write the binary data file
			err_fwrite(&numTagGroups, sizeof(numTagGroups), 1, outfd); // write the number of record filler for total number of records
			uint32_t nSorted = 0; err_fwrite(&nSorted, sizeof(nSorted), 1, outfd); // sort state: unsorted=0, chunk sorted=1, fully sorted=2
			err_fwrite(clusterRecords, sizeof(cpu_cluster_record_t), numTagGroups, outfd);
			
			t_timepointIO = realtime();
			t_diffIO += (t_timepointIO - t_timepointProcessing);
			// clean up
			for(i=0; i<n; ++i) free(bamRecs[i].data);
			free(bamRecs);
			free(readids);
			free(tagGroups);
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
		cpu_cluster_record_t *clusterKeys = NULL;
		clusterKeys = readClusterKeys(clusterKeyFilename.s, &nClusterKeys);
		t_timepointIO = realtime();
		t_diffProcessing = t_timepointIO - t_timepointProcessing;
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[M::%s] read %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
		}
		
		if (bwa_verbose>=4) {
			// let's keep the binary file for troubleshooting
		} else {
			// standard clean up!
			if (0!=unlink(clusterKeyFilename.s)) {
				fprintf(stderr, "[E::%s] fail to remove temporary file '%s': %s\n", __func__, clusterKeyFilename.s, strerror (errno));
			}
		}
		
		if (!clusterKeys) {
			fprintf(stderr, "[E::%s] Fail to read %u key records, terminating..\n", __func__, nClusterKeys);
		} else {
			if (CPU_DEDUP_DEBUG_READ_LIST==(CPU_DEDUP_DEBUG_READ_LIST&opt->outputSortedList)) {
				// for debugging purposes only
				t_timepointProcessing = realtime();
				fprintf(stdout, "\n\n===\tREAD LIST\n");
				report_iPETs_List (stdout, clusterKeys, nClusterKeys, in->header, -1);
				fprintf(stdout, "===\tEND : READ LIST\n");
				fflush(stdout);
				t_timepointIO = realtime();
				t_diffProcessing = t_timepointIO - t_timepointProcessing;
				if (bwa_verbose >= 3) {
					fprintf(stderr, "[M::%s] Writing read %u key records, i/o %.2f sec (%.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
				}
				// END - for debugging purposes only, remove from production code
			}
			
			// write sorted reads for juicer
			if (CPU_CLUSTER_JUICER==(opt->juicerOutput & CPU_CLUSTER_JUICER)) {
				uint32_t excludeSelfLigation = (CPU_CLUSTER_JUICER_SELF_LIGATION==(opt->juicerOutput & CPU_CLUSTER_JUICER_SELF_LIGATION)) ? 0 : 1;
				orderForJuicer (clusterKeys, nClusterKeys, excludeSelfLigation);
				report_Reads_Juicer(&juicerfd, clusterKeys, nClusterKeys, in->header, excludeSelfLigation);
				
				// let's close the juicer file asap for user to start using
				if (gzflush(juicerfd, Z_FINISH)) {
					fprintf(stderr, "[E::%s] fail to flush juicer data for gzip stream: %s\n", __func__, strerror (errno));
					ret = 1;
				}
				if (gzclose(juicerfd)) {
					fprintf(stderr, "[E::%s] fail to close juicer data stream: %s\n", __func__, strerror (errno));
				}
			}
			
			// let's extend the tags accordingly
			t_timepointProcessing = realtime();
			extend_PETs(opt, clusterKeys, nClusterKeys, in->header->target_len, in->header->n_targets, extension_worker);
			t_timepointIO = realtime();
			t_diffProcessing = t_timepointIO - t_timepointProcessing;
			if (bwa_verbose >= 3) {
				fprintf(stderr, "[M::%s] Extend %u PET records, processing %.2f sec ( %.2f min)..\n", __func__, nClusterKeys, t_diffProcessing, t_diffProcessing/60.0);
			}

			// let's partition the PETs into clusters that do not overlap
			uint32_t numUsables = 0;
			uint32_t numClusterBins = 0;
			firstPassClusterParition (opt, clusterKeys, nClusterKeys, &numUsables, &numClusterBins);
			
			// for DEBUGGING
			report_Total_iPET(opt, clusterKeys, nClusterKeys);
			
			report_iPETs_Distribution(opt, clusterKeys, numUsables, numClusterBins);
			
			if (CPU_DEDUP_DEBUG_SORTED_LIST==(CPU_DEDUP_DEBUG_SORTED_LIST&opt->outputSortedList) && (0==opt->cycle)) {
				// for debugging purposes only, remove from production code
				t_timepointProcessing = realtime();
				fprintf(stdout, "\n\n===\tSORTED LIST - Check Cycle #%d\n", 0);
				report_iPETs_List (stdout, clusterKeys, nClusterKeys, in->header, -1);
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
				nDiscordantPETs = cluster_iPETs(opt, clusterKeys, numUsables);
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
					report_iPETs_List (stdout, clusterKeys, nClusterKeys, in->header, -1);
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
				report_iPETs_List (stdout, clusterKeys, nClusterKeys, in->header, -1);
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
			if (CPU_CLUSTER_FORMAT_CHIASIG == (opt->outputFormat & CPU_CLUSTER_FORMAT_MASK))  {
				report_Clusters_List_CHIASIG(&clustersfd, &clustersTransfd, clusterKeys, numUsables, in->header);
			} else if (CPU_CLUSTER_FORMAT_BEDPE == (opt->outputFormat & CPU_CLUSTER_FORMAT_MASK))  {
				report_Clusters_List_BEDPE(&clustersfd, &clustersTransfd, clusterKeys, numUsables, in->header);
			}
			
			// release the sorted records as we no longer need them
			free(clusterKeys);
		}
		
		//-------------------------------------------------------------------------------------------------------------------
		
		samclose(in);
		
	}
	// END - PROCESSING
	
	free(fn_list);
	
	// close files, free and return
	// TODO: decide if we need additional closing due to re-opening
	//samclose(in);
	
	terminate_CPUBamBuffer(&bamsBuffer);
	
	//TODO: fclose(outfd);
	free(clusterKeyFilename.s);
	
	if (gzflush(clustersfd, Z_FINISH)) {
		fprintf(stderr, "[E::%s] fail to flush cluster data for gzip stream: %s\n", __func__, strerror (errno));
		ret = 1;
	}
	if (gzclose(clustersfd)) {
		fprintf(stderr, "[E::%s] fail to close cluster data stream: %s\n", __func__, strerror (errno));
	}
	if (CPU_CLUSTER_FORMAT_SEPARATE==(opt->outputFormat & CPU_CLUSTER_FORMAT_SEPARATE)) {
		if (gzflush(clustersTransfd, Z_FINISH)) {
			fprintf(stderr, "[E::%s] fail to flush cluster data for gzip stream: %s\n", __func__, strerror (errno));
			ret = 1;
		}
		if (gzclose(clustersTransfd)) {
			fprintf(stderr, "[E::%s] fail to close cluster data stream: %s\n", __func__, strerror (errno));
		}
		free(clustersTransFilename.s);
	}
	free(clustersFilename.s);

	// juicer file is closed the moment it is ready for use
	free(juicerFilename.s);
	
	cluster_opt_terminate(opt);
	free(opt);
	
	return ret;
}
