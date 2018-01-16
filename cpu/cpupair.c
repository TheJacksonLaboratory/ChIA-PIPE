// TODO:
// 1. dynamics allocation of read length
// 2. inject PG with PP chain

#include "zlib.h"
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
//WCH
//extern unsigned char nst_nt4_table[256];

#include "sam_header.h"
#include "sam.h"
#include "cputagcommon.h"
#include "cputagbam.h"
#include "cpuspan.h"

extern double G_t_real;
extern const char *CPU_version;

#define PAIR_F_REPEAT     0x01
#define PAIR_F_NOT_MAPEED 0x02

#define CPU_PAIR_OUTPUT_PREFIX "cppairs"

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
#define CPU_PAIR_OUTPUT_UU_FAILED     5
#define CPU_PAIR_OUTPUT_COUNT         (CPU_PAIR_OUTPUT_UU_FAILED+1)

#define G_MAX_QNAME_LEN 256

#define NUM_PASS_FAIL_STATES	(PF_STATES_FAIL+1)

typedef struct {
	long n_discordant;
	long n_differentChrom[NUM_PASS_FAIL_STATES];
	long n_sameChrom[NUM_PASS_FAIL_STATES];
	long n_sameChromLow[NUM_PASS_FAIL_STATES];
	long n_sameChromHigh[NUM_PASS_FAIL_STATES];
	long counts[CPU_MAPSTATE_COUNT][CPU_MAPSTATE_COUNT];
} pair_stat_t;


typedef struct {
	int chunk_size;
	int n_threads;
	int flag;               // see MEM_F_* macros
	
	int pairFlag;           // see PAIR_F_*macros
	int selfLigation;
	
	kstring_t outputPrefix;

	pairqc_opt_t pairqc_opt;
} pair_opt_t;

pair_opt_t *pair_opt_init()
{
	pair_opt_t *o;
	o = calloc(1, sizeof(pair_opt_t));
	
	o->chunk_size = 1000000;
	o->n_threads = 1;
	o->flag = 0;
	
	o->pairFlag = 0;
	o->selfLigation = G_SELF_LIGATION;
	
	o->pairqc_opt.minFractionalMAPQ = 0.95;
	//o->pairqc_opt.minMAPQ = 19;
	o->pairqc_opt.minMAPQ = -1;
	o->pairqc_opt.maxX1 = 10; // equavalent to mapq 19
	
	memset(&(o->outputPrefix), 0, sizeof(kstring_t));
	return o;
}

void pair_opt_terminate(pair_opt_t *o)
{
	free(o->outputPrefix.s);
}

static inline int init_CPPairs_Ouputs (const char *outPrefix, char *out_mode, bam_header_t *header, samfile_t **outfds) {
	// CPU_PAIR_OUTPUT_UU
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.UU.bam", outPrefix);
		outfds[CPU_PAIR_OUTPUT_UU] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_PAIR_OUTPUT_UU]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_PAIR_OUTPUT_UU+1;
		}
		free(filename.s);
	}
	// CPU_PAIR_OUTPUT_UxxU
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.UxxU.bam", outPrefix);
		outfds[CPU_PAIR_OUTPUT_UxxU] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_PAIR_OUTPUT_UxxU]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_PAIR_OUTPUT_UxxU+1;
		}
		free(filename.s);
	}
	// CPU_PAIR_OUTPUT_xx
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.xx.bam", outPrefix);
		outfds[CPU_PAIR_OUTPUT_xx] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_PAIR_OUTPUT_xx]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_PAIR_OUTPUT_xx+1;
		}
		free(filename.s);
	}
	// CPU_PAIR_OUTPUT_nn
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.nn.bam", outPrefix);
		outfds[CPU_PAIR_OUTPUT_nn] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_PAIR_OUTPUT_nn]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_PAIR_OUTPUT_nn+1;
		}
		free(filename.s);
	}
	
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.UU.discordant.bam", outPrefix);
		outfds[CPU_PAIR_OUTPUT_UU_DISCORDANT] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_PAIR_OUTPUT_UU_DISCORDANT]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_PAIR_OUTPUT_UU_DISCORDANT+1;
		}
		free(filename.s);
	}
	
	{
		kstring_t filename; memset(&filename, 0, sizeof(kstring_t));
		ksprintf(&filename, "%s.UU.failed.bam", outPrefix);
		outfds[CPU_PAIR_OUTPUT_UU_FAILED] = samopen(filename.s, out_mode, header);
		if (!outfds[CPU_PAIR_OUTPUT_UU_FAILED]) {
			fprintf(stderr, "[E::%s] fail to open '%s' for writing\n", __func__, filename.s);
			free(filename.s);
			return CPU_PAIR_OUTPUT_UU_FAILED+1;
		}
		free(filename.s);
	}
	
	return 0;
}

static inline int terminate_CPPairs_Ouputs (samfile_t **outfds) {
	// TODO: always synchronized with #defined
	int nFailed = 0;
	if (outfds[CPU_PAIR_OUTPUT_UU]) {
		samclose(outfds[CPU_PAIR_OUTPUT_UU]);
	}
	if (outfds[CPU_PAIR_OUTPUT_UxxU]) {
		samclose(outfds[CPU_PAIR_OUTPUT_UxxU]);
	}
	if (outfds[CPU_PAIR_OUTPUT_xx]) {
		samclose(outfds[CPU_PAIR_OUTPUT_xx]);
	}
	if (outfds[CPU_PAIR_OUTPUT_nn]) {
		samclose(outfds[CPU_PAIR_OUTPUT_nn]);
	}
	if (outfds[CPU_PAIR_OUTPUT_UU_DISCORDANT]) {
		samclose(outfds[CPU_PAIR_OUTPUT_UU_DISCORDANT]);
	}
	if (outfds[CPU_PAIR_OUTPUT_UU_FAILED]) {
		samclose(outfds[CPU_PAIR_OUTPUT_UU_FAILED]);
	}
	
	return nFailed;
}

typedef struct {
	const pair_opt_t *opt;
	
	int64_t n_processed;
	
	bam1_t *bamRecs;
	cpu_readid_t *readids;
	
	int nTagGroups;
	tag_group_t *tagGroups;
	
	pair_stat_t *pair_stats;
	
	int (*pairQCCheck)(const bam1_t *,const bam1_t *,const pairqc_opt_t *);
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

static void pair_worker(void *data, int i, int tid)
{
	int j, k;
	worker_t *w = (worker_t*)data;
	// pairing of tags from read(s) to form interaction does not care if reads are from paired-end
	
	// perform pairing
	w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_nn;
	
	// STAGE 1: we set up the representative tags in reads for pairing selection
	// TODO: we ignore sceondary alignment for now
	int32_t readTags[2][2] = {{-1,-1},{-1,-1}};
	uint8_t readTagsMapState[2][2] = {{CPU_MAPSTATE_NA,CPU_MAPSTATE_NA},{CPU_MAPSTATE_NA,CPU_MAPSTATE_NA}};
	uint8_t numTagsSet = 0;
	for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
		if (BAM_FSECONDARY==(w->bamRecs[k].core.flag & BAM_FSECONDARY)) {
			w->tagGroups[i].secondaryTag = 1;
		} else {
			if (-1==readTags[w->readids[k].readId][w->readids[k].tagId]) {
				readTags[w->readids[k].readId][w->readids[k].tagId] = k;
				bam1_t *bam = &(w->bamRecs[k]);
				if (BAM_FUNMAP==(bam->core.flag & BAM_FUNMAP)) {
					readTagsMapState[w->readids[k].readId][w->readids[k].tagId] = CPU_MAPSTATE_NOT_MAPPED;
				} else {
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
	uint8_t anchor1MapState = CPU_MAPSTATE_NA;
	int anchor1ReadIndex = -1;
	if (CPU_MAPSTATE_UNIQUE==readTagsMapState[0][0]) {
		if (CPU_MAPSTATE_UNIQUE==readTagsMapState[1][1]) {
			// both unique
			if (w->bamRecs[readTags[0][0]].core.tid == w->bamRecs[readTags[1][1]].core.tid) {
				// TODO: we pick the 5'-most read for now as it should has more reliable bases
				//       in future we might consider the length, etc
				anchor1ReadIndex = readTags[0][0];
			} else {
				// different chromosome!
			}
		} else {
			// only R/1 tag 1 unique
			anchor1ReadIndex = readTags[0][0];
		}
		anchor1MapState = CPU_MAPSTATE_UNIQUE;
	} else if (CPU_MAPSTATE_UNIQUE==readTagsMapState[1][1]) {
		// only R/2 tag 2 unique
		anchor1ReadIndex = readTags[1][1];
		anchor1MapState = CPU_MAPSTATE_UNIQUE;
	} else {
		// other cases
		if (CPU_MAPSTATE_REPEAT==readTagsMapState[0][0] || CPU_MAPSTATE_REPEAT==readTagsMapState[1][1]) {
			// either repeat, considered as repeat
			anchor1MapState = CPU_MAPSTATE_REPEAT;
		} else {
			// not mapped, too many Ns, and not available
			anchor1MapState = CPU_MAPSTATE_NOT_MAPPED;
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
				anchor2ReadIndex = readTags[1][0];
			} else {
				// different chromosome!
			}
		} else {
			// only R/2 tag 1 unique
			anchor2ReadIndex = readTags[1][0];
		}
		anchor2MapState = CPU_MAPSTATE_UNIQUE;
	} else if (CPU_MAPSTATE_UNIQUE==readTagsMapState[0][1]) {
		// only R/1 tag 2 unique
		anchor2ReadIndex = readTags[0][1];
		anchor2MapState = CPU_MAPSTATE_UNIQUE;
	} else {
		// other cases
		if (CPU_MAPSTATE_REPEAT==readTagsMapState[1][0] || CPU_MAPSTATE_REPEAT==readTagsMapState[0][1]) {
			// either repeat, considered as repeat
			anchor2MapState = CPU_MAPSTATE_REPEAT;
		} else {
			// not mapped, too many Ns, and not available
			anchor2MapState = CPU_MAPSTATE_NOT_MAPPED;
		}
	}
	
	// STATISTICS: update the count
	w->pair_stats[tid].counts[anchor1MapState][anchor2MapState]++;
	
	// STATE 3: record interaction-pair
	if (CPU_MAPSTATE_UNIQUE==anchor1MapState && CPU_MAPSTATE_UNIQUE==anchor2MapState) {
		// [UU] both unique
		if (-1!=anchor1ReadIndex && -1!=anchor2ReadIndex) {
			// let's copy most of the structure
			w->tagGroups[i].bamRecs[0] = bam_dup1(&(w->bamRecs[anchor1ReadIndex]));
			w->tagGroups[i].bamRecs[1] = bam_dup1(&(w->bamRecs[anchor2ReadIndex]));
			bam1_t *anchor1Bam = w->tagGroups[i].bamRecs[0];
			bam1_t *anchor2Bam = w->tagGroups[i].bamRecs[1];
			
			//int mapqState = (anchor1Bam->core.qual>=w->opt->minMAPQ && anchor2Bam->core.qual>=w->opt->minMAPQ) ? PF_STATES_PASS : PF_STATES_FAIL;
			int mapqState = w->pairQCCheck(anchor1Bam, anchor2Bam, &(w->opt->pairqc_opt));
			w->tagGroups[i].outputClass = (PF_STATES_PASS==mapqState) ? CPU_PAIR_OUTPUT_UU : CPU_PAIR_OUTPUT_UU_FAILED;
			
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
				int32_t tlen = getIntraSpan(anchor1Bam->core.pos, bam_calend(&(anchor1Bam->core), bam1_cigar(anchor1Bam)), anchor2Bam->core.pos, bam_calend(&(anchor2Bam->core), bam1_cigar(anchor2Bam)));
				if (anchor1Bam->core.pos<anchor2Bam->core.pos) {
					anchor1Bam->core.isize = tlen;
					anchor2Bam->core.isize = -1*tlen;
				} else {
					anchor1Bam->core.isize = -1*tlen;
					anchor2Bam->core.isize = tlen;
				}
				if (tlen<w->opt->selfLigation) w->pair_stats[tid].n_sameChromLow[mapqState]++; // STATISTICS: update the count
				else w->pair_stats[tid].n_sameChromHigh[mapqState]++; // STATISTICS: update the count
				w->pair_stats[tid].n_sameChrom[mapqState]++; // STATISTICS: update the count
			} else {
				// inter-chromosomal, use 0 as template length
				anchor1Bam->core.isize = anchor2Bam->core.isize = 0; //TLEN
				w->pair_stats[tid].n_differentChrom[mapqState]++; // STATISTICS: update the count
			}
			
			// introduce a YT:Z tag so that we know the decision
			uint8_t *yt = bam_aux_get(anchor1Bam, G_PAIRING_TAG);
			if (0!=yt) bam_aux_del(anchor1Bam, yt); // have to remove
			bam_aux_append(anchor1Bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_UU);
			yt = bam_aux_get(anchor2Bam, G_PAIRING_TAG);
			if (0!=yt) bam_aux_del(anchor2Bam, yt); // have to remove
			bam_aux_append(anchor2Bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_UU);
		} else {
			// discordant for one anchor or both the anchors
			// we cannot use this, but we keep them in separate file for investigation
			// introduce a YT:Z tag so that we know the decision
			w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_UU_DISCORDANT;
			w->pair_stats[tid].n_discordant++; // STATISTICS: update the count
			for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
				bam1_t *bam = &(w->bamRecs[k]);
				uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
				if (0!=yt) bam_aux_del(bam, yt); // have to remove
				bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_UU_Discordant);
			}
		}
	} else if (CPU_MAPSTATE_UNIQUE==anchor1MapState || CPU_MAPSTATE_UNIQUE==anchor2MapState) {
		// [UxxU] only 1 unique, output all the tags
		// introduce a YT:Z tag so that we know the decision
		w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_UxxU;
		const char *ytCode = (CPU_MAPSTATE_UNIQUE==anchor1MapState) ? G_PAIRING_TAG_Ux : G_PAIRING_TAG_xU;
		for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
			bam1_t *bam = &(w->bamRecs[k]);
			uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
			if (0!=yt) bam_aux_del(bam, yt); // have to remove
			bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)ytCode);
		}
	} else {
		// [xx], [nn] separate into two classes, output all the tags
		if (CPU_MAPSTATE_REPEAT==anchor1MapState || CPU_MAPSTATE_REPEAT==anchor2MapState) {
			// [xx]
			// introduce a YT:Z tag so that we know the decision
			w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_xx;
			for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
				bam1_t *bam = &(w->bamRecs[k]);
				uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
				if (0!=yt) bam_aux_del(bam, yt); // have to remove
				bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_xx);
			}
		} else {
			// [nn]
			// introduce a YT:Z tag so that we know the decision
			w->tagGroups[i].outputClass = CPU_PAIR_OUTPUT_nn;
			for(j=0, k=w->tagGroups[i].readIndex; j<w->tagGroups[i].n_reads; ++j, ++k) {
				bam1_t *bam = &(w->bamRecs[k]);
				uint8_t *yt = bam_aux_get(bam, G_PAIRING_TAG);
				if (0!=yt) bam_aux_del(bam, yt); // have to remove
				bam_aux_append(bam, G_PAIRING_TAG, G_PAIRING_TAG_TYPE, G_PAIRING_TAG_LEN, (uint8_t *)G_PAIRING_TAG_nn);
			}
		}
	}
}

void pair_process_alns(const pair_opt_t *opt, int64_t n_processed, int n, bam1_t *bamRecs, cpu_readid_t *readids, int *pNumTagGroups, tag_group_t *tagGroups, CPUBamBuffer_t* bamsBuffer, pair_stat_t *pair_stats, pair_stat_t *pair_stat, int (*pairQCCheck)(const bam1_t *,const bam1_t *,const pairqc_opt_t *))
{
	int i, j, k;
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	
	if (n<=0) return;
		
	w.opt = opt;
	w.n_processed = n_processed; // TODO:
	// TODO: w.n_processed = (opt->flag&MEM_F_PE)? n_processed>>1 : n_processed;
	w.bamRecs = bamRecs;
	w.readids = readids;
	w.pair_stats = pair_stats;
	w.pairQCCheck = pairQCCheck;
	
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
	kt_for(opt->n_threads, pair_worker, &w, numTagGroups);
	
	// accumulate the statistics
	for(j=0; j<CPU_MAPSTATE_COUNT; ++j) {
		for(k=0; k<CPU_MAPSTATE_COUNT; ++k) {
			for(i=0; i<opt->n_threads; ++i) {
				pair_stat->counts[j][k] += pair_stats[i].counts[j][k];
			}
		}
	}
	for(i=0; i<opt->n_threads; ++i) {
		pair_stat->n_discordant += pair_stats[i].n_discordant;
		pair_stat->n_differentChrom[PF_STATES_PASS] += pair_stats[i].n_differentChrom[PF_STATES_PASS];
		pair_stat->n_differentChrom[PF_STATES_FAIL] += pair_stats[i].n_differentChrom[PF_STATES_FAIL];
		pair_stat->n_sameChrom[PF_STATES_PASS] += pair_stats[i].n_sameChrom[PF_STATES_PASS];
		pair_stat->n_sameChrom[PF_STATES_FAIL] += pair_stats[i].n_sameChrom[PF_STATES_FAIL];
		pair_stat->n_sameChromLow[PF_STATES_PASS] += pair_stats[i].n_sameChromLow[PF_STATES_PASS];
		pair_stat->n_sameChromLow[PF_STATES_FAIL] += pair_stats[i].n_sameChromLow[PF_STATES_FAIL];
		pair_stat->n_sameChromHigh[PF_STATES_PASS] += pair_stats[i].n_sameChromHigh[PF_STATES_PASS];
		pair_stat->n_sameChromHigh[PF_STATES_FAIL] += pair_stats[i].n_sameChromHigh[PF_STATES_FAIL];
	}
}

int main_pair(int argc, char *argv[])
{
	int c, n, i, j, k;
	int is_count = 0;
	pair_opt_t *opt;
	
	bam1_t *bamRecs;
	CPUBamBuffer_t bamsBuffer;
	
	double t_diff;
	double t_timepointIO, t_timepointProcessing;
	double t_diffIO, t_diffProcessing;
	
	int64_t n_processed = 0;
	int64_t n_pairProcessed = 0;

	int is_bamin = 1;
	char in_mode[5];
	char *fn_list = 0;
	samfile_t *in = 0;
	
	int compress_level = -1;
	int is_header = 0;
	int is_bamout = 1;
	char out_mode[5];
	samfile_t *outfds[CPU_PAIR_OUTPUT_COUNT] = {0};
	
	int ret = 0;
	
	int (*pairQCCheck)(const bam1_t *,const bam1_t *,const pairqc_opt_t *) = getPairQCStateByUniqueLoci;
	
	init_CPUBamBuffer(&bamsBuffer);
	
	opt = pair_opt_init();
	strcpy(in_mode, "r");
	strcpy(out_mode, "w");
	while ((c = getopt(argc, argv, "Srnpcs:t:O:q:a:")) >= 0) {
		if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'S') is_bamin = 0;
		else if (c == 'c') is_count = 1;
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'r') opt->pairFlag = PAIR_F_REPEAT;
		else if (c == 'n') opt->pairFlag = PAIR_F_NOT_MAPEED;
		else if (c == 'O') { ks_set(optarg, &(opt->outputPrefix)); }
		else if (c == 's') opt->selfLigation = atoi(optarg), opt->selfLigation = opt->selfLigation > 0 ? opt->selfLigation : G_SELF_LIGATION;
		else if (c == 'q') {
			if (strstr(optarg, ".")) {
				opt->pairqc_opt.minFractionalMAPQ = atof(optarg);
				opt->pairqc_opt.minMAPQ = -1;
				pairQCCheck = getPairQCStateByUniqueLoci;
			}
			else {
				opt->pairqc_opt.minMAPQ = atoi(optarg);
				opt->pairqc_opt.minFractionalMAPQ = -1.0;
				pairQCCheck = getPairQCStateByMapQ;
			}
		}
		else if (c == 'a') opt->pairqc_opt.maxX1 = atoi(optarg);
		else {
			pair_opt_terminate(opt);
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
		fprintf(stderr, "Usage: cpu pair [options] <in.sam/.bam>\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -S         sam input\n");
		fprintf(stderr, "       -p         paired-end data\n");
		fprintf(stderr, "       -O         output prefix\n");
		fprintf(stderr, "       -c         output counting rather than bam files\n");
		fprintf(stderr, "       -s INT     self-ligation distance [%d bp]\n", opt->selfLigation);
		fprintf(stderr, "       -q INT     minimum MAPQ [>=#] or second(score)/best(score)<%.2f\n", opt->pairqc_opt.minFractionalMAPQ);
		fprintf(stderr, "       -a INT     max. number of secondary alternative allowed [<=%d]\n", opt->pairqc_opt.maxX1);
		fprintf(stderr, "\n");
		fprintf(stderr, "       -r         include repeat for pairing consideration\n");
		fprintf(stderr, "       -n         include not mapped for pairing consideration\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		
		pair_opt_terminate(opt);
		free(opt);
		return 1;
	}

	// open file handlers
	if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open \"%s\" for reading.\n", __func__, argv[optind]);
		
		free(fn_list);
		pair_opt_terminate(opt);
		free(opt);
		
		return 1;
	}
	if (in->header == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to read the header from \"%s\".\n", __func__, argv[optind]);
		
		free(fn_list);
		// close files, free and return
		samclose(in);
		
		pair_opt_terminate(opt);
		free(opt);
		
		return 1;
	}

	if (!is_count) {
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
		if (init_CPPairs_Ouputs (opt->outputPrefix.s, out_mode, in->header, outfds)) {
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open for writing.\n", __func__);
			
			free(fn_list);
			// close files, free and return
			samclose(in);
			
			pair_opt_terminate(opt);
			free(opt);
			
			return 1;
		}
		
		{
			if (opt->n_threads > 1) samthreads(outfds[CPU_PAIR_OUTPUT_UU], opt->n_threads, 256);
			if (opt->n_threads > 1) samthreads(outfds[CPU_PAIR_OUTPUT_UxxU], opt->n_threads, 256);
			if (opt->n_threads > 1) samthreads(outfds[CPU_PAIR_OUTPUT_xx], opt->n_threads, 256);
			if (opt->n_threads > 1) samthreads(outfds[CPU_PAIR_OUTPUT_nn], opt->n_threads, 256);
			if (opt->n_threads > 1) samthreads(outfds[CPU_PAIR_OUTPUT_UU_DISCORDANT], opt->n_threads, 256);
			if (opt->n_threads > 1) samthreads(outfds[CPU_PAIR_OUTPUT_UU_FAILED], opt->n_threads, 256);
		}
	}
	
	// PROCESSING
	if (argc == optind + 1) { // convert/print the entire file
		pair_stat_t pair_stat; memset(&pair_stat, 0, sizeof(pair_stat_t));
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
			pair_stat_t *pair_stats = calloc(opt->n_threads, sizeof(pair_stat_t));
			int numTagGroups = 0;
			pair_process_alns(opt, n_processed, n, bamRecs, readids, &numTagGroups, tagGroups, &bamsBuffer, pair_stats, &pair_stat, pairQCCheck);
			t_timepointProcessing = realtime();
			t_diffProcessing = t_timepointProcessing - t_timepointIO;
			
			if (!is_count) {
				// write the bam results
				for(i=0; i<numTagGroups; ++i) {
					if (CPU_PAIR_OUTPUT_UU==tagGroups[i].outputClass) {
						samwrite(outfds[CPU_PAIR_OUTPUT_UU], tagGroups[i].bamRecs[0]);
						samwrite(outfds[CPU_PAIR_OUTPUT_UU], tagGroups[i].bamRecs[1]);
					} else if (CPU_PAIR_OUTPUT_UU_FAILED==tagGroups[i].outputClass) {
						samwrite(outfds[CPU_PAIR_OUTPUT_UU_FAILED], tagGroups[i].bamRecs[0]);
						samwrite(outfds[CPU_PAIR_OUTPUT_UU_FAILED], tagGroups[i].bamRecs[1]);
					} else {
						samfile_t *outfd = outfds[tagGroups[i].outputClass];
						for(j=0, k=tagGroups[i].readIndex; j<tagGroups[i].n_reads; ++j, ++k) {
							samwrite(outfd, &bamRecs[k]);
						}
					}
				}
				t_timepointIO = realtime();
				t_diffIO += (t_timepointIO - t_timepointProcessing);
			}
			
			// clean up
			free(pair_stats);
			for(i=0; i<n; ++i) free(bamRecs[i].data);
			free(bamRecs);
			free(readids);
			for(i=0; i<numTagGroups; ++i) {
				if (CPU_PAIR_OUTPUT_UU==tagGroups[i].outputClass || CPU_PAIR_OUTPUT_UU_FAILED==tagGroups[i].outputClass) {
					bam_destroy1(tagGroups[i].bamRecs[0]);
					bam_destroy1(tagGroups[i].bamRecs[1]);
				}
			}
			
			n_processed += (n-bamsBuffer.unprocessed);
			n_pairProcessed += numTagGroups;
			if (bwa_verbose >= 3) {
				t_diff = realtime() - G_t_real;
				fprintf(stderr, "[M::%s] %lld tags %lld pairs processed, %.0f tags/sec %.0f pairs/sec, %.2f min, i/o %.2f sec, processing %.2f sec..\n", __func__, n_processed, n_pairProcessed, 1.0*n_processed/t_diff, 1.0*n_pairProcessed/t_diff, t_diff/60.0, t_diffIO, t_diffProcessing);
			}
			
			t_timepointProcessing = realtime();
		}
		
		// report pairing statistics
		fprintf(stdout, "##CPU\t%s\n", CPU_version);
		{
			fprintf(stdout, "##COMMAND\t");
			for (i = 0; i < argc; ++i)
				fprintf(stdout, " %s", argv[i]);
			fprintf(stdout, "\n");
		}
		
		fprintf(stdout, ">>Library information\n");
		fprintf(stdout, "#Measure\tValue\n");
		fprintf(stdout, "Filename\t%s\n", argv[optind]);
		fprintf(stdout, ">>END\n");
		
		long lTotal = 0; for(i=0; i<CPU_MAPSTATE_COUNT; ++i) { for(j=0; j<CPU_MAPSTATE_COUNT; ++j) { lTotal += pair_stat.counts[i][j]; } }
		fprintf(stdout, ">>General Mapping Statistics\n");
		fprintf(stdout, "#Category\tCode\tCount\tPercent\n");
		long lValue = lTotal;
		long lDenominator = lTotal; if (0==lDenominator) lDenominator=1; // prevent divide by zero
		fprintf(stdout, "TOTAL\t\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		// report at general level
		lValue = pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_UNIQUE];
		fprintf(stdout, "Both ends mapped uniquely\tUU\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_REPEAT]+pair_stat.counts[CPU_MAPSTATE_REPEAT][CPU_MAPSTATE_UNIQUE]+pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_NOT_MAPPED]+pair_stat.counts[CPU_MAPSTATE_NOT_MAPPED][CPU_MAPSTATE_UNIQUE];
		fprintf(stdout, "One end mapped uniquely\tUx/xU\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_REPEAT]+pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_NOT_MAPPED];
		fprintf(stdout, "One end mapped uniquely, head\tUx\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_REPEAT][CPU_MAPSTATE_UNIQUE]+pair_stat.counts[CPU_MAPSTATE_NOT_MAPPED][CPU_MAPSTATE_UNIQUE];
		fprintf(stdout, "One end mapped uniquely, tail\txU\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_REPEAT][CPU_MAPSTATE_REPEAT]+pair_stat.counts[CPU_MAPSTATE_REPEAT][CPU_MAPSTATE_NOT_MAPPED]+pair_stat.counts[CPU_MAPSTATE_NOT_MAPPED][CPU_MAPSTATE_REPEAT];
		fprintf(stdout, "Both ends do not mapped uniquely\txx\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_NOT_MAPPED][CPU_MAPSTATE_NOT_MAPPED];
		fprintf(stdout, "Both ends do not mapped\tnn\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		fprintf(stdout, ">>END\n");
		
		// UU breakdown
		fprintf(stdout, ">>Unique Interactions Statistics\n");
		fprintf(stdout, "#Category\tCount\tPercent\n");
		lValue = pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_UNIQUE];
		fprintf(stdout, "Total\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.n_discordant;
		fprintf(stdout, "Discordant anchor\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		
		if (-1==opt->pairqc_opt.minFractionalMAPQ) {
			fprintf(stdout, "mapq>=%d\tCount\tPercent\tPassed\t%%Passed\tFailed\t%%Failed\n", opt->pairqc_opt.minMAPQ);
		} else {
			fprintf(stdout, "second/best<%.2f\tCount\tPercent\tPassed\t%%Passed\tFailed\t%%Failed\n", opt->pairqc_opt.minFractionalMAPQ);
		}
		lValue = pair_stat.n_differentChrom[PF_STATES_PASS]+pair_stat.n_differentChrom[PF_STATES_FAIL];
		fprintf(stdout, "Different chromosomes\t%ld\t%.2f%%\t%ld\t%.2f%%\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator,
				pair_stat.n_differentChrom[PF_STATES_PASS], pair_stat.n_differentChrom[PF_STATES_PASS]*100.0f/lDenominator,
				pair_stat.n_differentChrom[PF_STATES_FAIL], pair_stat.n_differentChrom[PF_STATES_FAIL]*100.0f/lDenominator);
		lValue = pair_stat.n_sameChrom[PF_STATES_PASS]+pair_stat.n_sameChrom[PF_STATES_FAIL];
		fprintf(stdout, "Same chromosome\t%ld\t%.2f%%\t%ld\t%.2f%%\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator,
				pair_stat.n_sameChrom[PF_STATES_PASS], pair_stat.n_sameChrom[PF_STATES_PASS]*100.0f/lDenominator,
				pair_stat.n_sameChrom[PF_STATES_FAIL], pair_stat.n_sameChrom[PF_STATES_FAIL]*100.0f/lDenominator);
		// adjust for same chromosome breakdown
		lDenominator = pair_stat.n_sameChrom[PF_STATES_PASS]+pair_stat.n_sameChrom[PF_STATES_FAIL]; if (0==lDenominator) lDenominator=1; // prevent divide by zero
		kstring_t condition = {0,0,0};
		if (0==(opt->selfLigation%1000)) {
			ksprintf(&condition, "%d Kb", opt->selfLigation/1000);
		} else if (0==(opt->selfLigation%100)) {
			ksprintf(&condition, "%.1f Kb", opt->selfLigation/100.0);
		} else {
			ksprintf(&condition, "%d b", opt->selfLigation);
		}
		lValue = pair_stat.n_sameChromLow[PF_STATES_PASS]+pair_stat.n_sameChromLow[PF_STATES_FAIL];
		fprintf(stdout, "Same chromosome (<%s)\t%ld\t%.2f%%\t%ld\t%.2f%%\t%ld\t%.2f%%\n", condition.s, lValue, lValue*100.0f/lDenominator,
				pair_stat.n_sameChromLow[PF_STATES_PASS], pair_stat.n_sameChromLow[PF_STATES_PASS]*100.0f/lDenominator,
				pair_stat.n_sameChromLow[PF_STATES_FAIL], pair_stat.n_sameChromLow[PF_STATES_FAIL]*100.0f/lDenominator);
		lValue = pair_stat.n_sameChromHigh[PF_STATES_PASS]+pair_stat.n_sameChromHigh[PF_STATES_FAIL];
		fprintf(stdout, "Same chromosome (>=%s)\t%ld\t%.2f%%\t%ld\t%.2f%%\t%ld\t%.2f%%\n", condition.s, lValue, lValue*100.0f/lDenominator,
				pair_stat.n_sameChromHigh[PF_STATES_PASS], pair_stat.n_sameChromHigh[PF_STATES_PASS]*100.0f/lDenominator,
				pair_stat.n_sameChromHigh[PF_STATES_FAIL], pair_stat.n_sameChromHigh[PF_STATES_FAIL]*100.0f/lDenominator);
		free(condition.s);
		fprintf(stdout, ">>END\n");
		
		// report details
		lDenominator = lTotal; if (0==lDenominator) lDenominator=1; // prevent divide by zero
		fprintf(stdout, ">>Details Mapping Statistics\n");
		fprintf(stdout, "#Category\tCode\tCount\tPercent\n");
		lValue = pair_stat.counts[CPU_MAPSTATE_NOT_MAPPED][CPU_MAPSTATE_NOT_MAPPED];
		fprintf(stdout, "Unmappable\tNN\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_UNIQUE];
		fprintf(stdout, "Uniquely Mapped\tUU\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_REPEAT][CPU_MAPSTATE_REPEAT];
		fprintf(stdout, "Multiply Mapped\tMM\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_REPEAT];
		fprintf(stdout, "Multiply Mapped\tUM\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_REPEAT][CPU_MAPSTATE_UNIQUE];
		fprintf(stdout, "Multiply Mapped\tMU\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_UNIQUE][CPU_MAPSTATE_NOT_MAPPED];
		fprintf(stdout, "Head Only\tUN\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_REPEAT][CPU_MAPSTATE_NOT_MAPPED];
		fprintf(stdout, "Head Only\tMN\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_NOT_MAPPED][CPU_MAPSTATE_UNIQUE];
		fprintf(stdout, "Tail Only\tNU\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		lValue = pair_stat.counts[CPU_MAPSTATE_NOT_MAPPED][CPU_MAPSTATE_REPEAT];
		fprintf(stdout, "Tail Only\tNM\t%ld\t%.2f%%\n", lValue, lValue*100.0f/lDenominator);
		fprintf(stdout, ">>END\n");
		
	}
	// END - PROCESSING
	
	free(fn_list);
	
	// close files, free and return
	samclose(in);
		
	if (!is_count) {
		terminate_CPPairs_Ouputs(outfds);
	}
	terminate_CPUBamBuffer(&bamsBuffer);
	
	pair_opt_terminate(opt);
	free(opt);
	
	return ret;
}


