#ifndef CPU_TAG_BAM_H
#define CPU_TAG_BAM_H

#include "bam.h"

#define G_MAPPROGRAM_TAG      "XL"
#define G_MAPPROGRAM_TAG_TYPE 'Z'
// need to include the null-terminator, thus 3+1
#define G_MAPPROGRAM_TAG_LEN  4
#define G_MAPPROGRAM_TAG_BWAALN  "aln"
#define G_MAPPROGRAM_TAG_BWAMEM  "mem"
#define G_TOPHIT_TAG				"X0"
#define G_ALTERNATIVEHIT_TAG		"X1"
#define G_TOPHIT_SCORE_TAG			"Y0"
#define G_ALTERNATIVEHIT_SCORE_TAG	"Y1"


#define G_PAIRING_TAG      "YT"
#define G_PAIRING_TAG_TYPE 'Z'
// need to include the null-terminator, thus 2+1
#define G_PAIRING_TAG_LEN  3
#define G_PAIRING_TAG_UU   "UU"
#define G_PAIRING_TAG_UU_Discordant   "ud"
#define G_PAIRING_TAG_Ux   "Ux"
#define G_PAIRING_TAG_xU   "xU"
#define G_PAIRING_TAG_xx   "xx"
#define G_PAIRING_TAG_nn   "nn"

#define G_SELF_LIGATION    3000

typedef struct {
	int     readIndex;
	uint8_t secondaryTag:1, n_reads:7;
	
	uint8_t releaseBamRecs:1, finalLinker:3, outputClass:4;
	bam1_t  *bamRecs[2];
	uint32_t endPoss[2];
} tag_group_t;

typedef struct {
	int unprocessed;
	bam1_t *unprocessedBamRecs;
	int lastBlock;
} CPUBamBuffer_t;

static inline void init_CPUBamBuffer(CPUBamBuffer_t *buffer)
{
	buffer->unprocessed = 0;
	buffer->unprocessedBamRecs = NULL;
	buffer->lastBlock = 0;
}

static inline void add_CPUBamBuffer(bam1_t *bams, int nBams, CPUBamBuffer_t *buffer)
{
	int i;
	int finalCount = buffer->unprocessed + nBams;
	bam1_t *newBamRecs = realloc(buffer->unprocessedBamRecs, finalCount*sizeof(bam1_t));
	if (!newBamRecs) {
		// TODO: error! out of memory
		exit(ENOMEM);
	} else {
		buffer->unprocessedBamRecs = newBamRecs;
		newBamRecs = &(buffer->unprocessedBamRecs[buffer->unprocessed]);
		memset(newBamRecs, 0, nBams*sizeof(bam1_t));
		for(i=0; i<nBams; ++i) {
			bam_copy1(&(newBamRecs[i]), &(bams[i]));
		}
		buffer->unprocessed = finalCount;
	}
}

static inline void terminate_CPUBamBuffer(CPUBamBuffer_t *buffer)
{
	int i;
	for(i=0; i<buffer->unprocessed; ++i) {
		free(buffer->unprocessedBamRecs[i].data);
	}
	
	buffer->unprocessed = 0;
	if (buffer->unprocessedBamRecs) {
		free(buffer->unprocessedBamRecs);
		buffer->unprocessedBamRecs = NULL;
	}
	buffer->lastBlock = 0;
}

//static inline bam1_t *readCPUBam(int *nread, int n_, samfile_t *sbin, int nUnprocessed, bam1_t *unprocessedBamRecs, int *lastBlock)
static inline bam1_t *readCPUBam(int *nread, int n_, samfile_t *sbin, CPUBamBuffer_t *buffer)

{
	int i;
	int start = buffer->unprocessed;
	bam1_t *bams = 0;
	bams = realloc(bams, n_ * sizeof(bam1_t));
	memset(bams, 0, n_*sizeof(bam1_t));
	*nread = 0;
	
	buffer->lastBlock = 0;
	if (buffer->unprocessed>0) {
		// copy from cached buffer first!
		for(i=0; i<buffer->unprocessed; ++i) {
			*nread = i+1;
			bam_copy1(&bams[i], &(buffer->unprocessedBamRecs[i]));
			free(buffer->unprocessedBamRecs[i].data);
		}
		buffer->unprocessed = 0;
		free(buffer->unprocessedBamRecs); buffer->unprocessedBamRecs = NULL;
	}
	
	// read from the file
	for(i=start; i<n_; ++i) {
		int r = samread(sbin, &bams[i]);
		if (r>0) {
			*nread = i+1;
		} else if (-1==r) {
			// EOF
			buffer->lastBlock = 1;
			break;
		} else if (r < -1) {
			if (bwa_verbose >= 2) fprintf(stderr, "[W::%s] truncated file.\n", __func__);
			buffer->lastBlock = 1;
			break;
		}
	}
	if (0==*nread) {
		free(bams); bams=0;
	}
	return bams;
}


#endif // CPU_TAG_BAM_H
