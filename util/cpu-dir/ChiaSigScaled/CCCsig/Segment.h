// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen

#ifndef GUARD_Segment
#define GUARD_Segment

#include <string>
#include <map> // QuantileMapper
#include <vector> // QuantileMapper
#include <stdlib.h> //abs
#include <stdint.h> //uint16_t
#include <fstream>

#include "bar.h"

//#define USE_SIMD_OPERATIONS 1
//#define FALSE_POSITIVE_SERIAL 1
//#define DUMP_CHR_NCHG 1
//#define DUMP_CHROMOSOME "chr18"
//#define DUMP_CHROMOSOME "chr1"

typedef uint16_t ChromosomeIndexType;

struct Segment {
	Segment();
	Segment(ChromosomeIndexType c, uint32_t s, uint32_t e);
	
	ChromosomeIndexType chr;
	uint32_t start;
	uint32_t end;
	uint32_t pos;
	
	bool operator < (const Segment& str) const {
		return (chr == str.chr ? pos < str.pos : chr < str.chr);
	}
	bool operator == (const Segment& str) const {
		return (chr == str.chr and start == str.start and end == str.end);
	}
};

struct SegmentMin {
	SegmentMin();
	SegmentMin(ChromosomeIndexType c, uint32_t s, uint32_t e);
	SegmentMin(const Segment &);
	
	ChromosomeIndexType chr;
	uint32_t start;
	uint32_t end;
};

struct Interaction {
	Interaction();
	Interaction(uint32_t c, bool m, double p);
	uint32_t mask:1, count:31;
	double pvalue;
};

struct DeltaCountSums {
	//DeltaCountSums():ulTotalDeltas(0),ulTotalMasked(0),nMaxSpan(0),pDeltas(NULL),pSums(NULL),deltaPresences(NULL) { }
	DeltaCountSums():ulTotalDeltas(0),nMaxSpan(0),pDeltas(NULL),pSums(NULL),deltaPresences(NULL) { }
	~DeltaCountSums() { terminate(); }
	
	void init(uint32_t maxSpan)
	{
		nMaxSpan = maxSpan;
		pDeltas = (uint32_t*) calloc(nMaxSpan+1, sizeof(uint32_t));
		pSums = (uint32_t*) calloc(nMaxSpan+1, sizeof(uint32_t));
		deltaPresences = barcreate(nMaxSpan+1);
	}
	
	void terminate()
	{
		if(pDeltas) { free(pDeltas); pDeltas=NULL; }
		if(pSums) { free(pSums); pSums=NULL; }
		if(deltaPresences) { bardestroy(deltaPresences); deltaPresences=NULL; }
		ulTotalDeltas = 0;
		//ulTotalMasked = 0;
		nMaxSpan = 0;
	}
	
	uint64_t ulTotalDeltas;
	//uint64_t ulTotalMasked;
	uint32_t nMaxSpan;
	uint32_t* pDeltas;
	uint32_t* pSums;
	BIT_ARRAY *deltaPresences;
};

struct DeltaStatistics {
	DeltaStatistics();
	DeltaStatistics(uint64_t s, uint64_t c);
	double getMean() { return (sum*1.0/count); }
	uint64_t sum;
	uint64_t count;
};

struct QuantileMapper {
	QuantileMapper() { }
	void computeDeltaCountQuantile(uint64_t, uint32_t*, uint32_t, std::ofstream*);
	std::vector<uint32_t>& getQuantileDeltas();
	std::vector<uint32_t>& getQuantiles();
	uint32_t getDeltaQuantile(uint32_t);
	uint32_t getDeltaQuantileIdx(uint32_t);
	std::vector<uint32_t> m_deltas;
	std::vector<uint32_t> m_quantiles;
};

#endif
