// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen


#ifndef GUARD_CCCMatrix
#define GUARD_CCCMatrix

#include "Segment.h"
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <limits>
#include <iostream>

#include "../spline/spline.h" // Spline-smoothing

#include "bar.h"

const uint32_t MAXDELTA = std::numeric_limits<uint32_t>::max();

typedef int32_t SegmentKey;
enum IndexType {NONINDEX, INTRAINDEX, ROWINDEX, COLINDEX};

#define CACHE_FLAG_N       0x00000001
#define CACHE_FLAG_EXPECTN 0x00000002

typedef std::map<SegmentKey, uint32_t > SegmentKeyMapType;
typedef std::vector<SegmentKeyMapType > SegmentKeyMapMatrixType;

typedef std::vector<SegmentKey > SegmentKeyRowType;
typedef std::vector<SegmentKeyRowType > SegmentKeyMatrixType;

typedef std::vector<uint32_t > SegmentIterType;
typedef std::vector<uint32_t > CountRowType;
typedef std::vector<CountRowType > CountMatrixType;
typedef std::vector<double > PvalueRowType;
typedef std::vector<PvalueRowType > PvalueMatrixType;
typedef std::vector<bool > MaskRowType;
typedef std::vector<MaskRowType > MaskMatrixType;

typedef std::vector<double > ExpectRowType;
typedef std::vector<uint32_t > MarginalRowType;

class CCCMatrixInteraction {
public:
	CCCMatrixInteraction();
	//CCCMatrixInteraction(std::set<Segment>&, std::set<Segment>&, Interaction defaultvalue, bool, std::string& chr);
	~CCCMatrixInteraction();
	void LoadInteractions(std::set<Segment> &rows, std::set<Segment> &cols, Interaction defaultvalue, bool isintra, std::string& chr);
	void InteractionLoaded(uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/);
	void printInteractionMatrix(std::ostream&, bool hideMask=false);
	void printPvalueMatrix(std::ostream&, bool hideMask=false);

	void setElement(Segment&, Segment&, uint32_t);
	void getRowSums(CountRowType &);
	
	// TODO: decide whether these are external aux function or members
	void getStatisticsPerDelta(std::vector<DeltaStatistics >&, uint32_t*, uint32_t minDelta=8000, uint32_t maxDelta=MAXDELTA, bool masking=false);

	uint64_t getNumberOfPositives(double alpha);
	uint64_t getNumberOfPositives(double alpha, uint32_t cutoff);
	void getPositivesBound(uint32_t cutoff, uint64_t& cutoffMaxPositive);
	uint64_t maskByAlpha(double alpha, uint32_t cutoff, bool dump=false);
	void calculatePvalues(uint32_t minDelta=8000, uint32_t maxDelta=MAXDELTA);
	void printPositives(std::ostream& ostr, std::string&, double alpha, uint32_t cutoff, bool printNij, bool printAll=false);
	// TODO: END - decide whether these are external aux function or members
	
	uint64_t maskDeltasCountSum(uint32_t*, uint32_t*, BIT_ARRAY*, uint32_t, uint32_t minDelta /* = 0 */, uint32_t maxDelta /*=MAXDELTA*/);
	uint64_t getDeltasCountSum(uint32_t*, uint32_t*, BIT_ARRAY*, uint32_t, uint32_t minDelta /* = 0 */, uint32_t maxDelta /*=MAXDELTA*/, bool masking /* =false */);
	uint64_t getDeltasCount(uint32_t*, BIT_ARRAY*, uint32_t, uint32_t minDelta /* = 0 */, uint32_t maxDelta /*=MAXDELTA*/, bool masking /* =false */);
	uint32_t getN(uint32_t minDelta=0, uint32_t maxDelta=MAXDELTA) { return cachedN; }
	bool isIntra;
	
	// TODO: we are keeping the expectation values
	void setDeltaToExpectationMapper(spline& smoothDelta2ExpMapper, uint32_t minDelta=8000, uint32_t maxDelta=MAXDELTA);
	double getExpectN(uint32_t minDelta=0, uint32_t maxDelta=MAXDELTA) { return cachedExpectN; }
	double getExpectElement(SegmentKey rowKey, SegmentKey colKey) {
		// TODO:
		int delta = abs(orderedSegmentMidPoints[colKey]-orderedSegmentMidPoints[rowKey]);
		if(delta>=minExpDelta and delta<=maxExpDelta) return pDelta2ExpSpline->cachedAt(delta);
		else return 0.0;
	}
	void printExpectationMatrix(std::ostream&, bool hideMask=false);
	void getExpectRowSums(std::vector<double > &);
	// TODO: END - we are keeping the expectation values
	
	// TODO: for FDR
	double estimateFalsePositives(double alpha, uint32_t cutoff, uint32_t minDelta=8000, uint32_t maxDelta=MAXDELTA);
	uint32_t getSegmentCount();
	// TODO: END - for FDR
		
	// TODO: assess caching effect
	void cacheRowSums();
	// TODO: END - assess caching effect
	
	unsigned long getInteractionCount();
	void setThreads(int n);
	
	// for deltas managed memory
	uint32_t getWidestInteractionSpan();
	// END - for deltas managed memory
	
	void setDepthThreshold (float t) { fDepthThreshold=t; }
	
private:
	Interaction defaultVal;
	std::map<Segment, SegmentKey > segmentToKey;
	std::vector<SegmentMin> orderedSegments;
	uint32_t* orderedSegmentMidPoints;
	
	SegmentKeyMapMatrixType mat_skToIndex;
	SegmentKeyMatrixType mat_sk;
	CountMatrixType mat_count;
	PvalueMatrixType mat_pvalue;
	MaskMatrixType mat_mask;

	// TODO: for reporting
	uint64_t m_TotalInteractions;
	
	// TODO: assess caching effect
	ExpectRowType m_getExp;
	MarginalRowType m_getMarginal;
	// TODO: END - assess caching effect

	
	// TODO: we are keeping the expectation values
	spline* pDelta2ExpSpline;
	uint32_t minExpDelta;
	uint32_t maxExpDelta;
	// TODO: END - we are keeping the expectation values
	
	int nThreads;
	
	// TODO: more caching flags
	uint32_t cacheFlag;
	uint32_t cachedN;
	void calculate_cache_N(uint32_t minDelta=0, uint32_t maxDelta=MAXDELTA);
	double cachedExpectN;
	void calculate_cache_ExpectN(uint32_t minDelta=0, uint32_t maxDelta=MAXDELTA);
	
	// TODO: for easier debugging
	std::string m_chr;

	std::vector<SegmentKey > segmentKeyStarts;
	std::vector<SegmentKey > segmentKeyEnds;
	std::vector<SegmentKey > segmentPrefixKeyStarts;
	std::vector<SegmentKey > segmentPrefixKeyEnds;
	std::vector<uint32_t > rowMasked;

	SegmentIterType segmentIterStarts;
	SegmentIterType segmentIterEnds;
	
	float fDepthThreshold;
};

#endif
