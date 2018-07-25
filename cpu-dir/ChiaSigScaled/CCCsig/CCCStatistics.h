// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen

#ifndef GUARD_CCCStatistics
#define GUARD_CCCStatistics

#include "CCCMatrix.h"
#include "Segment.h"

#include <map>
#include <vector>

#include <algorithm>
#include <numeric> 

#include "../spline/spline.h" // Spline-smoothing

void setThreads(int n);

void printPositives(std::ostream& ostr, std::string&, CCCMatrixInteraction&, double, uint32_t, bool printNij=false, bool printAll=false);

// Utils:
template<typename T> std::vector<T> getSequence(T, T, uint32_t);

std::pair<double, double> getAlphaRange(std::map<double, double> , double);

std::map<double, double>  estimateFDRAlpha(const char *szDataFile, std::vector<std::string>&, std::vector<std::pair<ChromosomeIndexType, uint64_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, double, uint32_t, bool noEarlyTermination=false, uint32_t cutoff=3, uint32_t deltaCutoff=8000, uint32_t maxDelta=MAXDELTA);

void computeInteractionsPvalues(std::vector<std::string>&, std::vector<std::pair<ChromosomeIndexType, uint64_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, spline& smoothDelta2ExpMapper, uint32_t deltaCutoff, uint32_t minDelta=8000, uint32_t maxDelta=MAXDELTA);

void estimateDeltaSpline(const char *szDataFile, spline&, std::vector<double>&, DeltaCountSums&, std::vector<std::string>&, std::vector<std::pair<ChromosomeIndexType, uint64_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, uint32_t minDelta=8000, uint32_t maxDelta=MAXDELTA, bool masking=false, uint32_t percentiles=1000, bool dumpSplines=false);

void estimateResources(std::string&, std::vector<std::string>&, std::vector<std::pair<ChromosomeIndexType, uint64_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, uint32_t minDelta=8000, uint32_t maxDelta=MAXDELTA, bool masking=false);

void getChromosomesDataSizeDesc(std::vector<std::string>& chromosomes, std::map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, std::vector<std::pair<ChromosomeIndexType, uint64_t > >&);

void sweepQunatilesDeltaSpline(std::string&, std::vector<std::string>&, std::vector<std::pair<ChromosomeIndexType, uint64_t > >&, std::map<ChromosomeIndexType, CCCMatrixInteraction>&, std::vector<uint32_t >&, uint32_t minDelta=8000, uint32_t maxDelta=MAXDELTA, bool masking=false);

#endif
