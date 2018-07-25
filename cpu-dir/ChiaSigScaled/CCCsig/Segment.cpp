// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen


#include "Segment.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <assert.h>

using namespace std;

extern "C" {
	double cputime();
	double realtime();
}

Segment::Segment()
:chr(0),start(0),end(0),pos(0)
{
}

Segment::Segment(ChromosomeIndexType c, uint32_t s, uint32_t e)
:chr(c), start(s), end(e), pos((s+e)/2)
{
}

SegmentMin::SegmentMin()
:chr(0),start(0),end(0)
{
}

SegmentMin::SegmentMin(ChromosomeIndexType c, uint32_t s, uint32_t e)
:chr(c), start(s), end(e)
{
}

SegmentMin::SegmentMin(const Segment &a)
:chr(a.chr), start(a.start), end(a.end)
{
}


Interaction::Interaction()
:count(0), mask(0), pvalue(1)
{
}

Interaction::Interaction(uint32_t c, bool m, double p)
:count(c), mask(m ? 1 : 0), pvalue(p)
{
}

DeltaStatistics::DeltaStatistics()
: sum(0), count(0)
{
}

DeltaStatistics::DeltaStatistics(uint64_t s, uint64_t c)
: sum(s), count(c)
{
}

// WCH: this version uses keeps the counting of the deltas instead of vector of deltas
void QuantileMapper::computeDeltaCountQuantile(uint64_t vecSize, uint32_t* vec, uint32_t nQuant, ofstream* os)
{
	m_deltas.clear(); m_deltas.reserve(nQuant+1);
	m_quantiles.clear(); m_quantiles.reserve(nQuant+1);
	
	// Function that computes quantiles of size quantSize and returns
	// a map from each integer element to its quantile.
	// Obs: If a given for delta is found in two adjacent divisions,
	// then the smallest delta is used.
	
	double t_diff;
	double t_start;
	
	uint64_t quantSize = vecSize/nQuant;
	//cout << "vecSize & quantSize: " <<  vecSize << " " << quantSize << endl;
	fprintf(stderr, "[M::%s:%d] vecSize=%llu, nQuant=%u, quantSize=%llu\n", __func__, __LINE__, vecSize, nQuant, quantSize);
	assert(vecSize > quantSize);
	
	if (os) {
		*os << endl;
		*os << "Delta Count Qunatile" << endl;
		*os << "vecSize:" << vecSize << endl;
		*os << "nQuant:" << nQuant << endl;
		*os << "quantSize:" << quantSize << endl;
		*os << "s\te\tqsize\tmid\tsdelta\tedelta\tmiddelta\tmedian\tlast" << endl;
	}

	t_start = realtime();
	uint64_t nTotalQuantSize = quantSize*nQuant;
	int nProcessedWRTdata = 0;
	if (nTotalQuantSize>vecSize) nProcessedWRTdata = 1;
	else if (nTotalQuantSize<vecSize) nProcessedWRTdata = -1;
	// 0==nProcessedWRTdata : exact bin size matched
	// nProcessedWRTdata < 0 : i.e. quantSize*nQuant < vecSize
	//                         we have some spill over, which we will just used the last chunk media
	// nProcessedWRTdata > 0 : i.e. quantSize*nQuant > vecSize
	//                         last chunk will have insufficient element in bin
	uint64_t nEffQuant = (nProcessedWRTdata > 0) ? nQuant - 1 : nQuant;
	uint32_t median = 0;
	uint64_t i = 0, j=0;

	// loop thru' once just for the indexes so that we can remember the values that we need
	// then, we look up values from these needed position for the actual computation
	// we should also keep the 1st and last deltas
	vector<uint64_t> indices; //indices.push_back(0);
	for (uint64_t nQuantIndex=0; nQuantIndex < nEffQuant; ++nQuantIndex) {
		indices.push_back(i); // before debugging dump
		uint64_t midpos=i+(quantSize)/2;
		if (quantSize % 2) {
			indices.push_back(midpos);
		} else {
			indices.push_back(midpos-1);
			indices.push_back(midpos);
		}
		
		j = i + (quantSize-1);
		indices.push_back(j);
		
		i+=quantSize;
	}
	// Last remaindig chunk of values:
	if (0!=nProcessedWRTdata) {
		indices.push_back(i); // before debugging dump
		uint64_t tmp = quantSize;
		quantSize=(vecSize-i); // quantSize of final chunk
		j = i + (quantSize-1);
		
		if (0>nProcessedWRTdata) {
			// handle the spill over
			// WCH: IMPORTANT: we are mimicking the original implementation!!!
			//      Since the map is "full", set to the last value
			
			// optimized implementation
			indices.push_back(j);
		} else /*if (0<nProcessedWRTdata)*/ {
			// handle insufficient data in last bin
			uint64_t midpos=i+(quantSize)/2;
			if (quantSize % 2) {
				indices.push_back(midpos);
			} else {
				indices.push_back(midpos-1);
				indices.push_back(midpos);
			}
			
			indices.push_back(j);
		}
		quantSize = tmp;
	}
	indices.push_back(vecSize-1);
	
	// we now go thru' the deltasCount vector to cache the values at the required indices
	sort(indices.begin(), indices.end());
	map<uint64_t, uint32_t > vecMap; // <prefix sum> --> <delta>
	uint64_t nCurrCount = 0;
	uint64_t indicesIdx = 0;
	for(uint32_t i=0; nCurrCount<vecSize && indicesIdx<indices.size(); ++i) {
		if (vec[i]>0) {
			nCurrCount += vec[i];
			// (nCurrCount-1) for getting the last index, as index is 0-based
			while ((nCurrCount-1)>=indices[indicesIdx]) {
				// passed a band
				vecMap[indices[indicesIdx]] = i;
				indicesIdx++;
				if (indicesIdx>=indices.size()) {
					break;
				}
			}
		}
	}
	// END - hacked
	/*
	for(vector<unsigned long>::iterator it=indices.begin(); it!=indices.end(); ++it) {
		cerr << (*it) << "\t" << vecMap[*it] << endl;
	}
	 */
	
	// re-initialized as the previous section has tamed the variable
	i = 0;
	// let's build the look up table
	for (uint64_t nQuantIndex=0; nQuantIndex < nEffQuant; ++nQuantIndex) {
		uint64_t midpos=i+(quantSize)/2;
		median = (quantSize % 2) ? vecMap[midpos] : (vecMap[midpos-1]+vecMap[midpos])/2;
		
		j = i + (quantSize-1);
		if (0==m_deltas.size()) {
			m_deltas.push_back(vecMap[j]);
			m_quantiles.push_back(median);
		} else if (m_deltas.back()!=vecMap[j]) {
			m_deltas.push_back(vecMap[j]);
			m_quantiles.push_back(median);
		}
		
		if (os) {
			*os << i << "\t" << j << "\t" << quantSize << "\t" << midpos << "\t" << vecMap[i] << "\t" << vecMap[j] << "\t";
			if (quantSize % 2) *os << vecMap[midpos];
			else  *os << vecMap[midpos-1] << "," << vecMap[midpos];
			*os << "\t" << median << "\tn.a." << endl;
		}
		
		i+=quantSize;
	}
	
	uint32_t lastval=median;
	// Last remaindig chunk of values:
	if (0!=nProcessedWRTdata) {
		quantSize=(vecSize-i); // quantSize of final chunk
		j = i + (quantSize-1);
		
		if (0>nProcessedWRTdata) {
			// handle the spill over
			// WCH: IMPORTANT: we are mimicking the original implementation!!!
			//      Since the map is "full", set to the last value
			
#if 0
			// proper implementation
			if (m_deltas.back()!=vecMap[j]) {
				m_deltas.push_back(vecMap[j]);
				m_quantiles.push_back(lastval);
			}
#else
			// optimized implementation
			if (0==m_deltas.size()) {
				m_deltas.push_back(vecMap[j]);
				m_quantiles.push_back(lastval);
			} else {
				m_deltas.back() = vecMap[j];
			}
#endif
			
			if (os) {
				*os << i << "\t" << j << "\t" << quantSize << "\tsee prev\t" << vecMap[i] << "\t" << vecMap[j] << "\tsee prev\t" << lastval << "\t*" << endl;
			}
			
		} else /*if (0<nProcessedWRTdata)*/ {
			// handle insufficient data in last bin
			uint64_t midpos=i+(quantSize)/2;
			median = (quantSize % 2) ? vecMap[midpos] : (vecMap[midpos-1]+vecMap[midpos])/2;
			
			if (0==m_deltas.size()) {
				m_deltas.push_back(vecMap[j]);
				m_quantiles.push_back(median);
			} else if (m_deltas.back()!=vecMap[j]) {
				m_deltas.push_back(vecMap[j]);
				m_quantiles.push_back(median);
			}
			
			if (os) {
				*os << i << "\t" << j << "\t" << quantSize << "\t" << midpos << "\t" << vecMap[i] << "\t" << vecMap[j] << "\t";
				if (quantSize % 2) *os << vecMap[midpos];
				else *os << vecMap[midpos-1] << "," << vecMap[midpos];
				*os << "\t" << median << "\tn.a." << endl;
			}
		}
	}
	
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] setup took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
}

vector<uint32_t>& QuantileMapper::getQuantileDeltas()
{
	return m_deltas;
}

vector<uint32_t>& QuantileMapper::getQuantiles()
{
	return m_quantiles;
}

uint32_t QuantileMapper::getDeltaQuantile(uint32_t delta)
{
	// check look up table
	vector<uint32_t>::iterator it = lower_bound(m_deltas.begin(), m_deltas.end(), delta);
	
	// TODO: check if we need to difference to be divided by the size of the element?!!!
	uint64_t idx = it-m_deltas.begin(); // max(it-m_deltas.begin(), 0)
	return m_quantiles[idx];
}

uint32_t QuantileMapper::getDeltaQuantileIdx(uint32_t delta)
{
	// check look up table
	vector<uint32_t>::iterator it = lower_bound(m_deltas.begin(), m_deltas.end(), delta);
	
	// TODO: check if we need to difference to be divided by the size of the element?!!!
	uint64_t idx = it-m_deltas.begin(); // max(it-m_deltas.begin(), 0)
	return (uint32_t)idx;
}


