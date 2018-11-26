// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen


#include "CCCMatrix.h"
#include <assert.h>
#include <iostream>
#include <fstream>  // for logging
#include <numeric>  // for accumulate
#include "../stocc/stocc.h" // Non-central hypergeometric and Poisson
#include <string.h> // for memcpy

void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);

using namespace std;

const double fSmallValue = 1E-20;
const double dOmegaTablePrecision = 1E-53; //dOmegaPrecision * 0.001

const double dMinOmega=1E-01;

double getPvalNCHGNew(bool dump, uint32_t nij, uint32_t ni, uint32_t nj, uint32_t n, double odds)
{
	nij--;
	odds = max(odds, dMinOmega); //max(odds, minOmega)
	
#ifdef DUMP_CHR_NCHG
	if (dump) fprintf(stderr, "getPvalNCHGNew(nij:%d, ni:%d, nj:%d, N:%d, odds:%g, precision:%g)\n", nij, ni, nj, n, odds, precision);
#endif
	
	// TODO: in p-value calculation loop, it is 1e-20
	CFishersNCHypergeometric xhyp(nj, ni, n, odds/*, dOmegaPrecision*/);
	int BufferLength = xhyp.getTableLength();
#ifdef DUMP_CHR_NCHG
	if (dump) fprintf(stderr, "getPvalNCHGNew.MakeTable(len:%d)\n", BufferLength);
#endif
	double buffer[BufferLength];
	uint32_t x1, x2;
	double factor = 1. / xhyp.MakeTable(buffer, BufferLength, &x1, &x2, dOmegaTablePrecision); // make table of probability values
	uint32_t xmean = (uint32_t)(xhyp.mean() + 0.5);           // Round mean
	
	assert(xmean >= x1);
	assert(xmean <= x2);
	
	// NOTE: no need to re-assign, just use nij directly
	double dval;
	if (nij <= xmean) { // Left tail
		if (nij < x1) {
			dval = 1.; //p = 0.; // Outside table
		} else {
			double sum; int x;
			for (x = x1, sum = 0; x <= nij; x++) sum = buffer[x-x1] += sum; // Make left tail of table cumulative
			dval = 1. - buffer[nij-x1] * factor; //p = buffer[nij-x1] * factor; // Probability from table
		}
		// NOTE: false==lower_tail, thus, we NEED to invert above!!!
	} else {      // Right tail
		if (nij >= x2) {
			dval = 0.; //p = 0.; // Outside table
		} else {
			double sum; int x;
			for (x = x2, sum = 0; x > nij; x--) sum = buffer[x-x1] += sum;  // Make right tail of table cumulative from the right
			dval = buffer[nij-x1+1] * factor;  //p = buffer[nij-x1+1] * factor; // Probability from table
		}
		// NOTE: false==lower_tail, thus, we need NOT invert above
	}
#ifdef DUMP_CHR_NCHG
	if (dump) fprintf(stderr, "getPvalNCHGNew(nij:%d, ni:%d, nj:%d, N:%d, odds:%g, precision:%g) BufferLength:%d factor:%g returns %g\n", nij, ni, nj, n, odds, precision, BufferLength, factor, dval);
#endif
	return dval;
}
// END - WCH : OPTIMZED

// collapse of getQuantileNCHG()+getPvalNCHGNew()
double getQuantilePvalNCHG(bool dump, uint32_t nij, double p, uint32_t ni, uint32_t nj, uint32_t n, double odds)
{
	
#ifdef DUMP_CHR_NCHG
	if (dump)
		fprintf(stderr, "getQuantilePvalNCHG.getQuantileNCHG(p:%g, ni:%d, nj:%d, n:%d, odds:%g, precision:%g, lower_tail:%s)\n", p, ni, nj, n, odds, precision, lower_tail?"true":"false");
#endif
	
#if 0
	assert(p >= 0 and p <=1);
	assert(odds>=dMinOmega);
#endif
	
	CFishersNCHypergeometric xhyp(nj, ni, n, odds/*, dOmegaPrecision*/);
	int BufferLength = xhyp.getTableLength();
	if (0==BufferLength) return 0.;
	
#ifdef DUMP_CHR_NCHG
	if (dump)
		fprintf(stderr, "getQuantileNCHG.MakeTable(len:%d)\n", BufferLength);
#endif
	double buffer[BufferLength];
	uint32_t x1,x2;
	double factor = xhyp.MakeTable(buffer, BufferLength, &x1, &x2, dOmegaTablePrecision); // make table of probability values
	double secbuffer[BufferLength];
	memcpy(secbuffer, buffer, BufferLength*sizeof(double));
	
	uint32_t x;
	double sum;
	for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum; // Make table cumulative
	
	// NOTE: OPTIMIZED out as true==lower_tail  ==> this inversion does not need to happen
	p *= factor; // Table is scaled by factor
	
	// Binary search in table:
	unsigned int a, b, c; // Used in binary search
	a = 0; b = x2 - x1 + 1;
	while (a < b) {
		c = (a + b) / 2;
		if (p <= buffer[c]) {
			b = c;
		}
		else {
			a = c + 1;
		}
	}
	x = x1 + a;
	if (x > x2) x = x2;           // Prevent values > xmax that occur because of small imprecisions
#ifdef DUMP_CHR_NCHG
	if (dump)
		fprintf(stderr, "getQuantilePvalNCHG.getQuantileNCHG(p:%g, ni:%d, nj:%d, n:%d, odds:%g, precision:%g, lower_tail:%s) BufferLength:%d factor:%g returns %d\n", p, ni, nj, n, odds, precision, lower_tail?"true":"false", BufferLength, factor, x);
#endif
	
	// TODO: continue with getPvalNCHGNew()
	
	nij = max(nij, 1+x);
	nij--; //nij-1
	
	//odds>dMinOmega, used cached values
	factor = 1. / factor;
#ifdef DUMP_CHR_NCHG
	if (dump) fprintf(stderr, "getQuantilePvalNCHG.getPvalNCHGNew(nij:%d, ni:%d, nj:%d, N:%d, odds:%g, precision:%g)\n", nij, ni, nj, n, odds, precision);
#endif
	int xmean = (int)(xhyp.mean() + 0.5);           // Round mean
	
	assert(xmean >= x1);
	assert(xmean <= x2);
	
	// NOTE: no need to re-assign, just use nij directly
	double dval;
	if (nij <= xmean) { // Left tail
		if (nij < x1) {
			dval =  1.; //p = 0.; // Outside table
		} else {
			for (x = x1, sum = 0; x <= nij; x++) sum = secbuffer[x-x1] += sum; // Make left tail of table cumulative
			dval = 1. - secbuffer[nij-x1] * factor; //p = buffer[nij-x1] * factor; // Probability from table
		}
		// NOTE: false==lower_tail, thus, we NEED to invert above!!!
	} else {      // Right tail
		if (nij >= x2) {
			dval = 0.; //p = 0.; // Outside table
		} else {
			for (x = x2, sum = 0; x > nij; x--) sum = secbuffer[x-x1] += sum;  // Make right tail of table cumulative from the right
			dval = secbuffer[nij-x1+1] * factor;  //p = buffer[nij-x1+1] * factor; // Probability from table
		}
		// NOTE: false==lower_tail, thus, we need NOT invert above
	}
#ifdef DUMP_CHR_NCHG
	if (dump) fprintf(stderr, "getQuantilePvalNCHG.getPvalNCHGNew(nij:%d, ni:%d, nj:%d, N:%d, odds:%g, precision:%g) BufferLength:%d factor:%g returns %g\n", nij, ni, nj, n, odds, precision, BufferLength, factor, dval);
#endif
	return dval;
}

// original version
int32_t getQuantileNCHG(bool dump, double p, uint32_t ni, uint32_t nj, uint32_t n, double odds)
{
	
#ifdef DUMP_CHR_NCHG
	if (dump)
		fprintf(stderr, "getQuantileNCHG(p:%g, ni:%d, nj:%d, n:%d, odds:%g, precision:%g, lower_tail:%s)\n", p, ni, nj, n, odds, precision, lower_tail?"true":"false");
#endif
	
#if 0
	assert(p >= 0.0 and p <=1.0);
#endif
	
	CFishersNCHypergeometric xhyp(nj, ni, n, odds/*, dOmegaPrecision*/);
	uint32_t BufferLength = xhyp.getTableLength();
	if (0==BufferLength) return -1;
	
#ifdef DUMP_CHR_NCHG
	if (dump)
		fprintf(stderr, "getQuantileNCHG.MakeTable(len:%d)\n", BufferLength);
#endif
	double buffer[BufferLength];
	uint32_t x1,x2;
	double factor = xhyp.MakeTable(buffer, BufferLength, &x1, &x2, dOmegaTablePrecision); // make table of probability values
	
	uint32_t x;
	double sum;
	for (x = x1, sum = 0; x <= x2; x++) sum = buffer[x-x1] += sum; // Make table cumulative
	
	// NOTE: OPTIMIZED out as true==lower_tail  ==> this inversion does not need to happen
	//if (!lower_tail) p = 1. - p; // Invert if right tail
	p *= factor; // Table is scaled by factor
	
	// Binary search in table:
	uint32_t a, b, c; // Used in binary search
	a = 0; b = x2 - x1 + 1;
	while (a < b) {
		c = (a + b) / 2;
		if (p <= buffer[c]) {
			b = c;
		}
		else {
			a = c + 1;
		}
	}
	x = x1 + a;
	if (x > x2) x = x2;           // Prevent values > xmax that occur because of small imprecisions
#ifdef DUMP_CHR_NCHG
	if (dump)
		fprintf(stderr, "getQuantileNCHG(p:%g, ni:%d, nj:%d, n:%d, odds:%g, precision:%g, lower_tail:%s) BufferLength:%d factor:%g returns %d\n", p, ni, nj, n, odds, precision, lower_tail?"true":"false", BufferLength, factor, x);
#endif
	
	return x;
}


//------------------------------------------------------------------------------


CCCMatrixInteraction::CCCMatrixInteraction()
:isIntra(true),defaultVal(Interaction()),m_TotalInteractions(0),nThreads(1)
,orderedSegmentMidPoints(NULL)
,cacheFlag(0),cachedN(0),cachedExpectN(0.0), fDepthThreshold(0.25)
{
	//N=0;
}

CCCMatrixInteraction::~CCCMatrixInteraction()
{
	if (orderedSegmentMidPoints) {
		free(orderedSegmentMidPoints);
		orderedSegmentMidPoints = NULL;
	}
}

void CCCMatrixInteraction::LoadInteractions(set<Segment> &rows, set<Segment> &cols, Interaction defaultvalue, bool isintra, string& chr)
{
	isIntra = isintra;
	defaultVal = defaultvalue;
	m_chr = chr;
	
	// Sanity checks:
	assert(rows.size()>0);
	assert(cols.size()>0);
	if (isintra) {
		assert(rows.size()==cols.size());
	}
	else {
		// Obs: This could be made more sophisticated.
		// None of the elements in rows and cols should be allowed to be the same!
		assert(rows.size()!=cols.size());
	}
	
	assert(isintra); // WCH: we will have to handle trans-interaction differently (future work)
	// create the look up
	SegmentKey segmentCounter=0;
	for (set<Segment>::iterator it = rows.begin(); it != rows.end(); ++it) {
		segmentToKey[*it] = segmentCounter;
		segmentCounter++;
	}
	// store as row-major
	mat_skToIndex.resize(segmentCounter);
	mat_sk.resize(segmentCounter);
	mat_count.resize(segmentCounter);
	mat_pvalue.resize(segmentCounter);
	mat_mask.resize(segmentCounter);
}

void CCCMatrixInteraction::InteractionLoaded(uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/)
{
	// all the contacts have been loaded
	// we will switch over to crunching of data only
	if (0==orderedSegments.size()) {
		orderedSegments.resize(segmentToKey.size());
		orderedSegmentMidPoints = (uint32_t *) calloc(segmentToKey.size(), sizeof(uint32_t));
		if (!orderedSegmentMidPoints) {
			fprintf(stderr, "[E::%s:%d] Fail to allocate %lu unit of midpoint.\n", __func__, __LINE__, segmentToKey.size());
			exit(-1);
		}
		for (map<Segment, SegmentKey >::iterator itRow = segmentToKey.begin(); itRow != segmentToKey.end(); ++itRow) {
			orderedSegments[itRow->second] = SegmentMin(itRow->first);
			orderedSegmentMidPoints[itRow->second] = itRow->first.pos;
		}
		// we no longer needs the keys
		segmentToKey.clear();

		// need to reorder the vector data accordingly based on the map of segmentKey
		for(SegmentKey i=0; i<orderedSegments.size(); ++i) {
			//std::vector<std::map<SegmentKey, uint32_t > > mat_skToIndex;
			CountRowType counts; counts.resize(mat_count[i].size());
			PvalueRowType pvalues; pvalues.resize(mat_pvalue[i].size());
			MaskRowType masks; masks.resize(mat_mask[i].size());
			mat_sk[i].resize(mat_skToIndex[i].size());
			uint32_t j = 0;
			for(SegmentKeyMapType::iterator it=mat_skToIndex[i].begin(); it!=mat_skToIndex[i].end(); ++it, ++j) {
				counts[j] = mat_count[i][it->second];
				pvalues[j] = mat_pvalue[i][it->second];
				masks[j] = mat_mask[i][it->second];
				// TODO: decide if we need map or vector
				it->second = j;
				mat_sk[i][j] = it->first;
			}
			mat_count[i] = counts;
			mat_pvalue[i] = pvalues;
			mat_mask[i] = masks;
			
			// NOTE: we no longer need segmentKey-to-index
			//       all processing algo should traverse the list in a pass rather than binary-search
			mat_skToIndex[i].clear();
		}
		// NOTE: we no longer need segmentKey-to-index
		//       all processing algo should traverse the list in a pass rather than binary-search
		mat_skToIndex.clear();

		
		// let's remember masked count
		rowMasked.resize(orderedSegments.size(), 0);

		// let's cache the start and end position
		segmentKeyStarts.resize(orderedSegments.size(), (SegmentKey)(orderedSegments.size()));
		SegmentKey lastStartPos = 1;
		for(SegmentKey i=0; i<orderedSegments.size(); ++i) {
			uint32_t startPos = orderedSegmentMidPoints[i] + minDelta;
			for(SegmentKey j=lastStartPos; j<orderedSegments.size(); ++j) {
				if (orderedSegmentMidPoints[j]>=startPos) {
					segmentKeyStarts[i] = j;
					lastStartPos = j;
					break;
				}
			}
		}
		
		// let's cache the end position
		segmentKeyEnds.resize(orderedSegments.size(), (SegmentKey)(orderedSegments.size()));
		if (maxDelta<MAXDELTA) {
			// necessary to figure out where we are stopping
			// Let's work on the first row and propagate downward
			SegmentKey endStartKey = (int)orderedSegments.size() - 1;
			uint32_t endPos = orderedSegmentMidPoints[0] + maxDelta;
			for(SegmentKey j=endStartKey; j>0; --j) {
				if (orderedSegmentMidPoints[j]<=endPos) {
					segmentKeyStarts[0] = j+1;
					endStartKey = j;
					break;
				}
			}
			
			if (segmentKeyStarts[0]<orderedSegments.size()) {
				for(SegmentKey i=1; i<orderedSegments.size(); ++i) {
					uint32_t endPos = orderedSegmentMidPoints[i] + maxDelta;
					if (endPos<orderedSegmentMidPoints[i]) {
						// overflow occur!! Just include every points
						break;
					} else {
						for(SegmentKey j=endStartKey; j<orderedSegments.size(); ++j) {
							if (orderedSegmentMidPoints[j]>endPos) {
								segmentKeyStarts[i] = j;
								endStartKey = j;
								break;
							}
						}
					}
				}
			}
		}
		
		// let's compute the prefix segment start (i.e. the negative delta from lower triangle)
		// let's compute the prefix segment ends (i.e. the negative delta from lower triangle)
		segmentPrefixKeyStarts.resize(orderedSegments.size(), (SegmentKey)(orderedSegments.size()));
		segmentPrefixKeyEnds.resize(orderedSegments.size(), (SegmentKey)(orderedSegments.size()));
		for(SegmentKey i=1; i<orderedSegments.size(); ++i) {
			for(SegmentKey j=0; j<i; ++j) {
				if (segmentKeyStarts[j]<=i) {
					segmentPrefixKeyStarts[i] = j;
					break;
				}
			}
			if (segmentPrefixKeyStarts[i]!=orderedSegments.size()) {
				for(SegmentKey j=i-1; j>=segmentPrefixKeyStarts[i]; --j) {
					if (segmentKeyStarts[j]<=i) {
						segmentPrefixKeyEnds[i] = j+1;
						break;
					}
				}
			}
		}
		
#if 0
		// TODO: DEBUG
		for(SegmentKey i=0; i<orderedSegments.size(); ++i) {
			cerr << "i:" << i << " start:" << segmentKeyStarts[i] << " end:" << segmentKeyEnds[i] << " Pstart:" << segmentPrefixKeyStarts[i] << " Pend:" << segmentPrefixKeyEnds[i] << endl;
		}
#endif
		
		// cache the start & end iterators
		//std::vector<std::map<SegmentKey, Interaction >::iterator > segmentIterStarts;
		//std::vector<std::map<SegmentKey, Interaction >::iterator > segmentIterEnds;
		segmentIterStarts.resize(orderedSegments.size());
		segmentIterEnds.resize(orderedSegments.size());
		for(SegmentKey i=0; i<segmentKeyStarts.size(); ++i) {
			segmentIterStarts[i] = segmentIterEnds[i] = (uint32_t)mat_sk[i].size();
			if (segmentKeyStarts[i]<orderedSegments.size()) {
				// find iterator for start; exact key might not be part of interaction
				for(uint32_t j=0; j<mat_sk[i].size(); ++j) {
					if (segmentKeyStarts[i]<=mat_sk[i][j]) {
						segmentIterStarts[i] = j;
						break;
					}
				}
				if (segmentKeyEnds[i]<orderedSegments.size()) {
					// find iterator for end
					for(uint32_t j=segmentIterStarts[i]; j<mat_sk[i].size(); ++j) {
						if (segmentKeyEnds[i]>mat_sk[i][j]) {
							segmentIterEnds[i] = j;
							break;
						}
					}
				}
			}
		}
	}
}

// NO optimization needed
void CCCMatrixInteraction::printInteractionMatrix(ostream& os, bool hideMask /*=false*/) {
	// print column header
	os << "Count\t.\tColSum";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << m_getMarginal[segmentCounter];
	os << endl;
	
	os << ".\t.\tIndex";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << segmentCounter;
	os << endl;
	
	os << "RowSum\tIndex\tmidpt";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << orderedSegmentMidPoints[segmentCounter];
	os << endl;
	
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) {
		// prepare a row
		CountRowType aCountRow; aCountRow.resize(orderedSegments.size(), defaultVal.count);
		MaskRowType aMaskRow; aMaskRow.resize(orderedSegments.size(), (0!=defaultVal.mask));
		
		{
			// populate known columns
			SegmentKeyRowType& rowSegmentKeys = mat_sk[segmentCounter];
			CountRowType& currCountRow = mat_count[segmentCounter];
			MaskRowType& currMaskRow = mat_mask[segmentCounter];
			for(uint32_t i=0; i<rowSegmentKeys.size(); ++i) {
				aCountRow[rowSegmentKeys[i]] = currCountRow[i];
				aMaskRow[rowSegmentKeys[i]] = currMaskRow[i];
			}
		}

		// WCH: decided that we should only print UT so that we don't need to double the index structure size!!
		// print row header
		os << m_getMarginal[segmentCounter] << "\t" << segmentCounter << "\t" << orderedSegmentMidPoints[segmentCounter];
		// report
		for(SegmentKey col=0; col<=segmentCounter; ++col) os << "\t.";
		for(SegmentKey col=segmentCounter+1; col<orderedSegments.size(); ++col)
			if (aMaskRow[col] && hideMask) os << "\tNA";
			else os  << "\t" << aCountRow[col];
		os << endl;
	}
}

// NO optimization needed
void CCCMatrixInteraction::printPvalueMatrix(ostream& os, bool hideMask /*=false*/) {
	// print column header
	os << "Pvalue\t.\tColSum";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << m_getExp[segmentCounter];
	os << endl;
	
	os << ".\t.\tIndex";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << segmentCounter;
	os << endl;
	
	os << "RowSum\tIndex\tmidpt";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << orderedSegmentMidPoints[segmentCounter];
	os << endl;
	
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) {
		// prepare a row
		PvalueRowType aPvalueRow; aPvalueRow.resize(orderedSegments.size(), defaultVal.pvalue);
		MaskRowType aMaskRow; aMaskRow.resize(orderedSegments.size(), (0!=defaultVal.mask));
		
		{
			// populate known columns
			SegmentKeyRowType& rowSegmentKeys = mat_sk[segmentCounter];
			PvalueRowType& currPvalueRow = mat_pvalue[segmentCounter];
			MaskRowType& currMaskRow = mat_mask[segmentCounter];
			for(uint32_t i=0; i<rowSegmentKeys.size(); ++i) {
				aPvalueRow[rowSegmentKeys[i]] = currPvalueRow[i];
				aMaskRow[rowSegmentKeys[i]] = currMaskRow[i];
			}
		}
		
		// WCH: decided that we should only print UT so that we don't need to double the index structure size!!
		// print row header
		os << m_getExp[segmentCounter] << "\t" << segmentCounter << "\t" << orderedSegmentMidPoints[segmentCounter];
		// report
		for(SegmentKey col=0; col<=segmentCounter; ++col) os << "\t.";
		for(SegmentKey col=segmentCounter+1; col<orderedSegments.size(); ++col) {
			if (aMaskRow[col] && hideMask) os << "\tNA";
			else os  << "\t" << aPvalueRow[col];
		}
		os << endl;
	}
}

// NO optimization needed
// WCH: this is the only version to set the element in a matrix by "Segment" rather than "SegmentKey"
void CCCMatrixInteraction::setElement(Segment& row, Segment& col, uint32_t count) {
	//assert( this->isValidRow(row));
	//assert( this->isValidColumn(col));
	
	SegmentKey rowKey = segmentToKey.at(row);
	SegmentKey colKey = segmentToKey.at(col);
	
	assert(rowKey!=colKey);
	if (rowKey > colKey) {
		SegmentKey t = colKey; colKey = rowKey; rowKey = t;
	}
	
	m_TotalInteractions++;
	
	mat_skToIndex[rowKey][colKey] = (uint32_t) mat_count[rowKey].size();
	mat_count[rowKey].push_back(count);
	mat_pvalue[rowKey].push_back(defaultVal.pvalue);
	mat_mask[rowKey].push_back(0!=defaultVal.mask);
	
	if (isIntra) {
		m_TotalInteractions++;
		// WCH: we only store in UT format, thus no need to store the transposed element
	}
}

struct interaction_RowSum_worker_t {
	interaction_RowSum_worker_t(SegmentKeyMatrixType& ms, CountMatrixType& mc, SegmentIterType& sis, SegmentIterType& sie, CountRowType& rs)
	: mat_sk(ms), mat_count(mc), segmentIterStarts(sis), segmentIterEnds(sie), rowSums(rs)
	{}
	
	SegmentKeyMatrixType& mat_sk;
	CountMatrixType& mat_count;
	SegmentIterType& segmentIterStarts;
	SegmentIterType& segmentIterEnds;
	CountRowType& rowSums;
};


// [opt19.2: interaction_RowSum_worker {int}cell-dist,count-LT,UT]
static void interaction_RowSum_worker(void *data, int rowKey, int tid)
{
	// TODO: FIXME: we should arrange for no lock operation!!!
	
	interaction_RowSum_worker_t *w = (interaction_RowSum_worker_t*)data;
	uint32_t sum=0;
	SegmentKeyRowType& skRow = w->mat_sk[rowKey];
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		sum += countRow[i];
		__sync_fetch_and_add(&(w->rowSums[skRow[i]]), countRow[i]);
	}
	
	__sync_fetch_and_add(&(w->rowSums[rowKey]), sum);
}

void CCCMatrixInteraction::getRowSums(CountRowType &res) {
	res.clear();
	res.resize(orderedSegments.size());
	
	interaction_RowSum_worker_t w(mat_sk, mat_count, segmentIterStarts, segmentIterEnds, res);
	kt_for(nThreads, interaction_RowSum_worker, &w, (int)orderedSegments.size());
}


struct getDeltasCountSum_worker_t {
	getDeltasCountSum_worker_t(uint32_t* mps, SegmentKeyRowType& ss, SegmentKeyRowType& se, vector<uint32_t >& rm, SegmentKeyMatrixType& i2sk, CountMatrixType& mc, MaskMatrixType& mm, uint32_t* dc, SegmentIterType& sis, SegmentIterType& sie, uint32_t* ds, BIT_ARRAY *p, uint32_t ms, vector<uint64_t >& ns)
	: midpoints(mps), segmentKeyStarts(ss), segmentKeyEnds(se), rowMasked(rm), mat_sk(i2sk), mat_count(mc), mat_mask(mm), deltaCounts(dc), segmentIterStarts(sis), segmentIterEnds(sie), deltaSums(ds), presences(p), maxSpan(ms), totalDeltas(ns)
	{}
	
	vector<uint64_t >& totalDeltas;
	uint32_t* deltaCounts;
	uint32_t* deltaSums;
	BIT_ARRAY *presences;
	uint32_t maxSpan;
	
	uint32_t* midpoints;
	SegmentKeyRowType& segmentKeyStarts;
	SegmentKeyRowType& segmentKeyEnds;
	SegmentKeyMatrixType& mat_sk;
	CountMatrixType& mat_count;
	MaskMatrixType& mat_mask;
	
	SegmentIterType& segmentIterStarts;
	SegmentIterType& segmentIterEnds;
	
	vector<uint32_t >& rowMasked;
};

#if 0
//[opt16]
struct getDeltasCount_worker_t {
	getDeltasCount_worker_t(uint32_t* mps, SegmentKeyRowType& ss, SegmentKeyRowType& se, vector<uint32_t >& rm, SegmentKeyMatrixType& i2sk, MaskMatrixType& mm, uint32_t* r, BIT_ARRAY *p, uint32_t ms, vector<uint64_t >& ns)
	: midpoints(mps), segmentKeyStarts(ss), segmentKeyEnds(se), rowMasked(rm), mat_sk(i2sk), mat_mask(mm), res(r), presences(p), maxSpan(ms), totalDeltas(ns)
	{}
	
	vector<uint64_t >& totalDeltas;
	uint32_t* res;
	BIT_ARRAY *presences;
	uint32_t maxSpan;
	
	uint32_t* midpoints;
	SegmentKeyRowType& segmentKeyStarts;
	SegmentKeyRowType& segmentKeyEnds;
	SegmentKeyMatrixType& mat_sk;
	MaskMatrixType& mat_mask;
	
	vector<uint32_t >& rowMasked;
};
#endif

//[opt16]
// [opt19.2: getDeltasCount_worker {perm}cell-dist-UT]
static void getDeltasCountLoad_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		__sync_fetch_and_add((w->deltaCounts+delta), 1);
	}
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey];
	w->totalDeltas[tid] += totalDeltas;
}

static void getDeltasCountLoadCollate_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		__sync_fetch_and_add((w->deltaCounts+delta), 1);
		bitset_set_mt(w->presences->words, delta);
	}
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey];
	w->totalDeltas[tid] += totalDeltas;
}

// [opt19.2: getDeltasCount_worker {perm}cell-dist-UT]
static void getDeltasCount_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey];
	w->totalDeltas[tid] += totalDeltas;
}

// [opt19.2: getMaskedDeltasCountLoad_worker {int,perm}cell-dist,mask-UT]
static void getMaskedDeltasCountLoad_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		// we assume that it is presence and undo laster if it is masked
		__sync_fetch_and_add((w->deltaCounts+delta), 1);
	}
	SegmentKeyRowType& rowSK = w->mat_sk[rowKey];
	MaskRowType& rowMask = w->mat_mask[rowKey];
	for(uint32_t i=0; i<rowMask.size(); ++i) {
		if (rowMask[i]) {
			uint32_t delta=w->midpoints[rowSK[i]]-rowMidPoint;
			__sync_fetch_and_sub((w->deltaCounts+delta), 1);
			// NOTE: masked out cell MUST be marked for precomputation
			bitset_set_mt(w->presences->words, delta);
		}
	}
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey] - w->rowMasked[rowKey];
	w->totalDeltas[tid] += totalDeltas;
}

static void getMaskedDeltasCountLoadCollate_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		// we assume that it is presence and undo laster if it is masked
		__sync_fetch_and_add((w->deltaCounts+delta), 1);
		bitset_set_mt(w->presences->words, delta);
	}
	SegmentKeyRowType& rowSK = w->mat_sk[rowKey];
	MaskRowType& rowMask = w->mat_mask[rowKey];
	for(uint32_t i=0; i<rowMask.size(); ++i) {
		if (rowMask[i]) {
			uint32_t delta=w->midpoints[rowSK[i]]-rowMidPoint;
			__sync_fetch_and_sub((w->deltaCounts+delta), 1);
			// NOTE: all cells are already marked in the loop above; i.e. redundant here
			//bitset_set_mt(w->presences->words, delta);
		}
	}
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey] - w->rowMasked[rowKey];
	w->totalDeltas[tid] += totalDeltas;
}

// [opt19.2: getMaskedDeltasCount_worker {int,perm}cell-dist,mask-UT]
static void getMaskedDeltasCount_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	unsigned long totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey] - w->rowMasked[rowKey];
	w->totalDeltas[tid] += totalDeltas;
}

static void collate_Deltas_indicator_worker(void *data, int wordId, int tid)
{
	// we will check if there is any delta counter set
	// if so, we propagate it to the presence bit array
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	//#define bitset64_wrd(pos) ((pos) >> 6)
	//#define bitset64_idx(pos) ((pos) & 63)
	int rowKey = wordId * 64; //* bit_per_word
	int endIdx = 64;
	if (endIdx>(w->maxSpan-rowKey)) endIdx = w->maxSpan-rowKey;
	for(int i=0; i<endIdx; ++rowKey, ++i)
		if (w->deltaCounts[rowKey]>0) {
			//bitset_set(w->presences->words, rowKey); // NO multi-threading protection as each thread has their own words!
			bitset2_set(w->presences->words, wordId, i);
		}
}

//[opt16]
uint64_t CCCMatrixInteraction::getDeltasCount(uint32_t* res, BIT_ARRAY *presences, uint32_t maxSpan, uint32_t minDelta /* = 0 */, uint32_t maxDelta /*=MAXDELTA*/, bool masking /* =false */)
{
	vector<uint64_t > tmp; tmp.resize(nThreads);
	getDeltasCountSum_worker_t w(orderedSegmentMidPoints, segmentKeyStarts, segmentKeyEnds, rowMasked, mat_sk, mat_count, mat_mask, res, segmentIterStarts, segmentIterEnds, NULL, presences, maxSpan, tmp);
	if (masking) {
		if (res) {
			float estimateDepth = orderedSegments.size() * (orderedSegments.size()-1);
			estimateDepth /= maxSpan;
			if (estimateDepth>=fDepthThreshold) {
				kt_for(nThreads, getMaskedDeltasCountLoad_worker, &w, (int)orderedSegments.size());
				kt_for(nThreads, collate_Deltas_indicator_worker, &w, (int)presences->num_of_words);
			} else {
				kt_for(nThreads, getMaskedDeltasCountLoadCollate_worker, &w, (int)orderedSegments.size());
			}
		}
		else kt_for(nThreads, getMaskedDeltasCount_worker, &w, (int)orderedSegments.size());
	} else {
		if (res) {
			float estimateDepth = orderedSegments.size() * (orderedSegments.size()-1);
			estimateDepth /= maxSpan;
			if (estimateDepth>=fDepthThreshold) {
				kt_for(nThreads, getDeltasCountLoad_worker, &w, (int)orderedSegments.size());
				kt_for(nThreads, collate_Deltas_indicator_worker, &w, (int)presences->num_of_words);
			} else {
				kt_for(nThreads, getDeltasCountLoadCollate_worker, &w, (int)orderedSegments.size());
			}
		}
		else kt_for(nThreads, getDeltasCount_worker, &w, (int)orderedSegments.size());
	}
	
	uint64_t totalDeltas = accumulate(tmp.begin(), tmp.end(), (uint64_t)0);
	return totalDeltas;
}

// TODO:
static void getDeltasCountSumLoad_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		__sync_fetch_and_add((w->deltaCounts+delta), 1);
	}
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey];
	w->totalDeltas[tid] += totalDeltas;
	
	SegmentKeyRowType& skRow = w->mat_sk[rowKey];
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		uint32_t delta=w->midpoints[skRow[i]]-rowMidPoint;
		__sync_fetch_and_add((w->deltaSums+delta), countRow[i]);
	}
}

// TODO:
static void getDeltasCountSumLoadCollate_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		__sync_fetch_and_add((w->deltaCounts+delta), 1);
		bitset_set_mt(w->presences->words, delta);
	}
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey];
	w->totalDeltas[tid] += totalDeltas;
	
	SegmentKeyRowType& skRow = w->mat_sk[rowKey];
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		uint32_t delta=w->midpoints[skRow[i]]-rowMidPoint;
		__sync_fetch_and_add((w->deltaSums+delta), countRow[i]);
	}
}

// TODO:
static void getMaskedDeltasCountSumLoad_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		// we assume that it is presence and undo laster if it is masked
		__sync_fetch_and_add((w->deltaCounts+delta), 1);
	}
	SegmentKeyRowType& rowSK = w->mat_sk[rowKey];
	MaskRowType& rowMask = w->mat_mask[rowKey];
	for(uint32_t i=0; i<rowMask.size(); ++i) {
		if (rowMask[i]) {
			uint32_t delta=w->midpoints[rowSK[i]]-rowMidPoint;
			__sync_fetch_and_sub((w->deltaCounts+delta), 1);
			// NOTE: masked out cell MUST be marked for precomputation
			bitset_set_mt(w->presences->words, delta);
		}
	}
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey] - w->rowMasked[rowKey];
	w->totalDeltas[tid] += totalDeltas;
	
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		uint32_t delta=w->midpoints[rowSK[i]]-rowMidPoint;
		if (!rowMask[i])
			__sync_fetch_and_add((w->deltaSums+delta), countRow[i]);
	}
}

// TODO:
static void getMaskedDeltasCountSumLoadCollate_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		// we assume that it is presence and undo laster if it is masked
		__sync_fetch_and_add((w->deltaCounts+delta), 1);
		bitset_set_mt(w->presences->words, delta);
	}
	SegmentKeyRowType& rowSK = w->mat_sk[rowKey];
	MaskRowType& rowMask = w->mat_mask[rowKey];
	for(uint32_t i=0; i<rowMask.size(); ++i) {
		if (rowMask[i]) {
			uint32_t delta=w->midpoints[rowSK[i]]-rowMidPoint;
			__sync_fetch_and_sub((w->deltaCounts+delta), 1);
			// NOTE: all cells are already marked in the loop above; i.e. redundant here
			//bitset_set_mt(w->presences->words, delta);
		}
	}
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey] - w->rowMasked[rowKey];
	w->totalDeltas[tid] += totalDeltas;
	
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		uint32_t delta=w->midpoints[rowSK[i]]-rowMidPoint;
		if (!rowMask[i])
			__sync_fetch_and_add((w->deltaSums+delta), countRow[i]);
	}
}

uint64_t CCCMatrixInteraction::getDeltasCountSum(uint32_t* deltaCounts, uint32_t* deltaSums, BIT_ARRAY *presences, uint32_t maxSpan, uint32_t minDelta /* = 0 */, uint32_t maxDelta /*=MAXDELTA*/, bool masking /* =false */)
{
	vector<uint64_t > tmp; tmp.resize(nThreads);
	getDeltasCountSum_worker_t w(orderedSegmentMidPoints, segmentKeyStarts, segmentKeyEnds, rowMasked, mat_sk, mat_count, mat_mask, deltaCounts, segmentIterStarts, segmentIterEnds, deltaSums, presences, maxSpan, tmp);
	if (masking) {
		if (deltaCounts && deltaSums) {
			double estimateDepth = orderedSegments.size() * (orderedSegments.size()-1);
			estimateDepth /= maxSpan;
			if (estimateDepth>=fDepthThreshold) {
				kt_for(nThreads, getMaskedDeltasCountSumLoad_worker, &w, (int)orderedSegments.size());
				kt_for(nThreads, collate_Deltas_indicator_worker, &w, (int)presences->num_of_words);
			} else {
				kt_for(nThreads, getMaskedDeltasCountSumLoadCollate_worker, &w, (int)orderedSegments.size());
			}
		}
		else kt_for(nThreads, getMaskedDeltasCount_worker, &w, (int)orderedSegments.size());
	} else {
		if (deltaCounts && deltaSums) {
			double estimateDepth = orderedSegments.size() * (orderedSegments.size()-1);
			estimateDepth /= maxSpan;
			if (estimateDepth>=fDepthThreshold) {
				kt_for(nThreads, getDeltasCountSumLoad_worker, &w, (int)orderedSegments.size());
				kt_for(nThreads, collate_Deltas_indicator_worker, &w, (int)presences->num_of_words);
			} else {
				kt_for(nThreads, getDeltasCountSumLoadCollate_worker, &w, (int)orderedSegments.size());
			}
		}
		else kt_for(nThreads, getDeltasCount_worker, &w, (int)orderedSegments.size());
	}
	
	uint64_t totalDeltas = accumulate(tmp.begin(), tmp.end(), (uint64_t)0);
	return totalDeltas;
}

// TODO:
static void maskDeltasCountSumLoadCollate_worker(void *data, int rowKey, int tid)
{
	getDeltasCountSum_worker_t *w = (getDeltasCountSum_worker_t*)data;

	uint32_t rowMidPoint = w->midpoints[rowKey];
	SegmentKeyRowType& rowSK = w->mat_sk[rowKey];
	MaskRowType& rowMask = w->mat_mask[rowKey];
	for(uint32_t i=0; i<rowMask.size(); ++i) {
		if (rowMask[i]) {
			uint32_t delta=w->midpoints[rowSK[i]]-rowMidPoint;
			__sync_fetch_and_sub((w->deltaCounts+delta), 1);
			// NOTE: all cells are already marked in the loop above; i.e. redundant here
			//bitset_set_mt(w->presences->words, delta);
		}
	}
	
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		uint32_t delta=w->midpoints[rowSK[i]]-rowMidPoint;
		if (rowMask[i])
			__sync_fetch_and_sub((w->deltaSums+delta), countRow[i]);
	}
	
	uint32_t totalDeltas = w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey] - w->rowMasked[rowKey];
	w->totalDeltas[tid] += totalDeltas;
}

uint64_t CCCMatrixInteraction::maskDeltasCountSum(uint32_t* deltaCounts, uint32_t* deltaSums, BIT_ARRAY *presences, uint32_t maxSpan, uint32_t minDelta /* = 0 */, uint32_t maxDelta /*=MAXDELTA*/)
{
	vector<uint64_t > tmp; tmp.resize(nThreads);
	getDeltasCountSum_worker_t w(orderedSegmentMidPoints, segmentKeyStarts, segmentKeyEnds, rowMasked, mat_sk, mat_count, mat_mask, deltaCounts, segmentIterStarts, segmentIterEnds, deltaSums, presences, maxSpan, tmp);
	if (deltaCounts && deltaSums) {
		kt_for(nThreads, maskDeltasCountSumLoadCollate_worker, &w, (int)orderedSegments.size());
	}
	else kt_for(nThreads, getMaskedDeltasCount_worker, &w, (int)orderedSegments.size());
	
	uint64_t totalMaskedDeltas = accumulate(tmp.begin(), tmp.end(), (uint64_t)0);
	return totalMaskedDeltas;
}

struct getN_worker_t {
	getN_worker_t(CountMatrixType& mc, SegmentIterType& sis, SegmentIterType& sie, CountRowType& ts)
	: mat_count(mc), segmentIterStarts(sis), segmentIterEnds(sie), threadSums(ts)
	{}
	
	CountMatrixType& mat_count;
	SegmentIterType& segmentIterStarts;
	SegmentIterType& segmentIterEnds;
	CountRowType& threadSums;
};

// [opt19.2: getN_worker {int}cell-dist,count-UT]
static void getN_worker(void *data, int rowKey, int tid)
{
	getN_worker_t *w = (getN_worker_t*)data;
	uint32_t sum = 0;
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		sum += countRow[i];
	}
	w->threadSums[tid] += sum;
}

// [opt19.1]
void CCCMatrixInteraction::calculate_cache_N(uint32_t minDelta/*=0*/, uint32_t maxDelta/*=MAXDELTA*/) {
	if (CACHE_FLAG_N!=(CACHE_FLAG_N&cacheFlag)) {
		// TODO: WCH: we only affect count rather than p-value here
		// TODO: WCH: check that we did not scale the same reference twice if we are storing both upper and lower triangle!!!
		vector<uint32_t > res; res.resize(nThreads);
		getN_worker_t w(mat_count, segmentIterStarts, segmentIterEnds, res);
		kt_for(nThreads, getN_worker, &w, (int)orderedSegments.size());
		cachedN = accumulate(res.begin(), res.end(), (uint32_t)0);
		
		cacheFlag|=CACHE_FLAG_N;
	}
	
	//cerr << "getN: " << cachedN << endl;
}

struct getStatisticsPerDelta_worker_t {
	getStatisticsPerDelta_worker_t(uint32_t* mps, SegmentKeyRowType& ss, SegmentKeyRowType& se, SegmentKeyMatrixType& ms, CountMatrixType& mc, MaskMatrixType& mm, SegmentIterType& sis, SegmentIterType& sie, uint32_t* pqtm, vector<vector<DeltaStatistics > >& r)
	: midpoints(mps), segmentKeyStarts(ss), segmentKeyEnds(se), mat_sk(ms), mat_count(mc), mat_mask(mm), segmentIterStarts(sis), segmentIterEnds(sie), pQuantileMapper(pqtm), res(r)
	{}
	
	uint32_t* midpoints;
	SegmentKeyRowType& segmentKeyStarts;
	SegmentKeyRowType& segmentKeyEnds;
	SegmentKeyMatrixType& mat_sk;
	CountMatrixType& mat_count;
	MaskMatrixType& mat_mask;
	SegmentIterType& segmentIterStarts;
	SegmentIterType& segmentIterEnds;
	uint32_t* pQuantileMapper;
	vector<vector<DeltaStatistics > >& res;
};

// [opt19.2: getN_worker {int}cell-dist,count-UT]
static void getStatisticsPerDeltaNonMasked_worker(void *data, int rowKey, int tid)
{
	getStatisticsPerDelta_worker_t *w = (getStatisticsPerDelta_worker_t*)data;
	vector<DeltaStatistics >& tres = w->res[tid];
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for (SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		tres[w->pQuantileMapper[delta]].count++;
	}
	
	SegmentKeyRowType& skRow = w->mat_sk[rowKey];
	CountRowType& countRow = w->mat_count[rowKey];
	MaskRowType& maskRow = w->mat_mask[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		uint32_t delta=w->midpoints[skRow[i]]-rowMidPoint;
		if (maskRow[i]) {
			// we counted this entry in the 1st loop, so have to remove it
			tres[w->pQuantileMapper[delta]].count--;
		} else {
			tres[w->pQuantileMapper[delta]].sum+=countRow[i];
		}
	}
}

static void getStatisticsPerDelta_worker(void *data, int rowKey, int tid)
{
	getStatisticsPerDelta_worker_t *w = (getStatisticsPerDelta_worker_t*)data;
	vector<DeltaStatistics >& tres = w->res[tid];
	
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for (SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		tres[w->pQuantileMapper[delta]].count++;
	}

	SegmentKeyRowType& skRow = w->mat_sk[rowKey];
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		uint32_t delta=w->midpoints[skRow[i]]-rowMidPoint;
		tres[w->pQuantileMapper[delta]].sum+=countRow[i];
	}
}

// TODO: OPTIMIZE : multithreading
// [opt19.2: getStatisticsPerDelta {int,perm}cell-dist,mask,count-UT]
void CCCMatrixInteraction::getStatisticsPerDelta(vector<DeltaStatistics >& res, uint32_t* pQuantileMapper, uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/, bool masking/*=false*/)
{
	// WCH: original function include those which have count = 0
	//      we have to keep it the same so that the spline data points will be the same!!!
	// WCH: So, we cannot have a more efficient version which only compute what is necessary, i.e. count > 0
	vector<vector<DeltaStatistics > > tmp; tmp.resize(nThreads);
	for(uint32_t i=0; i<nThreads; ++i) tmp[i].resize(res.size());
	
	getStatisticsPerDelta_worker_t w(orderedSegmentMidPoints, segmentKeyStarts, segmentKeyEnds, mat_sk, mat_count, mat_mask, segmentIterStarts, segmentIterEnds, pQuantileMapper, tmp);
	if (masking) {
		kt_for(nThreads, getStatisticsPerDeltaNonMasked_worker, &w, (int)orderedSegments.size());
	} else {
		kt_for(nThreads, getStatisticsPerDelta_worker, &w, (int)orderedSegments.size());
	}
	
	for(uint32_t j=0; j<res.size(); ++j) {
		for(uint32_t i=0; i<nThreads; ++i) {
			res[j].sum += tmp[i][j].sum;
			res[j].count += tmp[i][j].count;
		}
	}
}

struct interaction_Positives_worker_t {
	interaction_Positives_worker_t(CountMatrixType& mc, PvalueMatrixType& mp, vector<uint64_t >& p)
	: mat_count(mc), mat_pvalue(mp), positives(p)
	{}
	
	double alpha;
	uint32_t cutoff;
	
	CountMatrixType& mat_count;
	PvalueMatrixType& mat_pvalue;
	vector<uint64_t >& positives;
};

// [opt19.2: interaction_Positives_A_worker {int}cell-pvalue-UT]
static void interaction_Positives_A_worker(void *data, int rowKey, int tid)
{
	interaction_Positives_worker_t *w = (interaction_Positives_worker_t*)data;
	uint64_t count = 0;
	PvalueRowType& pvalueRow = w->mat_pvalue[rowKey];
	for(uint32_t i=0; i<pvalueRow.size(); ++i)
		if (pvalueRow[i]<=w->alpha) count++;
	w->positives[tid]+=count;
}

// [opt19.2: interaction_Positives_AC_worker {int}cell-pvalue,count-UT]
static void interaction_Positives_AC_worker(void *data, int rowKey, int tid)
{
	interaction_Positives_worker_t *w = (interaction_Positives_worker_t*)data;
	uint64_t count = 0;
	PvalueRowType& pvalueRow = w->mat_pvalue[rowKey];
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=0; i<countRow.size(); ++i)
		if (countRow[i]>=w->cutoff and pvalueRow[i]<=w->alpha) count++;
	w->positives[tid]+=count;
}

// NO caller
uint64_t CCCMatrixInteraction::getNumberOfPositives(double alpha)
{
	assert(isIntra);
	if (0==alpha) return 0;
	
	vector<uint64_t> positives; positives.resize(nThreads);
	interaction_Positives_worker_t w(mat_count, mat_pvalue, positives);
	w.alpha = alpha;
	kt_for(nThreads, interaction_Positives_A_worker, &w, (int)orderedSegments.size());
	uint64_t res = accumulate(positives.begin(), positives.end(), (uint64_t)0);
	return res;
}

uint64_t CCCMatrixInteraction::getNumberOfPositives(double alpha, uint32_t cutoff)
{
	assert(isIntra);
	if (0==alpha) return 0;
	
	vector<uint64_t> positives; positives.resize(nThreads);
	interaction_Positives_worker_t w(mat_count, mat_pvalue, positives);
	w.alpha = alpha; w.cutoff = cutoff;
	kt_for(nThreads, interaction_Positives_AC_worker, &w, (int)orderedSegments.size());
	uint64_t res = accumulate(positives.begin(), positives.end(), (uint64_t)0);
	return res;
}

struct getPositivesBound_worker_t {
	getPositivesBound_worker_t(uint32_t c, CountMatrixType& mc, vector<uint64_t >& p)
	: cutoff(c), mat_count(mc), positives(p)
	{}
	
	uint32_t cutoff;
	CountMatrixType& mat_count;
	vector<uint64_t >& positives;
};

// [opt19.2: getAlphaBound_worker {int}cell-pvalue-UT]
static void getPositivesBound_worker(void *data, int rowKey, int tid)
{
	getPositivesBound_worker_t *w = (getPositivesBound_worker_t*)data;
	uint64_t positives = 0;
	CountRowType& countRow = w->mat_count[rowKey];
	for(uint32_t i=0; i<countRow.size(); ++i)
		if (countRow[i]>=w->cutoff) positives++;
	w->positives[tid] += positives;
}

void CCCMatrixInteraction::getPositivesBound(uint32_t cutoff, uint64_t& cutoffMaxPositive)
{
	vector<uint64_t> positives; positives.resize(nThreads);
	getPositivesBound_worker_t w(cutoff, mat_count, positives);
	kt_for(nThreads, getPositivesBound_worker, &w, (int)orderedSegments.size());
	
	cutoffMaxPositive = accumulate(positives.begin(), positives.end(), (uint64_t)cutoffMaxPositive);
}

struct maskByAlpha_worker_t {
	maskByAlpha_worker_t(CountMatrixType& mc, PvalueMatrixType& mp, MaskMatrixType& mm, vector<uint32_t >& rm, vector<uint64_t>& tm)
	: mat_count(mc), mat_pvalue(mp), mat_mask(mm), masked(rm), maskedThread(tm)
	{}
	
	double alpha;
	uint32_t cutoff;
	
	CountMatrixType& mat_count;
	PvalueMatrixType& mat_pvalue;
	MaskMatrixType& mat_mask;
	
	vector<uint32_t >& masked;
	vector<uint64_t >& maskedThread;
};

static void maskByAlpha_worker(void *data, int rowKey, int tid)
{
	maskByAlpha_worker_t *w = (maskByAlpha_worker_t*)data;
	
	uint32_t count = 0;
	CountRowType& countRow = w->mat_count[rowKey];
	PvalueRowType& pvalueRow = w->mat_pvalue[rowKey];
	MaskRowType& maskRow = w->mat_mask[rowKey];
	for(uint32_t i=0; i<countRow.size(); ++i)
		if (countRow[i]>=w->cutoff and pvalueRow[i]<=w->alpha) {
			maskRow[i] = true;
			count++;
		}
	w->masked[rowKey]=count;
	w->maskedThread[tid]+=count;
}

// OPTIMIZE: not a hotspot, skip optimization for now
// [opt19.2: maskByAlpha: {int}cell-pvalue,count,mask-UT]
uint64_t CCCMatrixInteraction::maskByAlpha(double alpha, uint32_t cutoff, bool dump/*=false*/)
{
	vector<uint64_t> tmpRes; tmpRes.resize(nThreads, 0);
	maskByAlpha_worker_t w(mat_count, mat_pvalue, mat_mask, rowMasked, tmpRes);
	w.alpha = alpha; w.cutoff=cutoff;
	kt_for(nThreads, maskByAlpha_worker, &w, (int)orderedSegments.size());
	uint64_t masked = accumulate(tmpRes.begin(), tmpRes.end(), (uint64_t)0);
	return masked;
}

struct calculate_Pvalues_worker_t {
	calculate_Pvalues_worker_t(uint32_t* mps, SegmentKeyMatrixType& ms, CountMatrixType& mc, PvalueMatrixType& mp, SegmentIterType& sis, SegmentIterType& sie, ExpectRowType& ge, MarginalRowType& gm, spline* sp)
	: midpoints(mps), mat_sk(ms), mat_count(mc), mat_pvalue(mp), segmentIterStarts(sis), segmentIterEnds(sie), getExp(ge), getMarginal(gm), delta2ExpSpline(sp), dump(false)
	{}
	
	uint32_t N2x;
	double L2x;
	
	uint32_t* midpoints;
	SegmentKeyMatrixType& mat_sk;
	CountMatrixType& mat_count;
	PvalueMatrixType& mat_pvalue;
	SegmentIterType& segmentIterStarts;
	SegmentIterType& segmentIterEnds;
	ExpectRowType& getExp;
	MarginalRowType& getMarginal;
	spline* delta2ExpSpline;
	
	bool dump;
};

// [opt19.2: calculate_Pvales_worker {int}cell-dist,pvalue-UT]
static void calculate_Pvalues_worker(void *data, int rowKey, int tid)
{
	calculate_Pvalues_worker_t *w = (calculate_Pvalues_worker_t*)data;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	SegmentKeyRowType& skRow = w->mat_sk[rowKey];
	CountRowType& countRow = w->mat_count[rowKey];
	PvalueRowType& pvalueRow = w->mat_pvalue[rowKey];
	for(uint32_t i=w->segmentIterStarts[rowKey]; i<w->segmentIterEnds[rowKey]; ++i) {
		SegmentKey colKey = skRow[i];
		uint32_t delta=w->midpoints[colKey] - rowMidPoint;
		//double Li = w->getExp[rowKey];
		//double Lj = w->getExp[colKey];
		//uint32_t Ni = w->getMarginal[rowKey];
		//uint32_t Nj = w->getMarginal[colKey];
		
		//uint32_t Nij = it->second.count;
		double Lij = max(w->delta2ExpSpline->cachedAt(delta), fSmallValue); // Arbitrary cutoff, but expectation cannot be negative or zero!
		double omegaDenominator = ((w->getExp[rowKey]-Lij)*(w->getExp[colKey]-Lij));
		if (omegaDenominator > fSmallValue ) { // Kind of arbitrary cutoff, but needs to be set at some low value to avoid Omega = Inf
			double Omega = (Lij*(w->L2x-w->getExp[rowKey]-w->getExp[colKey]+Lij)) / omegaDenominator;
			// TODO: check that the precision is set correctly?
			pvalueRow[i] = getPvalNCHGNew(w->dump, countRow[i], w->getMarginal[rowKey], w->getMarginal[colKey], w->N2x, Omega);
		}
		else {
			// NOTE: We need to keep this as the refined might have change the pvalue either way
			pvalueRow[i] = 1.0;
		}
	}
	
	/*
	 // WCH: but the default is already 1.0
	// the rest will have pvalue = 1.0
	for (map<SegmentKey, Interaction >::iterator it=w->matrix[rowKey].begin(); it!=w->segmentIterStarts[rowKey]; ++it) {
		it->second.pvalue = 1.0;
	}
	for (map<SegmentKey, Interaction >::iterator it=w->segmentIterEnds[rowKey]; it!=w->matrix[rowKey].end(); ++it) {
		it->second.pvalue = 1.0;
	}
	 */
	
}

// multi-threading version
void CCCMatrixInteraction::calculatePvalues(uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/)
{
	assert(isIntra);
	calculate_Pvalues_worker_t w(orderedSegmentMidPoints, mat_sk, mat_count, mat_pvalue, segmentIterStarts, segmentIterEnds, m_getExp, m_getMarginal, pDelta2ExpSpline);
	w.N2x = 2 * cachedN; w.L2x = 2.0 * cachedExpectN;
#if 0
	std::cerr << "getN: " << cachedN << std::endl;
	std::cerr << "getExpectN: " << cachedExpectN << std::endl;
#endif
#ifdef DUMP_CHROMOSOME
	if (DUMP_CHROMOSOME==m_chr) w.dump = true;
#endif
	kt_for(nThreads, calculate_Pvalues_worker, &w, (int)orderedSegments.size());
}

// NO optimization needed
void CCCMatrixInteraction::printPositives(ostream& ostr, string& chr, double alpha, uint32_t cutoff, bool printNij, bool printAll/*=false*/)
{
	for (SegmentKey rowKey=0; rowKey<orderedSegments.size(); ++rowKey) {
		SegmentKeyRowType& rowSegmentKeys = mat_sk[rowKey];
		CountRowType& countRow = mat_count[rowKey];
		PvalueRowType& pvalueRow = mat_pvalue[rowKey];
		for(uint32_t i=0; i<countRow.size(); ++i) {
			bool significant = (pvalueRow[i] <= alpha and countRow[i] >= cutoff);
			if (significant || printAll) {
				SegmentKey colKey = rowSegmentKeys[i];
				ostr << chr << "\t" << orderedSegments[rowKey].start << "\t" << orderedSegments[rowKey].end << "\t" <<  chr << "\t" << orderedSegments[colKey].start << "\t" << orderedSegments[colKey].end << "\t" << pvalueRow[i];
				if(printNij) {
					double expect = getExpectElement(rowKey, colKey);
					ostr << "\t" << countRow[i] << "\t" << ((expect<0) ? 0. : expect);
				}
				
				if (printAll) {
					ostr << "\t" << ((significant) ? "*" : ".");
				}
				ostr << endl;
			}
		}
	}
}


void CCCMatrixInteraction::setDeltaToExpectationMapper(spline& smoothDelta2ExpMapper, uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/)
{
	assert(isIntra);
	pDelta2ExpSpline = &smoothDelta2ExpMapper;
	minExpDelta = minDelta;
	maxExpDelta = maxDelta;

	calculate_cache_N(minDelta, maxDelta);

	cacheFlag &= (~CACHE_FLAG_EXPECTN);
	cachedExpectN = 0.0;
	calculate_cache_ExpectN(minDelta, maxDelta);
	
	// cache commonly used values
	cacheRowSums();
}

struct getExpectN_worker_t {
	getExpectN_worker_t(uint32_t* mps, SegmentKeyRowType& ss, SegmentKeyRowType& se, spline* sp, vector<double>& tr)
	: midpoints(mps), segmentKeyStarts(ss), segmentKeyEnds(se), delta2ExpSpline(sp), tmpRes(tr)
	{}
	
	uint32_t* midpoints;
	SegmentKeyRowType& segmentKeyStarts;
	SegmentKeyRowType& segmentKeyEnds;
	spline* delta2ExpSpline;
	vector<double>& tmpRes;
};

// [opt19.2: getExpectN_worker {perm}cell-dist-UT]
static void getExpectN_worker(void *data, int rowKey, int tid)
{
	getExpectN_worker_t *w = (getExpectN_worker_t*)data;
	double res = 0.0;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for (SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		res += w->delta2ExpSpline->cachedAt(w->midpoints[colKey]-rowMidPoint);
	}
	w->tmpRes[tid] += res;
	//w->tmpCount[tid] += w->segmentKeyEnds[rowKey] - w->segmentKeyStarts[rowKey];
}

// OPTIMIZE multi-threading
void CCCMatrixInteraction::calculate_cache_ExpectN(uint32_t minDelta/*=0*/, uint32_t maxDelta/*=MAXDELTA*/)
{
	if (CACHE_FLAG_EXPECTN!=(CACHE_FLAG_EXPECTN&cacheFlag)) {
		vector<double> tmpRes; tmpRes.resize(nThreads, 0.0);
		getExpectN_worker_t w(orderedSegmentMidPoints, segmentKeyStarts, segmentKeyEnds, pDelta2ExpSpline, tmpRes);
		kt_for(nThreads, getExpectN_worker, &w, (int)orderedSegments.size());
		cachedExpectN = accumulate(tmpRes.begin(), tmpRes.end(), (double)0.0);
		
		cacheFlag|=CACHE_FLAG_EXPECTN;
	}
	
#if 0
	cerr << "getExpectN: " << cachedExpectN << " count:" << count << endl;
#endif
}

void CCCMatrixInteraction::printExpectationMatrix(ostream& os,bool hideMask /*=false*/) {
	// print column header
	os << "Expect\t.\tColSum";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << m_getExp[segmentCounter];
	os << endl;
	
	os << ".\t.\tIndex";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << segmentCounter;
	os << endl;
	
	os << "RowSum\tIndex\tmidpt";
	for(SegmentKey segmentCounter=0; segmentCounter<orderedSegments.size(); ++segmentCounter) os << "\t" << orderedSegmentMidPoints[segmentCounter];
	os << endl;
	
	for (SegmentKey rowKey=0; rowKey<orderedSegments.size(); ++rowKey) {
		// WCH: decided that we should only print UT so that we don't need to double the index structure size!!
		os << m_getExp[rowKey] << "\t" << rowKey << "\t" << orderedSegmentMidPoints[rowKey];
		for (SegmentKey colKey=0; colKey<=rowKey; ++colKey) os << "\t.";
		for (SegmentKey colKey=rowKey+1; colKey<orderedSegments.size(); ++colKey) {
			uint32_t delta=orderedSegmentMidPoints[colKey]-orderedSegmentMidPoints[rowKey];
			if(delta>=minExpDelta and delta<=maxExpDelta) {
				double res = pDelta2ExpSpline->cachedAt(delta);
				os << "\t" << res;
			} else {
				os << "\t0";
			}
		}
		os << endl;
	}
}

struct expectation_RowSum_worker_t {
	expectation_RowSum_worker_t(uint32_t* mps, SegmentKeyRowType& ss, SegmentKeyRowType& se, SegmentKeyRowType& sps, SegmentKeyRowType& spe, spline* s, vector<double >& rs)
	: midpoints(mps), segmentKeyStarts(ss), segmentKeyEnds(se), segmentPrefixKeyStarts(sps), segmentPrefixKeyEnds(spe), pDelta2ExpSpline(s), rowSums(rs)
	{}
	
	uint32_t* midpoints;
	SegmentKeyRowType& segmentKeyStarts;
	SegmentKeyRowType& segmentKeyEnds;
	SegmentKeyRowType& segmentPrefixKeyStarts;
	SegmentKeyRowType& segmentPrefixKeyEnds;
	spline* pDelta2ExpSpline;
	vector<double >& rowSums;
};

// [opt19.2: expectation_RowSum_worker {perm}cell-dist-LT,UT]
static void expectation_RowSum_worker(void *data, int rowKey, int tid)
{
	expectation_RowSum_worker_t *w = (expectation_RowSum_worker_t*)data;
	
	double sum = 0.0f;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for (SegmentKey colKey=w->segmentPrefixKeyStarts[rowKey]; colKey<w->segmentPrefixKeyEnds[rowKey]; ++colKey) {
		sum += w->pDelta2ExpSpline->cachedAt(rowMidPoint-w->midpoints[colKey]);
	}
	for (SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		sum += w->pDelta2ExpSpline->cachedAt(w->midpoints[colKey]-rowMidPoint);
	}
	w->rowSums[rowKey] = sum;
}

void CCCMatrixInteraction::getExpectRowSums(vector<double > &res)
{
	res.clear();
	res.resize(orderedSegments.size());
	expectation_RowSum_worker_t w(orderedSegmentMidPoints, segmentKeyStarts, segmentKeyEnds, segmentPrefixKeyStarts, segmentPrefixKeyEnds, pDelta2ExpSpline, res);
	kt_for(nThreads, expectation_RowSum_worker, &w, (int)orderedSegments.size());
}


// multithreading version
struct interaction_FalsePositives_worker_t {
	interaction_FalsePositives_worker_t(uint32_t* mps, SegmentKeyRowType& ss, SegmentKeyRowType& se, vector<double >& ge, vector<uint32_t >& gm, spline* sp, vector<double >& rs)
	: midpoints(mps), segmentKeyStarts(ss), segmentKeyEnds(se), getExp(ge), getMarginal(gm), delta2ExpSpline(sp), rowSums(rs), dump(false)
	{}
	
	double alpha;
	uint32_t cutoff;
	
	uint32_t N2x;
	double L2x;

	uint32_t* midpoints;
	SegmentKeyRowType& segmentKeyStarts;
	SegmentKeyRowType& segmentKeyEnds;
	vector<double >& getExp;
	vector<uint32_t >& getMarginal;
	spline* delta2ExpSpline;
	vector<double >& rowSums;
	
	bool dump;
};

// [opt19.2: interaction_FalsePositives_worker {perm}cell-dist-UT]
static void interaction_FalsePositives_worker(void *data, int rowKey, int tid)
{
	interaction_FalsePositives_worker_t *w = (interaction_FalsePositives_worker_t*)data;
	double res = 0.0;
	uint32_t rowMidPoint = w->midpoints[rowKey];
	for(SegmentKey colKey=w->segmentKeyStarts[rowKey]; colKey<w->segmentKeyEnds[rowKey]; ++colKey) {
		uint32_t delta=w->midpoints[colKey]-rowMidPoint;
		//double Li = w->getExp[rowKey];
		//double Lj = w->getExp[colKey];
		//uint32_t Ni = w->getMarginal[rowKey];
		//uint32_t Nj = w->getMarginal[colKey];
		
		double Lij = max(w->delta2ExpSpline->cachedAt(delta), fSmallValue);
		
		double omegaDenominator = ((w->getExp[rowKey]-Lij)*(w->getExp[colKey]-Lij));
		if (omegaDenominator > fSmallValue ) { // Arbitrary cutoff to avoid Omega = Inf
			double Omega = (Lij*(w->L2x-w->getExp[rowKey]-w->getExp[colKey]+Lij)) / omegaDenominator;
			if (Omega<dMinOmega) {
				// need to perform separately
				int32_t quant = getQuantileNCHG(w->dump, 1-w->alpha, w->getMarginal[rowKey], w->getMarginal[colKey], w->N2x, Omega);
				if (quant >= 0) { // quant == -1 indicates memory issues
					uint32_t cut = (w->cutoff>quant) ? w->cutoff : 1+quant;  //max(w->cutoff,1+quant)
					res += getPvalNCHGNew(w->dump, cut, w->getMarginal[rowKey], w->getMarginal[colKey], w->N2x, Omega);
				}
			} else {
				// can re-use the FNCHG table
				double pval = getQuantilePvalNCHG(false/*w->dump*/, w->cutoff, 1-w->alpha, w->getMarginal[rowKey], w->getMarginal[colKey], w->N2x, Omega);
				if (pval>0.) res += pval;
			}
		}
	}
	
#ifdef FALSE_POSITIVE_SERIAL
	w->rowSums[rowKey] = res;
#else
	w->rowSums[tid] += res;
#endif
}

double CCCMatrixInteraction::estimateFalsePositives(double alpha, uint32_t cutoff, uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/)
{
	assert(isIntra);
	if(alpha == 0) {
		return 0;
	}
	
	vector<double> tmpSums;
#ifdef FALSE_POSITIVE_SERIAL
	tmpSums.resize(orderedSegments.size());
#else
	tmpSums.resize(nThreads);
#endif
	interaction_FalsePositives_worker_t w(orderedSegmentMidPoints, segmentKeyStarts, segmentKeyEnds, m_getExp, m_getMarginal, pDelta2ExpSpline, tmpSums);
	w.alpha = alpha; w.cutoff = cutoff;
	w.N2x = 2*cachedN; w.L2x = 2.0*cachedExpectN;
#if 0
	std::cerr << "getN: " << cachedN << std::endl;
	std::cerr << "getExpectN: " << cachedExpectN << std::endl;
#endif
#ifdef DUMP_CHROMOSOME
	if (DUMP_CHROMOSOME==m_chr) w.dump = true;
#endif
	kt_for(nThreads, interaction_FalsePositives_worker, &w, (int)orderedSegments.size());
	double res = accumulate(tmpSums.begin(), tmpSums.end(), (double)0.0);
	
	//cout << "FP:" <<  res << endl;
	return res;
}

uint32_t CCCMatrixInteraction::getSegmentCount()
{
	// FIXME: data type size
	return (uint32_t) orderedSegments.size();
}


void CCCMatrixInteraction::cacheRowSums()
{
	getExpectRowSums(m_getExp);
	getRowSums(m_getMarginal);
}

unsigned long CCCMatrixInteraction::getInteractionCount()
{
	return m_TotalInteractions;
}

void CCCMatrixInteraction::setThreads(int n)
{
	nThreads = n;
}

uint32_t CCCMatrixInteraction::getWidestInteractionSpan()
{
	if (orderedSegments.size()>0) {
		return orderedSegmentMidPoints[orderedSegments.size()-1] - orderedSegmentMidPoints[0] + 1;
	} else {
		return 0;
	}
}

