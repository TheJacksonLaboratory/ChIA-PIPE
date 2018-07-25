// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen


#include "CCCStatistics.h"
#include <iomanip>
#include <iostream>
#include <fstream>

#include <assert.h>
#include <errno.h>

#include "bar.h"

using namespace std;

extern "C" {
	double cputime();
	double realtime();
}

// TODO: WCH [pending]
void printPositives(ostream& ostr, string& chr, CCCMatrixInteraction& cpemat, double alpha, uint32_t cutoff, bool printNij /*=false*/, bool printAll /*=false*/) {
	assert(cpemat.isIntra);
	cpemat.printPositives(ostr, chr, alpha, cutoff, printNij, printAll);
}


// Utils:
// WCH: loop does not terminate
template<typename T> vector<T> getSequence(T from, T to, uint32_t length) {
#if 0
	assert(to > from);
#else
	assert(to >= from);
#endif
	vector<T> res;
	T stepSize=(to-from)/(length-1);
	T i = from;
	for(uint32_t c=0; c<length; c++) {
		res.push_back(i);
		i+=stepSize;
	}
	return res;
}

pair<double, double> getAlphaRange(map<double, double> alpha2FDR, double minFDR) {
	typedef map<double, double>::iterator it2;
	//WCH
	double alpha=0, FDR=0;
	pair<double, double> res;
	res.first = 0;
	res.second = 0;
	for(it2 iter = alpha2FDR.begin(); iter != alpha2FDR.end(); iter++) {
		alpha = iter->first;
		FDR = iter->second;
		//cout << "alpha/FDR: " << alpha << " " << FDR << endl;
		fprintf(stderr, "[M::%s:%d] Re-estimating alpha/FDR=%.6e, FDR=%.6e%s\n", __func__, __LINE__, alpha, FDR, (FDR<minFDR)?" (*)":"");
		if(FDR > minFDR) {
			if(res.first != 0 ) {
				res.second = alpha;
				fprintf(stderr, "[M::%s:%d] success! FDR(%.6e)>minFDR(%.6e)\n", __func__, __LINE__, FDR, minFDR);
				return res; // success! minFDR is within the range
			}
			else { // minFDR is not within the range. C
				res.first = 0;
				res.second = alpha;
				fprintf(stderr, "[M::%s:%d] minFDR(%.6e) is not within range\n", __func__, __LINE__, minFDR);
				return res; // the right alpha must be somewhere between 0 and minFDR
			}
		}
		res.first = alpha;
	}
	// At this point, the right alpha must be somewhere between alpha and minFDR
	
	fprintf(stderr, "[M::%s:%d] right alpha is between alpha(%.6e) and minFDR(%.6e)\n", __func__, __LINE__, alpha, minFDR);
	res.second = minFDR;
	return res;
}

// TODO: WCH [coded]
static int G_nThreads = 1;
void setThreads(int n) {
	G_nThreads = n;
}

void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);

bool chromSegmentCountDesc(const pair<ChromosomeIndexType, uint32_t >& a, const pair<ChromosomeIndexType, uint32_t >& b) {
	return b.second < a.second;
}

// this routine will estimate the largest alpha such that the FDR<fdr
map<double, double>  estimateFDRAlpha(const char *szDataFile, vector<string>& chromosomes, vector<pair<ChromosomeIndexType, uint64_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, double fdr, uint32_t refinesteps, bool noEarlyTermination/*=false*/, uint32_t cutoff /*=3*/, uint32_t deltaCutoff /*=8000*/, uint32_t maxDelta /*=MAXDELTA*/) {
	
	double t_diff;
	double t_start = realtime();
	double t_start_iter;

	map<double, pair<uint64_t, double> > tmp;
	map<double, double> res;
	
	double alphaLower = 0.0;
	double alphaUpper = fdr / 100.0;
	double alpha = alphaUpper;
	uint32_t numStep = 1;
	uint64_t maxAlphaPositives = 0;
	
	uint64_t numPositives = 0;
	double numFalsePositives = 0.0;
	bool maxAlphaPositivesIncreasing = true;
	uint64_t cutoffMaxPositive = 0;
	
	string szTmpDataFile(szDataFile); szTmpDataFile.append(".tmp");
	
	// let's attempt to restore if we may
	bool bResumedOperation = false;
	{
		ifstream fdrPersistence(szDataFile, ios::in);
		if (fdrPersistence.is_open()) {
			double savedFDR;
			fdrPersistence >> savedFDR;
			uint32_t savedNumStep;
			fdrPersistence >> savedNumStep;
			
			// TODO: there is a need to set alpha, alphaLower, alphaUpper, maxAlphaPositives
			uint32_t numFDRs = 0;
			while (true) {
				// read in the line of data
				double savedAlpha;
				uint64_t savedNumPositives;
				double savedNumFalsePositives;
				double savedFDR;
				fdrPersistence >> savedAlpha >> savedNumPositives >> savedNumFalsePositives >> savedFDR;
				if (fdrPersistence.eof()) break;
				// process the line of data
				tmp[savedAlpha] = make_pair(savedNumPositives, savedNumFalsePositives);
				// TODO: if there is NO savedFDR
				double alphaFDR = (0==savedNumPositives) ? 0.0 : savedNumFalsePositives / savedNumPositives;
				res[savedAlpha] = alphaFDR;
				numFDRs++;
			}
			
			// TODO: let's set up the progressive logic to make sure that we do NOT perform an additional run!
			if (numFDRs>0) {
				if (numFDRs != savedNumStep) {
					savedNumStep = numFDRs;
				}
				
				double savedAlphaLower = 0.0;
				double savedAlphaUpper = fdr; //WCH: original is fdr / 100.0;
				uint64_t savedMaxAlphaPositives = 0;
				for(map<double, pair<uint64_t, double> >::iterator iter = tmp.begin(); iter != tmp.end(); iter++) {
					double savedAlpha = iter->first;
					uint64_t savedNumPositives = iter->second.first;
					double savedNumFalsePositives = iter->second.second;
					double savedFDR = (0==savedNumPositives) ? 0.0 : savedNumFalsePositives / savedNumPositives;
					if (savedFDR<fdr) {
						if (savedAlpha>savedAlphaLower) {
							savedAlphaLower = savedAlpha;
							if (savedNumPositives>savedMaxAlphaPositives)
								savedMaxAlphaPositives = savedNumPositives;
							else
								maxAlphaPositivesIncreasing = false;
						}
					} else {
						if (savedAlpha<savedAlphaUpper) {
							savedAlphaUpper = savedAlpha;
							if (savedNumPositives<=savedMaxAlphaPositives)
								maxAlphaPositivesIncreasing = false;
							break;
						}
					}
					
				}
				
				// set up the bounds and conditions
				alphaLower = savedAlphaLower;
				alphaUpper = savedAlphaUpper;
				numStep = savedNumStep;
				maxAlphaPositives = savedMaxAlphaPositives;
				alpha = 0.5 * (alphaLower + alphaUpper);
				
				bResumedOperation = true;
				if (maxAlphaPositivesIncreasing) {
					// there are more work to do
					fprintf(stderr, "[M::%s:%d] estimateFDR %u iteraction results RESTORED.\n", __func__, __LINE__, numStep);
				} else {
					// result improvement impossible
					fprintf(stderr, "[M::%s:%d] estimateFDR %u iteraction results RESTORED. NO further iteraction needed.\n", __func__, __LINE__, numStep);
				}
			}
			
		}
		fdrPersistence.close();
	}
	
	// get the first FDR estimation
	if (!bResumedOperation) {
		t_start_iter = realtime();
		numPositives = 0;
		numFalsePositives = 0.0;
		for(vector<pair<ChromosomeIndexType, uint64_t > >::iterator itChrom = chromosomesRows.begin(); itChrom!=chromosomesRows.end(); ++itChrom) {
			double t_chrStart = realtime();
			numPositives += cmpmat[itChrom->first].getNumberOfPositives(alpha, cutoff);
			numFalsePositives += cmpmat[itChrom->first].estimateFalsePositives(alpha, cutoff, deltaCutoff, maxDelta);
			t_diff = realtime() - t_chrStart;
			fprintf(stderr, "[M::%s:%d] estimateFDR %s FP:%.2f P:%llu rtmin:%.2f..\n", __func__, __LINE__, chromosomes[itChrom->first].c_str(), numFalsePositives, numPositives, t_diff/60.0);
		}
		
		tmp[alpha] = make_pair(numPositives, numFalsePositives);
		double alphaFDR = (0==numPositives) ? 0.0 : numFalsePositives / numPositives;
		res[alpha] = alphaFDR;
		
		t_diff = realtime() - t_start_iter;
		fprintf(stderr, "[M::%s:%d] step:%d alpha:%.6e FDR:%.6e FP:%.2f P:%llu rtmin:%.2f%s\n", __func__, __LINE__, numStep, alpha, alphaFDR, numFalsePositives, numPositives, t_diff/60.0, ((alphaFDR<fdr)?" (*)":""));
		
		// let's attempt to persist for resumability
		{
			ofstream fdrPersistence(szTmpDataFile.c_str(), ios::out);
			fdrPersistence << scientific << setprecision(6) << fdr << endl;
			fdrPersistence << numStep << endl;
			for(map<double, pair<uint64_t, double> >::iterator iter = tmp.begin(); iter != tmp.end(); iter++) {
				double savedAlpha = iter->first;
				uint64_t savedNumPositives = iter->second.first;
				double savedNumFalsePositives = iter->second.second;
				double savedFDR = (0==savedNumPositives) ? 0.0 : savedNumFalsePositives / savedNumPositives;
				fdrPersistence << scientific << setprecision(6) << savedAlpha << "\t" << savedNumPositives << "\t" << scientific << setprecision(6) << savedNumFalsePositives << "\t" << scientific << setprecision(6) << savedFDR << endl;
			}
			fdrPersistence.close();
			if (0!=remove(szDataFile)) {
				if (2!=errno)
					fprintf(stderr, "[M::%s:%d] step:%d FAILED to remove earlier progress file, errno=%d.\n", __func__, __LINE__, numStep, errno);
			}
			if (0!=rename(szTmpDataFile.c_str(), szDataFile)) {
				fprintf(stderr, "[M::%s:%d] step:%d FAILED to rename latest progress file, errno=%d.\n", __func__, __LINE__, numStep, errno);
			} else {
				fprintf(stderr, "[M::%s:%d] step:%d progress saved.\n", __func__, __LINE__, numStep);
			}
		}
		
		if (alphaFDR >= fdr) {
			// we still have much larger FDR, so we need a more stringent alpha
			alphaUpper = alpha;
			alpha = 0.5 * (alphaLower + alphaUpper);
			
			// new smallest positives which FDR is exceeded
		} else if (alphaFDR < fdr) {
			// there are more positives that we can report, relax the alpha
			alphaLower = alpha;
			alphaUpper = fdr; // we set the Upper bound to a higher alpha
			alpha = 0.5 * (alphaLower + alphaUpper);
			
			// new lower bound of the attainable positives
			maxAlphaPositives = numPositives;
		}
	}
	
	if (numStep<=refinesteps && maxAlphaPositivesIncreasing) {
		// attempt to establish possible range to start iteration
		t_start_iter = realtime();
		// let's attempt to estimate a better alpha
		cutoffMaxPositive = 0;
		for(vector<pair<ChromosomeIndexType, uint64_t > >::iterator itChrom = chromosomesRows.begin(); itChrom!=chromosomesRows.end(); ++itChrom) {
			double t_chrStart = realtime();
			cmpmat[itChrom->first].getPositivesBound(cutoff, cutoffMaxPositive);
			t_diff = realtime() - t_chrStart;
			fprintf(stderr, "[M::%s:%d] estimateAlpha %s cutoff=%u positives=%llu rtmin:%.2f\n", __func__, __LINE__, chromosomes[itChrom->first].c_str(), cutoff, cutoffMaxPositive, t_diff/60.0);
		}
		t_diff = realtime() - t_start_iter;
		fprintf(stderr, "[M::%s:%d] estimateAlpha genome took %.2f\n", __func__, __LINE__, t_diff/60.0);
	}
	
	// iterate on FDR estimation until the specified number of steps or no better result possible
	numStep++;
	while(numStep<=refinesteps && maxAlphaPositivesIncreasing) {
		
		double t_start_iter = realtime();
		
		uint64_t numPositives = 0;
		double numFalsePositives = 0.0;
		
		for(vector<pair<ChromosomeIndexType, uint64_t > >::iterator itChrom = chromosomesRows.begin(); itChrom!=chromosomesRows.end(); ++itChrom) {
			double t_chrStart = realtime();
			numPositives += cmpmat[itChrom->first].getNumberOfPositives(alpha, cutoff);
			numFalsePositives += cmpmat[itChrom->first].estimateFalsePositives(alpha, cutoff, deltaCutoff, maxDelta);
			t_diff = realtime() - t_chrStart;
			fprintf(stderr, "[M::%s:%d] estimateFDR %s FP:%.2f P:%llu rtmin:%.2f..\n", __func__, __LINE__, chromosomes[itChrom->first].c_str(), numFalsePositives, numPositives, t_diff/60.0);
		}
		
		tmp[alpha] = make_pair(numPositives, numFalsePositives);
		double alphaFDR = (0==numPositives) ? 0.0 : numFalsePositives / numPositives;
		res[alpha] = alphaFDR;
		
		t_diff = realtime() - t_start_iter;
		fprintf(stderr, "[M::%s:%d] step:%u alpha:%.6e FDR:%.6e FP:%.2f P:%llu rtmin:%.2f%s\n", __func__, __LINE__, numStep, alpha, alphaFDR, numFalsePositives, numPositives, t_diff/60.0, ((alphaFDR<fdr)?" (*)":""));
		
		// let's attempt to persist for resumability
		{
			ofstream fdrPersistence(szTmpDataFile.c_str(), ios::out);
			fdrPersistence << scientific << setprecision(6) << fdr << endl;
			fdrPersistence << numStep << endl;
			for(map<double, pair<uint64_t, double> >::iterator iter = tmp.begin(); iter != tmp.end(); iter++) {
				double savedAlpha = iter->first;
				uint64_t savedNumPositives = iter->second.first;
				double savedNumFalsePositives = iter->second.second;
				double savedFDR = (0==savedNumPositives) ? 0.0 : savedNumFalsePositives / savedNumPositives;
				fdrPersistence << scientific << setprecision(6) << savedAlpha << "\t" << savedNumPositives << "\t" << scientific << setprecision(6) << savedNumFalsePositives << "\t" << scientific << setprecision(6) << savedFDR << endl;
			}
			fdrPersistence.close();
			if (0!=remove(szDataFile)) {
				if (2!=errno)
					fprintf(stderr, "[M::%s:%d] step:%d FAILED to remove earlier progress file, errno=%d.\n", __func__, __LINE__, numStep, errno);
			}
			if (0!=rename(szTmpDataFile.c_str(), szDataFile)) {
				fprintf(stderr, "[M::%s:%d] step:%d FAILED to rename latest progress file, errno=%d.\n", __func__, __LINE__, numStep, errno);
			} else {
				fprintf(stderr, "[M::%s:%d] step:%d progress saved.\n", __func__, __LINE__, numStep);
			}
		}
		
		// let's decide on the next alpha to try
		if (alphaFDR >= fdr) {
			// we still have much larger FDR, so we need a more stringent alpha
			alphaUpper = alpha;
			alpha = 0.5 * (alphaLower + alphaUpper);
		} else if (alphaFDR < fdr) {
			// there are more positives that we can report, relax the alpha
			alphaLower = alpha;
			alpha = 0.5 * (alphaLower + alphaUpper);
			// early termination
			if (numPositives>maxAlphaPositives) {
				maxAlphaPositives = numPositives;
				if (maxAlphaPositives==cutoffMaxPositive) {
					fprintf(stderr, "[M::%s:%d] #Positives reached maximum (%llu) attainable for cutoff>=%u at step#%u\n", __func__, __LINE__, cutoffMaxPositive, cutoff, numStep);
					
					if (!noEarlyTermination) maxAlphaPositivesIncreasing = false;
					if (numStep<refinesteps)
						fprintf(stderr, "[M::%s:%d] Results will not improve further. Terminating iterative process..\n", __func__, __LINE__);
				}
			} else {
				// TODO: for testing purpsoes; TO turn on afterward
				if (!noEarlyTermination) maxAlphaPositivesIncreasing = false;
				if (numStep<refinesteps)
					fprintf(stderr, "[M::%s:%d] Results will not improve further. Terminating iterative process..\n", __func__, __LINE__);
			}
		}
		
		numStep++;
	}
	
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] estimateFDR took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	
	return res;
}


struct delta_quantile_precompute_worker_t {
	delta_quantile_precompute_worker_t(uint32_t *pds, uint32_t cs, uint32_t qi)
	: pDeltas(pds), chunkStart(cs), quantileIndex(qi)
	{}
	
	uint32_t *pDeltas;
	uint32_t chunkStart;
	uint32_t quantileIndex;
};

struct delta_quantile_chunk_precompute_worker_t {
	delta_quantile_chunk_precompute_worker_t(uint32_t *pds, vector<uint32_t>& cs, vector<uint32_t>& ce)
	: pDeltas(pds), chunkStarts(cs), chunkEnds(ce)
	{}
	
	uint32_t *pDeltas;
	vector<uint32_t>& chunkStarts;
	vector<uint32_t>& chunkEnds;
};

/*
 static void delta_quantile_chunk_precompute_worker(void *data, int idx, int tid)
{
	delta_quantile_chunk_precompute_worker_t *w = (delta_quantile_chunk_precompute_worker_t*)data;
	for(uint32_t i=w->chunkStarts[idx]; i<=w->chunkEnds[idx]; ++i)
		w->pDeltas[i] = idx;
}
*/

struct delta_quantile_statistics_worker_t {
	delta_quantile_statistics_worker_t(uint32_t *pds, uint32_t *pss, vector<uint32_t>& cs, vector<uint32_t>& ce, vector<DeltaStatistics >& dss)
	: pDeltas(pds), pSums(pss), chunkStarts(cs), chunkEnds(ce), deltasStats(dss)
	{}
	
	uint32_t *pDeltas;
	uint32_t *pSums;
	vector<uint32_t>& chunkStarts;
	vector<uint32_t>& chunkEnds;
	
	vector<DeltaStatistics >& deltasStats;
};

static void delta_quantile_statistics_worker(void *data, int idx, int tid)
{
	delta_quantile_statistics_worker_t *w = (delta_quantile_statistics_worker_t*)data;
	for(uint32_t i=w->chunkStarts[idx]; i<=w->chunkEnds[idx]; ++i) {
		w->deltasStats[idx].count += w->pDeltas[i];
		w->deltasStats[idx].sum += w->pSums[i];
	}
}

// experimental function: under optimization testing
void estimateDeltaSpline(const char *szDataFile, spline& s, vector<double>& sx, DeltaCountSums& deltaCountSums, vector<string>& chromosomes, vector<pair<ChromosomeIndexType, uint64_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, uint32_t minDelta /*=8000*/, uint32_t maxDelta /*=MAXDELTA*/, bool masking /*=false*/, uint32_t percentiles /*=1000*/, bool dumpSpline/*=false*/) {

	ofstream fileDump;
	if (dumpSpline) fileDump.open(masking ? "deltas.refined.xls" : "deltas.exact.xls");
	
	double t_diff;
	// Get all deltas for all chromosomes:
	QuantileMapper quantileMapper;
	// deltas allocations
	if (0==deltaCountSums.nMaxSpan)
	{
		uint32_t nMaxSpan = 0;
		for(ChromosomeIndexType i=0; i<chromosomes.size(); ++i) {
			uint32_t maxChromSpan = cmpmat[i].getWidestInteractionSpan();
			if (maxChromSpan>nMaxSpan) nMaxSpan = maxChromSpan;
		}
		nMaxSpan++; // 0-based index, so need 1 more allocation
		fprintf(stderr, "[M::%s:%d] max span = %u\n", __func__, __LINE__, nMaxSpan);
		
		// NOTE: deltas stores count, the index is actually the span
		deltaCountSums.init(nMaxSpan);
	}
	
	double t_start = realtime();
	// attempt to restore from persistence
	// in event that persistence is incomplete, repeat the process
	bool bResumedOperation = false;
	if (!masking) {
		ifstream deltaPersistence(szDataFile, ios::in | ios::binary);
		if (deltaPersistence.is_open()) {
			// for sanity check
			uint32_t nMaxSpan; deltaPersistence.read((char *)(&nMaxSpan), sizeof(nMaxSpan));
			if (deltaPersistence) {
				if (deltaCountSums.nMaxSpan==nMaxSpan) {
					// read the rest of the summary
					deltaPersistence.read((char *)(&(deltaCountSums.ulTotalDeltas)), sizeof(deltaCountSums.ulTotalDeltas));
					if (deltaPersistence) {
						// read the rest of the elements value
						deltaPersistence.read((char *)(deltaCountSums.pDeltas), nMaxSpan*sizeof(*deltaCountSums.pDeltas));
						if (deltaPersistence) {
							deltaPersistence.read((char *)(deltaCountSums.pSums), nMaxSpan*sizeof(*deltaCountSums.pSums));
							if (deltaPersistence) {
								// sanity check on BIT_ARRAY
								bit_index_t num_of_bits;
								deltaPersistence.read((char *)(&num_of_bits), sizeof(num_of_bits));
								if (deltaPersistence) {
									if (deltaCountSums.deltaPresences->num_of_bits==num_of_bits) {
										word_addr_t num_of_words;
										deltaPersistence.read((char *)(&num_of_words), sizeof(num_of_words));
										if (deltaPersistence) {
											if (deltaCountSums.deltaPresences->num_of_words==num_of_words) {
												// okie to read in the data now!
												deltaPersistence.read((char *)(deltaCountSums.deltaPresences->words), deltaCountSums.deltaPresences->num_of_words*sizeof(*deltaCountSums.deltaPresences->words));
												if (deltaPersistence) {
													bResumedOperation = true;
												} else {
													// fail to read all the deltas presence bits
													fprintf(stderr, "[M::%s:%d] FAILED to restore Deltas' presence bits. Will recompute..\n", __func__, __LINE__);
												}
											} else {
												// bits configration does not agree!!
												fprintf(stderr, "[M::%s:%d] FAILED to restore. Inconsistent word configuration (recorded: %llu vs computed: %llu). Will recompute..\n", __func__, __LINE__, num_of_words, deltaCountSums.deltaPresences->num_of_words);
											}
										} else {
											// fail to read all the deltas presence word configuration
											fprintf(stderr, "[M::%s:%d] FAILED to restore Deltas' word configuration. Will recompute..\n", __func__, __LINE__);
										}
									} else {
										// bits configration does not agree!!
										fprintf(stderr, "[M::%s:%d] FAILED to restore. Inconsistent bit configuration (recorded: %llu vs computed: %llu). Will recompute..\n", __func__, __LINE__, num_of_bits, deltaCountSums.deltaPresences->num_of_bits);
									}
								} else {
									// fail to read all the deltas presence bit configuration
									fprintf(stderr, "[M::%s:%d] FAILED to restore Deltas' bit configuration. Will recompute..\n", __func__, __LINE__);
								}
							} else {
								// fail to read all the deltas sum
								fprintf(stderr, "[M::%s:%d] FAILED to restore Deltas' sum. Will recompute..\n", __func__, __LINE__);
							}
						} else {
							// fail to read all the deltas count
							fprintf(stderr, "[M::%s:%d] FAILED to restore Deltas' count. Will recompute..\n", __func__, __LINE__);
						}
					} else {
						// fail to read the total number of deltas
						fprintf(stderr, "[M::%s:%d] FAILED to restore TotalDeltas. Will recompute..\n", __func__, __LINE__);
					}
				} else {
					// MaxSpan does not agree!!
					fprintf(stderr, "[M::%s:%d] FAILED to restore. Inconsistent MaxSpan (recorded: %u vs computed: %u). Will recompute..\n", __func__, __LINE__, nMaxSpan, deltaCountSums.nMaxSpan);
				}
			} else {
				// fail to read the number of unique deltas
				fprintf(stderr, "[M::%s:%d] FAILED to restore MaxSpan. Will recompute..\n", __func__, __LINE__);
			}
			deltaPersistence.close();
			
			// if we cannot resume, it is safer to re-initialize
			if (!bResumedOperation) {
				deltaCountSums.terminate();
				deltaCountSums.init(nMaxSpan);
			} else {
				t_diff = realtime() - t_start;
				fprintf(stderr, "[M::%s:%d] #deltas(genome)=%llu RESTORED, %.2f min..\n", __func__, __LINE__, deltaCountSums.ulTotalDeltas, t_diff/60.0);
			}
		}
	}
	
	// work horse: getting all the possible deltas from the contact matrix
	if (!bResumedOperation) {
		t_start = realtime();
		deltaCountSums.ulTotalDeltas = 0;
		if (masking) {
			for(ChromosomeIndexType i=0; i<chromosomes.size(); ++i) {
				double t_subStart = realtime();
				ChromosomeIndexType nChrIdx = chromosomesRows[i].first;
				uint64_t ulDeltas = cmpmat[nChrIdx].maskDeltasCountSum(deltaCountSums.pDeltas, deltaCountSums.pSums, deltaCountSums.deltaPresences, deltaCountSums.nMaxSpan+1, minDelta, maxDelta);
				deltaCountSums.ulTotalDeltas += ulDeltas;
				t_diff = realtime() - t_subStart;
				fprintf(stderr, "[M::%s:%d] #deltas(%s)=%llu collected, %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), ulDeltas, t_diff/60.0);
			}
		} else {
			for(ChromosomeIndexType i=0; i<chromosomes.size(); ++i) {
				double t_subStart = realtime();
				ChromosomeIndexType nChrIdx = chromosomesRows[i].first;
				uint64_t ulDeltas = cmpmat[nChrIdx].getDeltasCountSum(deltaCountSums.pDeltas, deltaCountSums.pSums, deltaCountSums.deltaPresences, deltaCountSums.nMaxSpan+1, minDelta, maxDelta, masking);
				deltaCountSums.ulTotalDeltas += ulDeltas;
				t_diff = realtime() - t_subStart;
				fprintf(stderr, "[M::%s:%d] #deltas(%s)=%llu collected, %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), ulDeltas, t_diff/60.0);
			}
		}
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] #deltas(genome)=%llu collected, %.2f min..\n", __func__, __LINE__, deltaCountSums.ulTotalDeltas, t_diff/60.0);
	}
	
	// Let's save these delta summary into persistence
	if (!(masking || bResumedOperation)) {
		t_start = realtime();
		
		bool persisted = false;
		ofstream deltaPersistence(szDataFile, ios::out | ios::binary);
		// 1st object
		if (deltaPersistence.write((char *)(&(deltaCountSums.nMaxSpan)), sizeof(deltaCountSums.nMaxSpan))) {
			// 2nd object
			if (deltaPersistence.write((char *)(&(deltaCountSums.ulTotalDeltas)), sizeof(deltaCountSums.ulTotalDeltas))) {
				// 3rd object - deltas
				if (deltaPersistence.write((char *)(deltaCountSums.pDeltas), deltaCountSums.nMaxSpan*sizeof(*deltaCountSums.pDeltas))) {
					// 3rd object - sums
					if (deltaPersistence.write((char *)(deltaCountSums.pSums), deltaCountSums.nMaxSpan*sizeof(*deltaCountSums.pSums))) {
						// 3rd object - msking
						if (deltaPersistence.write((char *)(&(deltaCountSums.deltaPresences->num_of_bits)), sizeof(deltaCountSums.deltaPresences->num_of_bits))) {
							if (deltaPersistence.write((char *)(&(deltaCountSums.deltaPresences->num_of_words)), sizeof(deltaCountSums.deltaPresences->num_of_words))) {
								if (deltaPersistence.write((char *)(deltaCountSums.deltaPresences->words), deltaCountSums.deltaPresences->num_of_words*sizeof(*deltaCountSums.deltaPresences->words))) {
									persisted = true;
								} else {
									fprintf(stderr, "[M::%s:%d] FAILED to save Deltas' presence bits. Non-resumable in future..\n", __func__, __LINE__);
								}
							} else {
								fprintf(stderr, "[M::%s:%d] FAILED to save Deltas' word configuration. Non-resumable in future..\n", __func__, __LINE__);
							}
						} else {
							fprintf(stderr, "[M::%s:%d] FAILED to save Deltas' bit configuration. Non-resumable in future..\n", __func__, __LINE__);
						}
					} else {
						fprintf(stderr, "[M::%s:%d] FAILED to save Deltas' sum. Non-resumable in future..\n", __func__, __LINE__);
					}
				} else {
					fprintf(stderr, "[M::%s:%d] FAILED to save Deltas' count. Non-resumable in future..\n", __func__, __LINE__);
				}
			} else {
				fprintf(stderr, "[M::%s:%d] FAILED to save TotalDeltas. Non-resumable in future..\n", __func__, __LINE__);
			}
		} else {
			fprintf(stderr, "[M::%s:%d] FAILED to save MaxSpan. Non-resumable in future..\n", __func__, __LINE__);
		}
		deltaPersistence.close();

		if (persisted) {
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] #deltas(genome)=%llu PERSISTED, %.2f min..\n", __func__, __LINE__, deltaCountSums.ulTotalDeltas, t_diff/60.0);
		}
	}
	
	/*
	// DEBUG: WCH - for deltas checking
	{
		ofstream fileDump2;
		fileDump2.open(masking ? "deltas.refined.xls.baseline.opt19.4.14n" : "deltas.exact.xls.baseline.opt19.4.14n");
		//unsigned int z=0;
		unsigned long nCurrCount = 0;
		for(int i=0; nCurrCount<ulTotalDeltas; ++i) {
			if (pDeltas[i]>0) {
				//for(int x=0; x<pDeltas[i]; x++) {
					fileDump2 << i << "\t" << pDeltas[i] << endl;
					//z++;
					//if (1000000==z) break;
				//}
				//if (1000000==z) break;
				nCurrCount += pDeltas[i];
			}
		}
		fileDump2.close();
	}
	// END - DEBUG: WCH - for deltas checking
	*/
	
	// NOTE: pDeltas needs to hold the count for quantile computation to work
	// NOTE: pDelta = indicator + count; but ONLY count is needed here!!!
	double t_subStart = realtime();
	quantileMapper.computeDeltaCountQuantile(deltaCountSums.ulTotalDeltas, deltaCountSums.pDeltas, percentiles, (dumpSpline ? &fileDump : NULL));
	t_diff = realtime() - t_subStart;
	fprintf(stderr, "[M::%s:%d] delta quantiles computation deltas took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	
	if (dumpSpline)
	{
		fileDump << endl;
		fileDump << "Deltas statistics" << endl;
		fileDump << "quantileUpper\tquantileRep" << endl;
		vector<uint32_t>& quantileDeltas = quantileMapper.getQuantileDeltas();
		for(vector<uint32_t>::iterator it=quantileDeltas.begin(); it!=quantileDeltas.end(); ++it) {
			fileDump << (*it) << "\t" << quantileMapper.getDeltaQuantile(*it) << endl;
		}
	}
	
	t_start = realtime();
	vector<DeltaStatistics > deltasStats;
	{
		vector<uint32_t>& quantileDeltas = quantileMapper.getQuantileDeltas();
		deltasStats.resize(quantileDeltas.size());
		// multithreading version
		vector<uint32_t> chunkStarts; chunkStarts.resize(quantileDeltas.size());
		chunkStarts[0] = minDelta;
		for(uint32_t idx=1; idx<quantileDeltas.size(); ++idx) {
			chunkStarts[idx] = quantileDeltas[idx-1]+1;
		}
		// NOTE: pDelta = indicator + count ==> indicator + qunatile
		delta_quantile_statistics_worker_t w(deltaCountSums.pDeltas, deltaCountSums.pSums, chunkStarts, quantileDeltas, deltasStats);
		kt_for(G_nThreads, delta_quantile_statistics_worker, &w, (int)quantileDeltas.size());
	}
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] delta quantiles statistics computation (genome) took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	
	if (dumpSpline)
	{
		fileDump << endl;
		fileDump << "Deltas statistics" << endl;
		fileDump << "deltaRep\tsum\tcount\tmean" << endl;
		vector<uint32_t>& quantiless = quantileMapper.getQuantiles();
		for(uint32_t idx=0; idx<quantiless.size(); ++idx) {
			uint32_t deltaQuantileRep = quantiless[idx];
			fileDump << deltaQuantileRep << "\t" << deltasStats[idx].sum << "\t" << deltasStats[idx].count << "\t" << deltasStats[idx].getMean() << endl;
		}
	}
	
	//vector<double> delta_x; // x-object for spline generation
	sx.clear();
	vector<double> exp_y; // y-object for spline generation
	{
		vector<uint32_t>& quantiless = quantileMapper.getQuantiles();
		for(uint32_t idx=0; idx<quantiless.size(); ++idx) {
			sx.push_back((double)(quantiless[idx]));
			exp_y.push_back((double)(deltasStats[idx].getMean()));
		}
	}

	//spline s;
	// If x is not sorted, an error will occur.
	//s.set_points(delta_x,exp_y);    // currently it is required that X is already sorted
	s.set_points(sx,exp_y);    // currently it is required that X is already sorted

	// NOTE: pDeltas needs to indicate if the precomputation for a particular delta value should be carried out
	// pre-compute the appropriate deltas
	t_subStart = realtime();
	s.precompute(deltaCountSums.nMaxSpan, deltaCountSums.deltaPresences);
	t_diff = realtime() - t_subStart;
	fprintf(stderr, "[M::%s:%d] delta expectation pre-computation took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	//free(pDeltas);
	//free(pSums);
	//bardestroy(deltaPresences);
	
	if (dumpSpline) {
		fileDump << endl;
		fileDump << "Spline data points" << endl;
		fileDump << "x\ty\ta\tb\tc" <<  endl;
		for(vector<double>::size_type e = 0 ; e<s.m_x.size(); ++e) {
			fileDump << s.m_x[e] << "\t" << s.m_y[e] << "\t" << s.m_a[e] << "\t" << s.m_b[e] << "\t" << s.m_c[e] << endl;
		}
		
		fileDump.close();
	}
}

void estimateResources(string& filename, vector<string>& chromosomes, vector<pair<ChromosomeIndexType, uint64_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/, bool masking/*=false*/)
{
	double t_diff;
	double t_start = realtime();
	
	// Get all deltas for all chromosomes:
	vector<uint64_t> chromDeltaCounts;
	uint64_t ulTotalDeltas = 0;
	for(ChromosomeIndexType i=0; i<chromosomes.size(); ++i) {
		double t_subStart = realtime();
		ChromosomeIndexType nChrIdx = chromosomesRows[i].first;
		uint64_t ulDeltas = cmpmat[nChrIdx].getDeltasCount(NULL, NULL, 0, minDelta, maxDelta, masking);
		ulTotalDeltas += ulDeltas;
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] #deltas(%s)=%llu collected, %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), ulDeltas, t_diff/60.0);
	}
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] #deltas(genome)=%llu collected, %.2f min..\n", __func__, __LINE__, ulTotalDeltas, t_diff/60.0);
	
	// loop as table grid!!!
	{
		string ofile(filename); ofile.append(".resource.xls");
		ofstream fileDump;
		fileDump.open(ofile.c_str());
		
		unsigned long totalSegments = 0;
		unsigned long totalInteractions = 0;
		unsigned long totalDeltas = 0;
		fileDump << "#filename\t" << filename.c_str() << endl;
		fileDump << "#chr\tanchors\tinteractions\tdeltas" << endl;
		for(ChromosomeIndexType i=0; i<chromosomes.size(); ++i) {
			totalSegments += cmpmat[i].getSegmentCount();
			totalInteractions += cmpmat[i].getInteractionCount();
			totalDeltas += chromDeltaCounts[i];
			fileDump << chromosomes[i];
			fileDump << "\t" << cmpmat[i].getSegmentCount();
			fileDump << "\t" << cmpmat[i].getInteractionCount();
			fileDump << "\t" << chromDeltaCounts[i];
			fileDump << endl;
		}
		fileDump << "GRAND TOTAL";
		fileDump << "\t" << totalSegments;
		fileDump << "\t" << totalInteractions;
		fileDump << "\t" << totalDeltas;
		fileDump << endl;
		
		fileDump.close();
	}
}

void getChromosomesDataSizeDesc(vector<string>& chromosomes, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, vector<pair<ChromosomeIndexType, uint64_t > >& chromosomesRows)
{
	chromosomesRows.clear();
	chromosomesRows.reserve(chromosomes.size());
	for(ChromosomeIndexType i=0; i<chromosomes.size(); ++i) {
		chromosomesRows.push_back(make_pair(i, cmpmat[i].getSegmentCount()));
	}
	sort(chromosomesRows.begin(), chromosomesRows.end(), chromSegmentCountDesc);
}

void sweepQunatilesDeltaSpline(string& filename, vector<string>& chromosomes, vector<pair<ChromosomeIndexType, uint64_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, vector<uint32_t >& percentilesList, uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/, bool masking/*=false*/)
{
	double t_diff;
	
	// Get all deltas for all chromosomes:
	vector<map<ChromosomeIndexType, DeltaStatistics > > chromDeltasStats;
	chromDeltasStats.resize(chromosomes.size());
	
	std::sort(percentilesList.begin(), percentilesList.end(), std::greater<unsigned int>());
	
	// deltas allocations
	uint32_t nMaxSpan = 0;
	{
		for(ChromosomeIndexType i=0; i<chromosomes.size(); ++i) {
			uint32_t maxChromSpan = cmpmat[i].getWidestInteractionSpan();
			if (maxChromSpan>nMaxSpan) nMaxSpan = maxChromSpan;
		}
		nMaxSpan++; // 0-based index, so need 1 more allocation
		//cerr << "Largest span " << nMaxSpan << endl;
	}
	uint32_t* pDeltas = (uint32_t*) calloc(nMaxSpan+1, sizeof(uint32_t));
	BIT_ARRAY *deltaPresences = barcreate(nMaxSpan+1);
	
	vector<vector<precomputeType > > cells;
	cells.resize(percentilesList.size());
	
	//int dmax = (maxDelta == MAXDELTA) ? deltas.back() / 100 : maxDelta;
	uint32_t dmax = (maxDelta == MAXDELTA) ? nMaxSpan / 100 : maxDelta;
	vector<double> d = getSequence((double)minDelta, (double)dmax, 2000);
	for(unsigned int percentileIdx = 0; percentileIdx<percentilesList.size(); ++percentileIdx) {
		cells[percentileIdx].resize(d.size());
	}
	
	
	double t_start = realtime();
	uint64_t ulTotalDeltas = 0;
	for(ChromosomeIndexType i=0; i<chromosomes.size(); ++i) {
		double t_subStart = realtime();
		ChromosomeIndexType nChrIdx = chromosomesRows[i].first;
		uint64_t ulDeltas = cmpmat[nChrIdx].getDeltasCount(pDeltas, deltaPresences, nMaxSpan+1, minDelta, maxDelta, masking);
		ulTotalDeltas += ulDeltas;
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] #deltas(%s)=%llu collected, %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), ulDeltas, t_diff/60.0);
	}
	t_diff = realtime() - t_start;
	fprintf(stderr, "[M::%s:%d] #deltas(genome)=%llu collected, %.2f min..\n", __func__, __LINE__, ulTotalDeltas, t_diff/60.0);

	for(unsigned int percentileIdx = 0; percentileIdx<percentilesList.size(); ++percentileIdx) {
		double t_subStart = realtime();
		QuantileMapper qtm;
		spline s;
		qtm.computeDeltaCountQuantile(ulTotalDeltas, pDeltas, percentilesList[percentileIdx], NULL);
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] delta quantiles(%u) computation deltas took %.2f min..\n", __func__, __LINE__, percentilesList[percentileIdx], t_diff/60.0);
		
		t_subStart = realtime();
		map<uint32_t, DeltaStatistics > deltasStats;
		{
			vector<uint32_t>& quantileDeltas = qtm.getQuantileDeltas();
			for(vector<uint32_t>::iterator it=quantileDeltas.begin(); it!=quantileDeltas.end(); ++it) {
				uint32_t deltaQuantile = qtm.getDeltaQuantile(*it);
				deltasStats[deltaQuantile] = DeltaStatistics();
			}
		}
		for(ChromosomeIndexType i=0; i<chromosomesRows.size(); i++) {
			double t_subStart = realtime();
			int nChrIdx = chromosomesRows[i].first;
#if 0
			//
			// TODO: IMPORTANT : FIX AFTER CHECKING THE SINGLE SWEEP OPERATION ABOVE
			//
			cmpmat[nChrIdx].getStatisticsPerDelta(deltasStats, &qtm, minDelta, maxDelta, masking);
#endif
			t_diff = realtime() - t_subStart;
			fprintf(stderr, "[M::%s:%d] delta quantiles statistics computation (%s), %.2f min..\n", __func__, __LINE__, chromosomes[nChrIdx].c_str(), t_diff/60.0);
		}
		t_diff = realtime() - t_subStart;
		fprintf(stderr, "[M::%s:%d] delta quantiles(%u) statistics computation deltas took %.2f min..\n", __func__, __LINE__, percentilesList[percentileIdx], t_diff/60.0);
		
		vector<double> delta_x; // x-object for spline generation
		vector<double> exp_y; // y-object for spline generation
		for(map<uint32_t, DeltaStatistics>::iterator iter = deltasStats.begin(); iter != deltasStats.end(); iter++) {
			delta_x.push_back((double)(iter->first));
			exp_y.push_back((double)(iter->second.getMean()));
			//cout << iter->first << " " << iter->second << endl;
		}
		
		// If x is not sorted, an error will occur.
		s.set_points(delta_x,exp_y);    // currently it is required that X is already sorted

		// keep results
		for(unsigned int i=0; i < d.size(); i++) {
			cells[percentileIdx][i] = s.at(d[i]) ;
		}
	}
	
	free(pDeltas);
	bardestroy(deltaPresences);
	
	// write result in a table grid
	if (percentilesList.size()>0)
	{
		string ofile(filename); ofile.append(".quantiles.sweep.xls");
		ofstream fileDump;
		fileDump.open(ofile.c_str());
		
		// print header
		fileDump << "delta";
		for(unsigned int percentileIdx = 0; percentileIdx<percentilesList.size(); ++percentileIdx) {
			fileDump << "\tsmoothed-" << percentilesList[percentileIdx] ;
		}
		fileDump << endl;

		for(unsigned int i=0; i < d.size(); i++) {
			fileDump << d[i];
			for(unsigned int percentileIdx = 0; percentileIdx<percentilesList.size(); ++percentileIdx) {
				fileDump << "\t" << cells[percentileIdx][i];
			}
			fileDump << endl;
		}
		fileDump.close();
	}
}


void computeInteractionsPvalues(vector<string>& chromosomes, vector<pair<ChromosomeIndexType, uint64_t > >&chromosomesRows, map<ChromosomeIndexType, CCCMatrixInteraction>& cmpmat, spline& s, uint32_t deltaCutoff, uint32_t minDelta/*=8000*/, uint32_t maxDelta/*=MAXDELTA*/)
{
	// [exp{chromsome}, pval{chromosome}] is faster than [chromosome{exp, pval}]
	for(ChromosomeIndexType j=0; j<chromosomesRows.size(); j++) {
		int nChrIdx = chromosomesRows[j].first;
		string chr=chromosomes[nChrIdx];
		double t_loopStart = realtime();
		cmpmat[nChrIdx].setDeltaToExpectationMapper(s, minDelta, MAXDELTA);
		double t_diff = realtime() - t_loopStart;
		fprintf(stderr, "[M::%s:%d] createExpectationMatrix %s N:%u L:%g rtmin:%.2f\n", __func__, __LINE__, chr.c_str(), cmpmat[nChrIdx].getN(), cmpmat[nChrIdx].getExpectN(), t_diff/60.0);
	}
	for(ChromosomeIndexType j=0; j<chromosomesRows.size(); j++) {
		int nChrIdx = chromosomesRows[j].first;
		string chr=chromosomes[nChrIdx];
		double t_loopStart = realtime();
		cmpmat[nChrIdx].calculatePvalues(deltaCutoff, MAXDELTA);// During refinement, all P-values for all deltas are considered!
		
		double t_diff = realtime() - t_loopStart;
		fprintf(stderr, "[M::%s:%d] calculatePvalues %s rtmin:%.2f\n", __func__, __LINE__, chr.c_str(), t_diff/60.0);
	}
}


// Add the possible template classes explicitly, to allow linker to find them:
// These are the only allowed:
template vector<double> getSequence(double from, double to, uint32_t length);
template vector<uint32_t> getSequence(uint32_t from, uint32_t to, uint32_t length);
//template vector<float> getSequence(float from, float to, uint32_t length);
