// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen

#define PROGRAM_NAME "ChiaSig"
#define PROGRAM_VERSION "v1.19.45"

#include "../CCCsig/CCCDataReader.h"
#include "../CCCsig/Segment.h"
#include "../CCCsig/CCCMatrix.h"
#include "../CCCsig/CCCStatistics.h"
#include "../stocc/stocc.h" // Non-central hypergeometric and Poisson

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "../tclap-1.2.1/include/tclap/CmdLine.h" // Command-line parsing

using namespace std;

// WCH:
#include <sys/time.h>
#include <sys/resource.h>


#define STREAMDELTAS 1

const double dOmegaPrecision = 1E-50;

double G_t_real;

#ifdef __cplusplus
extern "C" {
#endif
	double cputime();
	double realtime();
#ifdef __cplusplus
}
#endif


/*********
 * Timer *
 *********/

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
// END - WCH:

#define DEFAULT_NQUANT 500

void print_spline(ostream& ostr, spline s1, spline s2, vector<double> deltas, unsigned int nQuant=DEFAULT_NQUANT) {
	ostr << "#delta" << "\t" << "smoothed-" << nQuant << "\t" << "refined-" << nQuant << endl;
	for(unsigned int i=0; i < deltas.size(); i++) {
		ostr << deltas[i] << "\t" << s1.cachedAt(deltas[i]) << "\t" << s2.cachedAt(deltas[i]) << endl;
	}
}

void debug_report_Matrices(vector<string>& chromosomes, map<ChromosomeIndexType, CCCMatrixInteraction>& contactMatrices, const char* szStage, bool hideMask=false) {
	for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
		string chr=chromosomes[i];
#ifdef DUMP_CHROMOSOME
		if (chr!=DUMP_CHROMOSOME) continue;
#endif
		ofstream fileDump;
		string filename = chr+".cm."+szStage+".sparsed";
		fileDump.open(filename.c_str());
		fileDump << "=== " << chr.c_str() << " START refined contact matrix ===" << endl;
		contactMatrices[i].printInteractionMatrix(fileDump, hideMask);
		fileDump << "=== " << chr.c_str() << " END refined contact matrix ===" << endl;
		fileDump.close();
	}
	for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
		string chr=chromosomes[i];
#ifdef DUMP_CHROMOSOME
		if (chr!=DUMP_CHROMOSOME) continue;
#endif
		ofstream fileDump;
		string filename = chr+".em."+szStage+".sparsed";
		fileDump.open(filename.c_str());
		fileDump << "=== " << chr.c_str() << " START refined expectation matrix ===" << endl;
		contactMatrices[i].printExpectationMatrix(fileDump, hideMask);
		fileDump << "=== " << chr.c_str() << " END refined expectation matrix ===" << endl;
		fileDump.close();
	}
	for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
		string chr=chromosomes[i];
#ifdef DUMP_CHROMOSOME
		if (chr!=DUMP_CHROMOSOME) continue;
#endif
		ofstream fileDump;
		string filename = chr+".pm."+szStage+".sparsed";
		fileDump.open(filename.c_str());
		fileDump << "=== " << chr.c_str() << " START refined pvalue matrix ===" << endl;
		contactMatrices[i].printPvalueMatrix(fileDump, hideMask);
		fileDump << "=== " << chr.c_str() << " END refined pvalue matrix ===" << endl;
		fileDump.close();
	}
}
	
int main(int argc, char** argv) {
	G_t_real = realtime();
	double t_diff;
	
	try {
		TCLAP::CmdLine cmd("ChiaSig is a program to find significant interactions in ChIA-PET datasets, using the Non-central Hypergeometric (NCHG) distribution. For details about this statistical method, see Paulsen et al. 'A statistical model of ChIA-PET data for accurate detection of chromatin 3D interactions', Nucleic acids research. (2014).\n Copyright (2014) Jonas Paulsen" , ' ', "0.93");
		
		TCLAP::SwitchArg dumpMatricesSwitch("","dump","Save both initial and refined matrices for debugging", false);
		cmd.add( dumpMatricesSwitch );
		
		TCLAP::SwitchArg dumpInitialMatricesSwitch("","dumpinitial","Save initial matrices for debugging", false);
		cmd.add( dumpInitialMatricesSwitch );
		
		TCLAP::SwitchArg dumpRefinedMatricesSwitch("","dumprefined","Save refined matrices for debugging", false);
		cmd.add( dumpRefinedMatricesSwitch );
		
		TCLAP::SwitchArg dumpSplineSwitch("","dumpspline","Save spline details for debugging", false);
		cmd.add( dumpSplineSwitch );
		
		TCLAP::SwitchArg resourceEstimationArg("","resource","Only report essential resource metrices.", false);
		cmd.add( resourceEstimationArg);
		
		TCLAP::MultiArg<unsigned int> quantilesSweepArg("","sweep","Specify (multiple) number of quantiles to use for estimating the expected number of interactions given the genomic distance. (See -n option.) This output a tabulated list of results for comparison to select the best quantiles.", false, "int");
		cmd.add( quantilesSweepArg);
		
		TCLAP::SwitchArg noEarlyTerminationSwitch("","noearly","Do NOT terminate estimation loop early", false);
		cmd.add( noEarlyTerminationSwitch );
		
		TCLAP::ValueArg<unsigned int> percentilesArg("n","quantiles","Number of quantiles to use for estimating the expected number of interactions given the genomic distance (default is 500). Higher numbers (1000 and more) are appropriate for very large datasets, with many interacting anchors, while lower numbers (some hundred) will be appropriate for smaller datasets.",false,DEFAULT_NQUANT,"int");
		cmd.add( percentilesArg );
		
		TCLAP::ValueArg<double> depthThresholdArg("","depthThreshold","Depth threshold to switch to span scanning algorithm (default is 0.25).",false,0.25,"float");
		cmd.add( depthThresholdArg );
		
		TCLAP::SwitchArg skipFDRArg("f","skipFDR","Skip FDR correction, report raw P-values only",false);
		cmd.add( skipFDRArg );
		
		
		TCLAP::ValueArg<double> alphaArg("a","alpha","False discovery rate for selecting significant interactions (default is 0.05).",false,0.05,"double");
		cmd.add( alphaArg );
		
		
		TCLAP::ValueArg<uint32_t> cutoffArg("c","cutoff","Minimum allowed number of interactions for a given pair of anchors, to consider it significant (default is 3)",false,3,"uint32_t");
		cmd.add( cutoffArg );
		
		TCLAP::ValueArg<double> refinealphaArg("A","refinealpha","False discovery rate used for refinement (default is 0.01).",false,0.01,"double");
		cmd.add( refinealphaArg );
		
		
		TCLAP::ValueArg<uint32_t> refinecutoffArg("C","refinecutoff","Minimum allowed number of interactions for a given pair of anchors, to consider it significant during refinement (default is 3)",false,3,"uint32_t");
		cmd.add( refinecutoffArg );
		
		TCLAP::ValueArg<uint32_t> stepsArg("s","steps","Number of steps used to calculate false discovery rate (default is 22). Higher numbers will give more accurate estimation of the false discovery rate, but will also be slower.",false,22,"uint32_t");
		cmd.add( stepsArg );
		
		TCLAP::ValueArg<uint32_t> refinestepsArg("S","refinesteps","Number of steps used to calculate false discovery rate, during refinement (default is 22). Higher numbers will give more accurate estimation of the false discovery rate, but will also be slower.",false,22,"uint32_t");
		cmd.add( refinestepsArg );
		
		TCLAP::SwitchArg printdeltaSwitch("d","printdelta","Print estimated expectation values to stderr", false);
		cmd.add( printdeltaSwitch );
		
		
		TCLAP::SwitchArg printInteractionCountsSwitch("p","printcounts","Print observed and expected number of interactions for each significant interactions, in addition to P-values", false);
		cmd.add( printInteractionCountsSwitch );
		TCLAP::SwitchArg printAllInteractionsSwitch("P","printall","Print all interactions regardless of significance", false);
		cmd.add( printAllInteractionsSwitch );
		
		
		TCLAP::SwitchArg onlyprintdeltaSwitch("o","onlyprintdelta","Estimate the expectation values one time (no refinement), then print estimated expectation values to stdout and exit. No P-values are calculated. Suitable for exploring the effect of choosing different number of quantiles (-n).", false);
		cmd.add( onlyprintdeltaSwitch );
		
		//TCLAP::SwitchArg noRefinementSwitch("r","norefinement","Skip the refinement step when estimating expectation values. When using this option, expectation values are estimated from the raw (non-refined) data.", false);
		//cmd.add( onlyprintdeltaSwitch );

		
		TCLAP::ValueArg<uint32_t> mindeltaArg("m","mindelta","Minumum genomic distance (in bp) allowed between anchor pairs, below which interactions are excluded. (default is 8000 bp)",false,8000,"uint32_t");
		cmd.add( mindeltaArg );
		
		
		TCLAP::ValueArg<uint32_t> maxdeltaArg("M","maxdelta","Maximum genomic distance (in bp) allowed between anchor pairs, above which interactions are excluded. (Default is no maximum value)",false,MAXDELTA,"uint32_t");
		cmd.add( maxdeltaArg );
		
		TCLAP::ValueArg<string> selectChrArg("u","onlyChr","Use only the specified chromosome for P-value and FDR calculations. Default behaviour is to use all chromosomes. Even if a chromosome is specified using this argument, all chromosomes will be used for estimation of expected number of interactions.",false,"","string");
		cmd.add( selectChrArg );
		
		
		TCLAP::ValueArg<int> threadArg("t","thread","Number of threads (default is 1)",false,1,"int");
		cmd.add( threadArg );
		
		
		TCLAP::ValueArg<string> outputPrefixArg("","output","Output prefix.",false,"","string");
		cmd.add( outputPrefixArg );
		
		TCLAP::UnlabeledValueArg<string> nolabel( "filename", "Input file of the format chrA startA endA chrB startB endB I(A,B) where I(A,B) gives the number of interactions between A and B. This corresponds to the BEDPE format described here: http://bedtools.readthedocs.org/en/latest/content/general-usage.html#bedpe-format", true, "/dev/null", "filename"  );
		cmd.add( nolabel );
		
		cmd.parse( argc, argv );
		
		
		// WCH : to facilitate comparison, record the command invocation
		cerr << "# commands:"; for (int i=0; i<argc; ++i) cerr << " " << argv[i]; cerr << endl;

		
		bool resourceEstimation = resourceEstimationArg.getValue();
		vector<unsigned int> quantilesSweep = quantilesSweepArg.getValue();
		bool noEarlyTermination = noEarlyTerminationSwitch.getValue();
		
		bool dumpMatrices = dumpMatricesSwitch.getValue();
		bool dumpInitialMatrices = dumpInitialMatricesSwitch.getValue();
		bool dumpRefinedMatrices = dumpRefinedMatricesSwitch.getValue();
		bool dumpSplines = dumpSplineSwitch.getValue();
		
		unsigned int percentiles = percentilesArg.getValue();
		float depthThreshold = depthThresholdArg.getValue();
		uint32_t minDelta = mindeltaArg.getValue();
		uint32_t maxDelta = maxdeltaArg.getValue();
		bool printDelta = printdeltaSwitch.getValue();
		bool onlyprintDelta = onlyprintdeltaSwitch.getValue();
		bool printAllInteractions = printAllInteractionsSwitch.getValue();
		string fileName = nolabel.getValue();
		double alphaCut = alphaArg.getValue();
		uint32_t cutoff = cutoffArg.getValue();
		int nThreads = threadArg.getValue();
		
		uint32_t steps= stepsArg.getValue();
		uint32_t refinesteps= refinestepsArg.getValue();
		
		double refinealphaCut = refinealphaArg.getValue();
		uint32_t refinecutoff = refinecutoffArg.getValue();
		
		bool skipFDR = skipFDRArg.getValue();
		bool printNij = printInteractionCountsSwitch.getValue();
		//bool noRefinement = noRefinementSwitch.getValue();
		string selectChr = selectChrArg.getValue();
		string outputPrefix = outputPrefixArg.getValue();
		if (outputPrefix.empty()) {
			outputPrefix.append(fileName);
		}
		
		DeltaCountSums deltaCountSums;
		
		assert(alphaCut > 0 and alphaCut <= 1);
		assert(refinealphaCut > 0 and refinealphaCut <= 1);
		
		
		setThreads(nThreads);
		
		// WCH : to facilitate comparison, record the command invocation
		string outputSigfInteraction(outputPrefix); outputSigfInteraction.append(".sigf.interactions");
		ofstream fileSigfInteraction;
		fileSigfInteraction.open(outputSigfInteraction.c_str());
		fileSigfInteraction << "# version: " << PROGRAM_NAME << " " << PROGRAM_VERSION << endl;
		fileSigfInteraction << "# commands:"; for (int i=0; i<argc; ++i) fileSigfInteraction << " " << argv[i]; fileSigfInteraction << endl;
		fileSigfInteraction.close();
		
		// READING THE DATA:
		double t_start = realtime();
		uint32_t deltaCutoff=minDelta;
		map<ChromosomeIndexType, CCCMatrixInteraction> contactMatrices;
		vector<string> chromosomes;
		vector<pair<ChromosomeIndexType, uint64_t > > sizeOrderedChromosomes;
		{
			CCCDataReader dr(fileName, contactMatrices);
			dr.buildContactMatrices();
			dr.getChromosomes(chromosomes);
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] dr.buildContactMatrices() took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			//cout << "Data reading done! Read: " << chromosomes.size() << " chromosomes" << endl;
			
			t_start = realtime();
			// Build contact-matrices:
			assert(0==selectChr.size() || dr.isChromosomePresent(selectChr));
			for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
				string chr=chromosomes[i];
				//contactMatrices[i] = dr.getContactMatrix(i);
				contactMatrices[i].setThreads(nThreads);
				contactMatrices[i].setDepthThreshold(depthThreshold);
				contactMatrices[i].InteractionLoaded(minDelta, maxDelta);
				//contactMatrices[i].InteractionLoaded();
				//cout << "chr read: " << chr << " Num. interactions: " << contactMatrices[chr].getN() << endl;
			}
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] Build contact-matrices took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			
			// keep the chromosomes in descending order of data size
			getChromosomesDataSizeDesc(chromosomes, contactMatrices, sizeOrderedChromosomes);
		}
		
		if (resourceEstimation) {
			t_start = realtime();
			fprintf(stderr, "[M::%s:%d] Performing resource estimation..\n", __func__, __LINE__);
			estimateResources(fileName, chromosomes, sizeOrderedChromosomes, contactMatrices, deltaCutoff, MAXDELTA, false);
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] Resource estimation took %.2f min..\n", __func__, __LINE__, t_diff/60.0);

			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			return 0;
		}
		

		if (quantilesSweep.size()>0) {
			t_start = realtime();
			fprintf(stderr, "[M::%s:%d] Performing %lu quantiles sweep..\n", __func__, __LINE__, quantilesSweep.size());
			sweepQunatilesDeltaSpline(fileName, chromosomes, sizeOrderedChromosomes, contactMatrices, quantilesSweep, deltaCutoff, MAXDELTA, false);
			t_diff = realtime() - t_start;
			fprintf(stderr, "[M::%s:%d] %lu quantiles sweep took %.2f min..\n", __func__, __LINE__, quantilesSweep.size(), t_diff/60.0);
			
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			return 0;
		}
		
		// Estimate the spline-object, for the first time:
		t_start = realtime();
		spline s; s.setThreads(nThreads);
		vector<double> deltas;
		string outputFirstDeltaSpline(outputPrefix); outputFirstDeltaSpline.append(".first.deltaspline.chiasig");
		estimateDeltaSpline(outputFirstDeltaSpline.c_str(), s, deltas, deltaCountSums, chromosomes, sizeOrderedChromosomes, contactMatrices, deltaCutoff, MAXDELTA, false, percentiles, dumpSplines); // Spline is estimated on all delta>minDelta
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] estimateDeltaSpline took %.2f min..\n", __func__, __LINE__, t_diff/60.0);

		
		if(onlyprintDelta) {
			uint32_t dmax = (maxDelta == MAXDELTA) ? deltas.back() / 100 : maxDelta;
			vector<double> d = getSequence((double)deltaCutoff, (double)dmax, 2000);
			//print_spline(cout, s, s, d, percentiles);
			fileSigfInteraction.open(outputSigfInteraction.c_str(), ofstream::out | ofstream::app);
			print_spline(fileSigfInteraction, s, s, d, percentiles);
			fileSigfInteraction.close();
			
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			return 0;
		}

  
		t_start = realtime();
		CFishersNCHypergeometric::setAccuracy(dOmegaPrecision);
		// Calculating P-values (first time):
		computeInteractionsPvalues(chromosomes, sizeOrderedChromosomes, contactMatrices, s, deltaCutoff, minDelta, MAXDELTA);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] p-values calculation first time took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		
		
		if (dumpMatrices || dumpInitialMatrices)
			debug_report_Matrices(chromosomes, contactMatrices, "initial", false);

		
		// IF ONLY RAW P-VALUES SHOULD BE PRINTED::
		if(skipFDR) {
			fileSigfInteraction.open(outputSigfInteraction.c_str(), ofstream::out | ofstream::app);
			for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
				string chr=chromosomes[i];
				printPositives(fileSigfInteraction, chr, contactMatrices[i], 1.0, 0,printNij,printAllInteractions); // Everything is printed. Including 0s.
			}
			fileSigfInteraction.close();
			
			t_diff = realtime() - G_t_real;
			fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
			return 0;
		}
		
		
		// ----
		
		if (selectChr != "") {
			vector<string> selChr;
			selChr.push_back(selectChr);
			chromosomes = selChr;
		}
		
		
		t_start = realtime();
		// Masking, based on P-values, using the double procedure:
		// In original code, each "steps" will iterate 6 alpha values or 5 intervals
		string outputFirstAlpha(outputPrefix); outputFirstAlpha.append(".first.alpha.chiasig");
		map<double, double> FDRs = estimateFDRAlpha(outputFirstAlpha.c_str(), chromosomes, sizeOrderedChromosomes, contactMatrices, refinealphaCut, refinesteps, noEarlyTermination, refinecutoff, deltaCutoff, maxDelta);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] Masking based on p-values, estimateFDRAlpha took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		t_start = realtime();
		pair<double, double> alphaRange = getAlphaRange(FDRs, refinealphaCut);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] getAlphaRange(alpha:%.6e, FDR:%.6e) took %.2f min..\n", __func__, __LINE__, alphaRange.first, FDRs[alphaRange.first], t_diff/60.0);
		
		cerr << "# Masking P-value cutoff: " << scientific << setprecision(6) << alphaRange.first << ". FDR<=" << scientific << setprecision(6) << FDRs[alphaRange.first] << endl;
		fileSigfInteraction.open(outputSigfInteraction.c_str(), ofstream::out | ofstream::app);
		fileSigfInteraction << "# Masking P-value cutoff: " << scientific << setprecision(6) << alphaRange.first << ". FDR<=" << scientific << setprecision(6) << FDRs[alphaRange.first] << endl;
		fileSigfInteraction.close();

		
		// Masking data at the selected range:
		t_start = realtime();
		unsigned long genomeMasked = 0;
		for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
			double t_loopStart = realtime();
			string chr=chromosomes[i];
			//unsigned long masked = contactMatrices[i].maskByAlpha(alphaRange.first, cutoff, ("chr18" == chr));
			unsigned long masked = contactMatrices[i].maskByAlpha(alphaRange.first, cutoff, true);
			genomeMasked+=masked;
			t_diff = realtime() - t_loopStart;
			fprintf(stderr, "[M::%s:%d] #masked(%s)=%lu masked in maskByAlpha took %.2f min..\n", __func__, __LINE__, chr.c_str(), masked, t_diff/60.0);
		}
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] #masked(genome)=%lu took %.2f min..\n", __func__, __LINE__, genomeMasked, t_diff/60.0);
		
		t_start = realtime();
		//fdr: Re-estimating the delta-function:
		spline sRefined; sRefined.setThreads(nThreads);
		vector<double> deltasRefined;
		//string outputSecondDeltaSpline(outputPrefix); outputSecondDeltaSpline.append(".second.deltaspline.chiasig");
		string outputSecondDeltaSpline(outputFirstDeltaSpline);
		estimateDeltaSpline(outputSecondDeltaSpline.c_str(), sRefined, deltasRefined, deltaCountSums, chromosomes, sizeOrderedChromosomes, contactMatrices, deltaCutoff, MAXDELTA, true, percentiles, dumpSplines);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] estimateDeltaSpline took %.2f min..\n", __func__, __LINE__, t_diff/60.0);

		
		// Re-estimating expectations and P-values, now with updated expectations/lambdas:
		t_start = realtime();
		computeInteractionsPvalues(chromosomes, sizeOrderedChromosomes, contactMatrices, sRefined, deltaCutoff, minDelta, MAXDELTA);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] Re-estimating expectations and P-values, now with updated expectations/lambdas, loop took %.2f min..\n", __func__, __LINE__, t_diff/60.0);

		
		if (dumpMatrices || dumpRefinedMatrices)
			debug_report_Matrices(chromosomes, contactMatrices, "refined", true);

		
		// Re-estimating FDR on new P-values:
		t_start = realtime();
		string outputSecondAlpha(outputPrefix); outputSecondAlpha.append(".second.alpha.chiasig");
		// In original code, each "steps" will iterate 6 alpha values or 5 intervals
		FDRs = estimateFDRAlpha(outputSecondAlpha.c_str(), chromosomes, sizeOrderedChromosomes, contactMatrices, alphaCut, steps, noEarlyTermination, cutoff, deltaCutoff, maxDelta);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] Re-estimating FDR on new P-values, estimateFDRAlpha took %.2f min..\n", __func__, __LINE__, t_diff/60.0);
		t_start = realtime();
		alphaRange = getAlphaRange(FDRs, alphaCut);
		t_diff = realtime() - t_start;
		fprintf(stderr, "[M::%s:%d] Re-estimating FDR on new P-values, getAlphaRange(alpha:%.6e, FDR:%.6e) took %.2f min..\n", __func__, __LINE__, alphaRange.first, FDRs[alphaRange.first], t_diff/60.0);
		
		
		cerr << "# Refined P-value cutoff: " << scientific << setprecision(6) << alphaRange.first << ". FDR<=" << scientific << setprecision(6) << FDRs[alphaRange.first] << endl;
		fileSigfInteraction.open(outputSigfInteraction.c_str(), ofstream::out | ofstream::app);
		fileSigfInteraction << "# Refined P-value cutoff: " << scientific << setprecision(6) << alphaRange.first << ". FDR<=" << scientific << setprecision(6) << FDRs[alphaRange.first] << endl;
		
		//Print the final, significant interactions:
		for(ChromosomeIndexType i=0; i<chromosomes.size(); i++) {
			string chr=chromosomes[i];
			printPositives(fileSigfInteraction, chr, contactMatrices[i], alphaRange.first, cutoff, printNij,printAllInteractions);
		}
		
		// Print the estimated expectation values:
		if(printDelta) print_spline(cerr, s, sRefined, deltas, percentiles);
		
		fileSigfInteraction.close();
	}
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	
	t_diff = realtime() - G_t_real;
	fprintf(stderr, "[M::%s:%d] Total time %.2f min..\n", __func__, __LINE__, t_diff/60.0);
	
	return 0;
}

