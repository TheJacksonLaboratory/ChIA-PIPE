/*****************************   stocc.h   **********************************
 * Author:        Agner Fog
 * Date created:  2004-01-08
 * Last modified: 2013-09-20
 * Project:       randomc.h
 * Source URL:    www.agner.org/random
 *
 * Description:
 * This file contains function prototypes and class declarations for the C++
 * library of non-uniform random number generators. Most functions are fast and
 * accurate, even for extreme values of the parameters.
 *
 *
 * functions without classes:
 * ==========================
 *
 * void EndOfProgram(void);
 * System-specific exit code. You may modify this to make it fit your
 * user interface.
 *
 * void FatalError(const char * ErrorText);
 * Used for outputting error messages from the other functions and classes.
 * You may have to modify this function to make it fit your user interface.
 *
 * double Erf (double x);
 * Calculates the error function, which is the integral of the normal distribution.
 *
 * double LnFac(int32_t n);
 * Calculates the natural logarithm of the factorial of n.
 *
 *
 * class StochasticLib1:
 * ====================
 * This class can be derived from any of the uniform random number generators
 * defined in randomc.h. StochasticLib1 provides the following non-uniform random
 * variate generators:
 *
 * int Bernoulli(double p);
 * Bernoulli distribution. Gives 0 or 1 with probability 1-p and p.
 *
 * double Normal(double m, double s);
 * Normal distribution with mean m and standard deviation s.
 *
 * double NormalTrunc(double m, double s, double limit);
 * Truncated normal distribution with tails cut off at m +/- limit
 *
 * int32_t Poisson (double L);
 * Poisson distribution with mean L.
 *
 * int32_t Binomial (int32_t n, double p);
 * Binomial distribution. n trials with probability p.
 *
 * int32_t Hypergeometric (int32_t n, int32_t m, int32_t N);
 * Hypergeometric distribution. Taking n items out N, m of which are colored.
 *
 * void Multinomial (int32_t * destination, double * source, int32_t n, int colors);
 * void Multinomial (int32_t * destination, int32_t * source, int32_t n, int colors);
 * Multivariate binomial distribution.
 *
 * void MultiHypergeometric (int32_t * destination, int32_t * source, int32_t n, int colors);
 * Multivariate hypergeometric distribution.
 *
 * void Shuffle(int * list, int min, int n);
 * Shuffle a list of integers.
 *
 *
 * class StochasticLib2:
 * =====================
 * This class is derived from class StochasticLib1. It redefines the functions
 * Poisson, Binomial and HyperGeometric.
 * In StochasticLib1, these functions are optimized for being called with
 * parameters that vary. In StochasticLib2, the same functions are optimized
 * for being called repeatedly with the same parameters. If your parameters
 * seldom vary, then StochasticLib2 is faster. The two classes use different
 * calculation methods, both of which are accurate.
 *
 *
 * class StochasticLib3:
 * =====================
 * This class can be derived from either StochasticLib1 or StochasticLib2,
 * whichever is preferred. It contains functions for generating variates with
 * the univariate and multivariate Wallenius' and Fisher's noncentral
 * hypergeometric distributions.
 *
 * int32_t WalleniusNCHyp (int32_t n, int32_t m, int32_t N, double odds);
 * Sampling from Wallenius' noncentral hypergeometric distribution, which is
 * what you get when taking n items out N, m of which are colored, without
 * replacement, with bias.
 *
 * int32_t FishersNCHyp (int32_t n, int32_t m, int32_t N, double odds);
 * Sampling from Fisher's noncentral hypergeometric distribution which is the
 * conditional distribution of independent binomial variates given their sum n.
 *
 * void MultiWalleniusNCHyp (int32_t * destination, int32_t * source, double * weights, int32_t n, int colors);
 * Sampling from multivariate Wallenius' noncentral hypergeometric distribution.
 *
 * void MultiFishersNCHyp (int32_t * destination, int32_t * source, double * weights, int32_t n, int colors);
 * Sampling from multivariate Fisher's noncentral hypergeometric distribution.
 *
 *
 * Uniform random number generators (integer and float) are also available, as
 * these are inherited from the random number generator class that is the base
 * class of StochasticLib1.
 *
 *
 * class CWalleniusNCHypergeometric
 * ================================
 * This class implements various methods for calculating the probability
 * function and the mean and variance of the univariate Wallenius' noncentral
 * hypergeometric distribution. It is used by StochasticLib3 and can also be
 * used independently.
 *
 *
 * class CMultiWalleniusNCHypergeometric
 * =====================================
 * This class implements various methods for calculating the probability func-
 * tion and the mean of the multivariate Wallenius' noncentral hypergeometric
 * distribution. It is used by StochasticLib3 and can also be used independently.
 *
 *
 * class CMultiWalleniusNCHypergeometricMoments
 * ============================================
 * This class calculates the exact mean and variance of the multivariate
 * Wallenius' noncentral hypergeometric probability distribution.
 *
 *
 * class CFishersNCHypergeometric
 * ==============================
 * This class calculates the probability function and the mean and variance
 * of Fisher's noncentral hypergeometric distribution.
 *
 *
 * class CMultiFishersNCHypergeometric
 * ===================================
 * This class calculates the probability function and the mean and variance
 * of the multivariate Fisher's noncentral hypergeometric distribution.
 *
 *
 * source code:
 * ============
 * The code for EndOfProgram and FatalError is found in the file userintf.cpp.
 * The code for the functions in StochasticLib1 is found in the file stoc1.cpp.
 * The code for the functions in StochasticLib2 is found in the file stoc2.cpp.
 * The code for the functions in StochasticLib3 is found in the file stoc3.cpp.
 * The code for the functions in CWalleniusNCHypergeometric,
 * CMultiWalleniusNCHypergeometric and CMultiWalleniusNCHypergeometricMoments
 * is found in the file wnchyppr.cpp.
 * The code for the functions in CFishersNCHypergeometric and
 * CMultiFishersNCHypergeometric is found in the file fnchyppr.cpp
 * LnFac is found in stoc1.cpp.
 * Erf is found in wnchyppr.cpp.
 *
 *
 * Examples:
 * =========
 * The file ex-stoc.cpp contains an example of how to use this class library.
 *
 * The file ex-cards.cpp contains an example of how to shuffle a list of items.
 *
 * The file ex-lotto.cpp contains an example of how to generate a sequence of
 * random integers where no number can occur more than once.
 *
 * The file testbino.cpp contains an example of sampling from the binomial distribution.
 *
 * The file testhype.cpp contains an example of sampling from the hypergeometric distribution.
 *
 * The file testpois.cpp contains an example of sampling from the poisson distribution.
 *
 * The file testwnch.cpp contains an example of sampling from Wallenius noncentral hypergeometric distribution.
 *
 * The file testfnch.cpp contains an example of sampling from Fisher's noncentral hypergeometric distribution.
 *
 * The file testmwnc.cpp contains an example of sampling from the multivariate Wallenius noncentral hypergeometric distribution.
 *
 * The file testmfnc.cpp contains an example of sampling from the multivariate Fisher's noncentral hypergeometric distribution.
 *
 * The file evolc.zip contains examples of how to simulate biological evolution using this class library.
 *
 *
 * Documentation:
 * ==============
 * The file ran-instructions.pdf contains further documentation and
 * instructions for these random number generators.
 *
 * The file distrib.pdf contains definitions of the standard statistic distributions:
 * Bernoulli, Normal, Poisson, Binomial, Hypergeometric, Multinomial, MultiHypergeometric.
 *
 * The file sampmet.pdf contains theoretical descriptions of the methods used
 * for sampling from these distributions.
 *
 * The file nchyp.pdf, available from www.agner.org/random/, contains
 * definitions of the univariate and multivariate Wallenius and Fisher's
 * noncentral hypergeometric distributions and theoretical explanations of
 * the methods for calculating and sampling from these.
 *
 * Copyright 2004-2013 by Agner Fog.
 * GNU General Public License http://www.gnu.org/licenses/gpl.html
 *******************************************************************************/

#ifndef STOCC_H
#define STOCC_H

#include <math.h>
#include "randomc.h"

#ifdef R_BUILD
#include "stocR.h"           // Include this when building R-language interface
#endif


/***********************************************************************
 Choose which uniform random number generator to base these classes on
 ***********************************************************************/

// STOC_BASE defines which base class to use for the non-uniform
// random number generator classes StochasticLib1, 2, and 3.
#ifndef STOC_BASE
#ifdef R_BUILD
// Inherit from StocRBase when building for R-language interface
#define STOC_BASE StocRBase
#else
#define STOC_BASE CRandomMersenne     // C++ Mersenne Twister
// Or choose any other random number generator base class, for example:
//#include "randoma.h"
//#define STOC_BASE CRandomSFMTA      // Binary library SFMT generator
#endif
#endif

/***********************************************************************
 Other simple functions
 ***********************************************************************/

uint32_t NumSD (double accuracy);           // used internally for determining summation interval

/***********************************************************************
 Constants and tables
 ***********************************************************************/

// Maximum number of colors in the multivariate distributions
#ifndef MAXCOLORS
#define MAXCOLORS 32                // You may change this value
#endif

// constant for LnFac function:
static const int FAK_LEN = 1024;       // length of factorial table

// The following tables are tables of residues of a certain expansion
// of the error function. These tables are used in the Laplace method
// for calculating Wallenius' noncentral hypergeometric distribution.
// There are ERFRES_N tables covering desired precisions from
// 2^(-ERFRES_B) to 2^(-ERFRES_E). Only the table that matches the
// desired precision is used. The tables are defined in erfres.h which
// is included in wnchyppr.cpp.

// constants for ErfRes tables:
static const int ERFRES_B = 16;        // begin: -log2 of lowest precision
static const int ERFRES_E = 40;        // end:   -log2 of highest precision
static const int ERFRES_S =  2;        // step size from begin to end
static const int ERFRES_N = (ERFRES_E-ERFRES_B)/ERFRES_S+1; // number of tables
static const int ERFRES_L = 48;        // length of each table

// tables of error function residues:
extern "C" double ErfRes [ERFRES_N][ERFRES_L];

// number of std. deviations to include in integral to obtain desired precision:
extern "C" double NumSDev[ERFRES_N];

/***********************************************************************
 Class CFishersNCHypergeometric
 ***********************************************************************/

class CFishersNCHypergeometric {
	// This class contains methods for calculating the univariate Fisher's
	// noncentral hypergeometric probability function
public:
	CFishersNCHypergeometric(uint32_t n, uint32_t m, uint32_t N, double odds/*, double accuracy = 1E-8*/); // constructor
	double MakeTable(double * table, uint32_t MaxLength, uint32_t * xfirst, uint32_t * xlast, double cutoff = 0.); // make table of probabilities
	double mean(void);                                      // calculate approximate mean
	double variance(void);                                  // approximate variance
	uint32_t mode(void);                                     // calculate mode (exact)
	
	// WCH: for profiling to separate table length and actual table making
	// OPTIMIZE: stack setup and teardown cost is 33.7%
	inline uint32_t getTableLength() {
		// special cases
		if (xmin == xmax) return 1;
		if (odds <= 0.) {
			//if (n > N-m) FatalError("Not enough items with nonzero weight in  CWalleniusNCHypergeometric::MakeTable");
			return 1;
		}
		
		// Return UseTable and LengthNeeded
		uint32_t DesiredLength = xmax + 1 - xmin;     // max length of table
		if (DesiredLength > 200) {
			double sd = sqrt(variance()); // calculate approximate standard deviation
			// estimate number of standard deviations to include from normal distribution
			//uint32_t numSD = NumSD(accuracy);
			uint32_t i = (uint32_t)(CFishersNCHypergeometric::numSD * sd + 0.5);
			if (DesiredLength > i) DesiredLength = i;
		}
		return DesiredLength;
	}
	
protected:
	
	// parameters
	double odds;                                            // odds ratio
	//double accuracy;                                        // accuracy
	uint32_t n, m, N;                                        // Parameters
	uint32_t xmin, xmax;                                     // minimum and maximum of x
	
	// parameters used by subfunctions
	uint32_t xLast;
	
private:
	double cachedMean;
	void calculate_mode(void);                                     // calculate mode (exact)
	// cached
	uint32_t cachedMode;
	//uint32_t cached_N_minus_sum_nm;
	uint32_t cached_sum_nm;
	uint32_t cacheFlag;

	
public:
	static void setAccuracy(double a);
private:
	static double accuracy;
	static uint32_t numSD;
};

#endif
