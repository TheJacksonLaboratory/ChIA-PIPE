/*************************** fnchyppr.cpp **********************************
 * Author:        Agner Fog
 * Date created:  2002-10-20
 * Last modified: 2014-06-14
 * Project:       stocc.zip
 * Source URL:    www.agner.org/random
 *
 * Description:
 * Calculation of univariate and multivariate Fisher's noncentral hypergeometric
 * probability distribution.
 *
 * This file contains source code for the class oCFishersNCHypergeometric
 * and CMultiFishersNCHypergeometric defined in stocc.h.
 *
 * Documentation:
 * ==============
 * The file stocc.h contains class definitions.
 * The file ran-instructions.pdf contains further documentation and
 * instructions.
 *
 * Copyright 2002-2014 by Agner Fog.
 * GNU General Public License http://www.gnu.org/licenses/gpl.html
 *****************************************************************************/


#include <cstdio> //WCH: fprintf, stderr
#include <string.h>                    // memcpy function
#include "stocc.h"                     // class definition

#include <algorithm> // std::min
#include <assert.h>

// WCH: intentionally extracted only the functions that we really used
//      make it easier to focus on what are the ones needed for optimization

// extracted from : wnchyppr.cpp
static const double NumSD_fract[] = {
	2.699796e-03, 4.652582e-04, 6.334248e-05, 6.795346e-06, 5.733031e-07,
	3.797912e-08, 1.973175e-09, 8.032001e-11, 2.559625e-12, 6.381783e-14
};
static const int NumSD_numFracts = sizeof(NumSD_fract)/sizeof(*NumSD_fract);
int NumSD (double accuracy) {
	// Gives the length of the integration interval necessary to achieve
	// the desired accuracy when integrating/summating a probability
	// function, relative to the standard deviation
	// Returns an integer approximation to 2*NormalDistrFractile(accuracy/2)
	int i;
	for (i = 0; i < NumSD_numFracts; i++) {
		if (accuracy >= NumSD_fract[i]) break;
	}
	return i + 6;
}

/***********************************************************************
 Methods for class oCFishersNCHypergeometric
 ***********************************************************************/

oCFishersNCHypergeometric::oCFishersNCHypergeometric(int32_t n, int32_t m, int32_t N, double odds, double accuracy)
: n(n), m(m), N(N), odds(odds), accuracy(accuracy)
,xmax(std::min(n,m))
,cached_sum_nm_minus_N(n + m - N)
,xmin(std::max(n + m - N,0))
{
#if 0
	// TODO: OPTIMIZE: don't check?
	// check validity of parameters
	/*if (n < 0 || m < 0 || N < 0 || odds < 0. || n > N || m > N) {
	 FatalError("Parameter out of range in class oCFishersNCHypergeometric");
	 }*/
	if (n < 0) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class oCFishersNCHypergeometric, n < 0, %d < 0\n", __func__, __LINE__, n);
		FatalError("Parameter out of range in class oCFishersNCHypergeometric, n < 0");
	}
	if (m < 0) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class oCFishersNCHypergeometric, m < 0, %d < 0\n", __func__, __LINE__, m);
		FatalError("Parameter out of range in class oCFishersNCHypergeometric, m < 0");
	}
	if (N < 0) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class oCFishersNCHypergeometric, N < 0, %d < 0\n", __func__, __LINE__, N);
		FatalError("Parameter out of range in class oCFishersNCHypergeometric, N<0");
	}
	if (odds < 0.) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class oCFishersNCHypergeometric, odds < 0., %f < 0\n", __func__, __LINE__, odds);
		FatalError("Parameter out of range in class oCFishersNCHypergeometric, odds < 0.");
	}
	if (n > N) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class oCFishersNCHypergeometric, n > N, %d > %d\n", __func__, __LINE__, n, N);
		FatalError("Parameter out of range in class oCFishersNCHypergeometric, n > N");
	}
	if (m > N) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class oCFishersNCHypergeometric, m > N, %d > %d\n", __func__, __LINE__, m, N);
		FatalError("Parameter out of range in class oCFishersNCHypergeometric, m > N");
	}
#endif
#if 0
	if (accuracy < 0) accuracy = 0;
	if (accuracy > 1) accuracy = 1;
#else
#if 0
	// TODO: OPTIMIZE: don't check?
	// WCH: is this the intent?
	if (this->accuracy < 0) this->accuracy = 0;
	if (this->accuracy > 1) this->accuracy = 1;
#endif
#endif
	cachedMode = calculate_mode();

}


int32_t oCFishersNCHypergeometric::mode(void) {
	return cachedMode;
}

int32_t oCFishersNCHypergeometric::calculate_mode(void) {
	// Find mode (exact)
	// Uses the method of Liao and Rosen, The American Statistician, vol 55,
	// no 4, 2001, p. 366-369.
	// Note that there is an error in Liao and Rosen's formula.
	// Replace sgn(b) with -1 in Liao and Rosen's formula.
	
	double x;                           // mode
	
	if (odds == 1.) {
		// simple hypergeometric
		x = (m + 1.) * (n + 1.) / (N + 2.);
	}
	else {
		double A, B, C, D;                  // coefficients for quadratic equation
		//int32_t L = n + m - N;
		int32_t m1 = m+1, n1 = n+1;
		
		// calculate analogously to Cornfield mean
		A = 1. - odds;
		//B = (m1+n1)*odds - L;
		B = (m1+n1)*odds - cached_sum_nm_minus_N;
		C = -(double)m1*n1*odds;
		D = B*B -4*A*C;
		D = D > 0. ? sqrt(D) : 0.;
		x = (D - B)/(A+A);
	}
	return (int32_t)x;
}
	
double oCFishersNCHypergeometric::mean(void) {
	// Find approximate mean
	// Calculation analogous with mode
	if (odds == 1.) {                   // simple hypergeometric
		return double(m)*n/N;
	}
	double a, b;                        // temporaries in calculation
	double mean;                        // mean
	// calculate Cornfield mean
	//a = (n + m)*odds + (N-m-n);
	//a = (n + m)*odds - (n + m - N);
	a = (n + m)*odds - cached_sum_nm_minus_N;
	b = a*a - 4.*odds*(odds-1.)*m*n;
	b = b > 0. ? sqrt(b) : 0.;
	mean = (a-b)/(2.*(odds-1.));
	return mean;
}


double oCFishersNCHypergeometric::variance(void) {
	// find approximate variance (poor approximation)
	double my = mean(); // approximate mean
	// find approximate variance from Fisher's noncentral hypergeometric approximation
	double r1 = my * (m-my); double r2 = (n-my)*(my+N-n-m);
	if (r1 <= 0. || r2 <= 0.) return 0.;
	double var = N*r1*r2/((N-1)*(m*r2+(N-m)*r1));
	if (var < 0.) var = 0.;
	return var;
}

double oCFishersNCHypergeometric::MakeTable(double * table, int32_t MaxLength, int32_t * xfirst, int32_t * xlast, double cutoff) {
	// Makes a table of Fisher's noncentral hypergeometric probabilities.
	// Results are returned in the array table of size MaxLength.
	// The values are scaled so that the highest value is 1. The return value
	// is the sum, s, of all the values in the table. The normalized
	// probabilities are obtained by multiplying all values in the table by
	// 1/s.
	// The tails are cut off where the values are < cutoff, so that
	// *xfirst may be > xmin and *xlast may be < xmax.
	// The value of cutoff will be 0.01 * accuracy if not specified.
	// The first and last x value represented in the table are returned in
	// *xfirst and *xlast. The resulting probability values are returned in the
	// first (*xlast - *xfirst + 1) positions of table. If this would require
	// more than MaxLength values then the table is filled with as many
	// correct values as possible.
	//
	// The function will return the desired length of table when MaxLength = 0.
	
	double f;                           // probability function value
	double sum;                         // sum of table values
	double a1, a2, b1, b2;              // factors in recursive calculation of f(x)
	int32_t x;                          // x value
	int32_t x1, x2;                     // lowest and highest x
	int32_t i, i0, i1, i2;              // table index
	//int32_t mode = this->mode();        // mode
	//int32_t L = n + m - N;              // parameter
	//int32_t DesiredLength;              // desired length of table
	
	// limits for x
	
	//x1 = (cached_sum_nm_minus_N > 0) ? cached_sum_nm_minus_N : 0;               // xmin
	//x2 = (n < m) ? n : m;               // xmax
	x1 = xmin;
	x2 = xmax;
	
#if 1
	assert(MaxLength>0);
#endif
	
	// special cases
	if (x1 == x2) goto DETERMINISTIC;
	if (odds <= 0.) {
		if (n > N-m) FatalError("Not enough items with nonzero weight in  CWalleniusNCHypergeometric::MakeTable");
		x1 = 0;
	DETERMINISTIC:
		*xfirst = *xlast = x1;
		*table = 1.;
		return 1;
	}
	
	// place mode in the table
	if (cachedMode - x1 <= MaxLength/2) {
		// There is enough space for left tail
		i0 = cachedMode - x1;
	}
	else if (x2 - cachedMode <= MaxLength/2) {
		// There is enough space for right tail
		i0 = MaxLength - x2 + cachedMode - 1;
		if (i0 < 0) i0 = 0;
	}
	else {
		// There is not enough space for any of the tails. Place mode in middle of table
		i0 = MaxLength/2;
	}
	// Table start index
	i1 = i0 - cachedMode + x1;  if (i1 < 0) i1 = 0;
	
	// Table end index
	i2 = i0 + x2 - cachedMode;  if (i2 > MaxLength-1) i2 = MaxLength-1;
	
	// make center
	table[i0] = sum = f = 1.;
	
	// make left tail
	x = cachedMode;
	a1 = m + 1 - x;  a2 = n + 1 - x;
	b1 = x;  b2 = x - cached_sum_nm_minus_N; //b2 = x - L;
	for (i = i0 - 1; i >= i1; i--) {
		f *= b1 * b2 / (a1 * a2 * odds); // recursive formula
		a1++;  a2++;  b1--;  b2--;
		sum += table[i] = f;
		if (f < cutoff) {
			i1 = i;  break;               // cut off tail if < accuracy
		}
	}
	if (i1 > 0) {
		// move table down for cut-off left tail
		memcpy(table, table+i1, (i0-i1+1)*sizeof(*table));
		// adjust indices
		i0 -= i1;  i2 -= i1;  i1 = 0;
	}
	// make right tail
	x = cachedMode + 1;
	a1 = m + 1 - x;  a2 = n + 1 - x;
	b1 = x;  b2 = x - cached_sum_nm_minus_N; //b2 = x - L;
	f = 1.;
	for (i = i0 + 1; i <= i2; i++) {
		f *= a1 * a2 * odds / (b1 * b2); // recursive formula
		a1--;  a2--;  b1++;  b2++;
		sum += table[i] = f;
		if (f < cutoff) {
			i2 = i;  break;               // cut off tail if < accuracy
		}
	}
	// x limits
	*xfirst = cachedMode - (i0 - i1);
	*xlast  = cachedMode + (i2 - i0);
	
	return sum;
}

// WCH: for profiling purposes
int oCFishersNCHypergeometric::getTableLength()
{
	// special cases
	if (xmin == xmax) return 1;
	if (odds <= 0.) {
		//if (n > N-m) FatalError("Not enough items with nonzero weight in  CWalleniusNCHypergeometric::MakeTable");
		return 1;
	}
	
	// Return UseTable and LengthNeeded
	int32_t DesiredLength = xmax - xmin + 1;     // max length of table
	if (DesiredLength > 200) {
		double sd = sqrt(variance()); // calculate approximate standard deviation
		// estimate number of standard deviations to include from normal distribution
		int32_t i = (int32_t)(NumSD(accuracy) * sd + 0.5);
#if 0
		// WCH: DEBUG
		fprintf(stderr, "cFNCHG.getTableLength(DesiredLength:%d, i:%d, ni:%d, nj:%d) <---\n", DesiredLength, i, n, m);
#endif
		if (i>0) {
			if (DesiredLength > i) DesiredLength = i;
		} else {
			fprintf(stderr, "cFNCHG.getTableLength has sd==0 n:%u m:%u N:%u <---\n", n, m, N);
		}
	}
	return DesiredLength;
}

