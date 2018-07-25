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
 * This file contains source code for the class CFishersNCHypergeometric
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


/***********************************************************************
 Methods for class CFishersNCHypergeometric
 ***********************************************************************/

CFishersNCHypergeometric::CFishersNCHypergeometric(int32_t n, int32_t m, int32_t N, double odds, double accuracy) {
	// constructor
	// set parameters
	this->n = n;  this->m = m;  this->N = N;
	this->odds = odds;  this->accuracy = accuracy;
	
	// check validity of parameters
	/*if (n < 0 || m < 0 || N < 0 || odds < 0. || n > N || m > N) {
	 FatalError("Parameter out of range in class CFishersNCHypergeometric");
	 }*/
	if (n < 0) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class CFishersNCHypergeometric, n < 0, %d < 0\n", __func__, __LINE__, n);
		FatalError("Parameter out of range in class CFishersNCHypergeometric, n < 0");
	}
	if (m < 0) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class CFishersNCHypergeometric, m < 0, %d < 0\n", __func__, __LINE__, m);
		FatalError("Parameter out of range in class CFishersNCHypergeometric, m < 0");
	}
	if (N < 0) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class CFishersNCHypergeometric, N < 0, %d < 0\n", __func__, __LINE__, N);
		FatalError("Parameter out of range in class CFishersNCHypergeometric, N<0");
	}
	if (odds < 0.) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class CFishersNCHypergeometric, odds < 0., %f < 0\n", __func__, __LINE__, odds);
		FatalError("Parameter out of range in class CFishersNCHypergeometric, odds < 0.");
	}
	if (n > N) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class CFishersNCHypergeometric, n > N, %d > %d\n", __func__, __LINE__, n, N);
		FatalError("Parameter out of range in class CFishersNCHypergeometric, n > N");
	}
	if (m > N) {
		fprintf(stderr, "[M::%s:%d] Parameter out of range in class CFishersNCHypergeometric, m > N, %d > %d\n", __func__, __LINE__, m, N);
		FatalError("Parameter out of range in class CFishersNCHypergeometric, m > N");
	}
#if 0
	if (accuracy < 0) accuracy = 0;
	if (accuracy > 1) accuracy = 1;
#else
	// WCH: is this the intent?
	if (this->accuracy < 0) this->accuracy = 0;
	if (this->accuracy > 1) this->accuracy = 1;
#endif
	// initialize
	logodds = log(odds);  scale = rsum = 0.;
	ParametersChanged = 1;
	// calculate xmin and xmax
	xmin = m + n - N;  if (xmin < 0) xmin = 0;
	xmax = n;  if (xmax > m) xmax = m;
}


int32_t CFishersNCHypergeometric::mode(void) {
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
		int32_t L = m + n - N;
		int32_t m1 = m+1, n1 = n+1;
		
		// calculate analogously to Cornfield mean
		A = 1. - odds;
		B = (m1+n1)*odds - L;
		C = -(double)m1*n1*odds;
		D = B*B -4*A*C;
		D = D > 0. ? sqrt(D) : 0.;
		x = (D - B)/(A+A);
	}
	return (int32_t)x;
}


double CFishersNCHypergeometric::mean(void) {
	// Find approximate mean
	// Calculation analogous with mode
	if (odds == 1.) {                   // simple hypergeometric
		return double(m)*n/N;
	}
	double a, b;                        // temporaries in calculation
	double mean;                        // mean
	// calculate Cornfield mean
	a = (m+n)*odds + (N-m-n);
	b = a*a - 4.*odds*(odds-1.)*m*n;
	b = b > 0. ? sqrt(b) : 0.;
	mean = (a-b)/(2.*(odds-1.));
	return mean;
}


double CFishersNCHypergeometric::variance(void) {
	// find approximate variance (poor approximation)
	double my = mean(); // approximate mean
	// find approximate variance from Fisher's noncentral hypergeometric approximation
	double r1 = my * (m-my); double r2 = (n-my)*(my+N-n-m);
	if (r1 <= 0. || r2 <= 0.) return 0.;
	double var = N*r1*r2/((N-1)*(m*r2+(N-m)*r1));
	if (var < 0.) var = 0.;
	return var;
}

double CFishersNCHypergeometric::MakeTable(double * table, int32_t MaxLength, int32_t * xfirst, int32_t * xlast, double cutoff) {
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
	int32_t mode = this->mode();        // mode
	int32_t L = n + m - N;              // parameter
	int32_t DesiredLength;              // desired length of table
	
	// limits for x
	
	x1 = (L > 0) ? L : 0;               // xmin
	x2 = (n < m) ? n : m;               // xmax
	
	// special cases
	if (x1 == x2) goto DETERMINISTIC;
	if (odds <= 0.) {
		if (n > N-m) FatalError("Not enough items with nonzero weight in  CWalleniusNCHypergeometric::MakeTable");
		x1 = 0;
	DETERMINISTIC:
		if (MaxLength == 0) {
			if (xfirst) *xfirst = 1;
			return 1;
		}
		*xfirst = *xlast = x1;
		*table = 1.;
		return 1;
	}
	
	if (MaxLength <= 0) {
		// Return UseTable and LengthNeeded
		DesiredLength = x2 - x1 + 1;     // max length of table
		if (DesiredLength > 200) {
			double sd = sqrt(variance()); // calculate approximate standard deviation
			// estimate number of standard deviations to include from normal distribution
			i = (int32_t)(NumSD(accuracy) * sd + 0.5);
			if (DesiredLength > i) DesiredLength = i;
		}
		if (xfirst) *xfirst = 1;         // for analogy with CWalleniusNCHypergeometric::MakeTable
		return DesiredLength;
	}
	
	// place mode in the table
	if (mode - x1 <= MaxLength/2) {
		// There is enough space for left tail
		i0 = mode - x1;
	}
	else if (x2 - mode <= MaxLength/2) {
		// There is enough space for right tail
		i0 = MaxLength - x2 + mode - 1;
		if (i0 < 0) i0 = 0;
	}
	else {
		// There is not enough space for any of the tails. Place mode in middle of table
		i0 = MaxLength/2;
	}
	// Table start index
	i1 = i0 - mode + x1;  if (i1 < 0) i1 = 0;
	
	// Table end index
	i2 = i0 + x2 - mode;  if (i2 > MaxLength-1) i2 = MaxLength-1;
	
	// make center
	table[i0] = sum = f = 1.;
	
	// make left tail
	x = mode;
	a1 = m + 1 - x;  a2 = n + 1 - x;
	b1 = x;  b2 = x - L;
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
	x = mode + 1;
	a1 = m + 1 - x;  a2 = n + 1 - x;
	b1 = x;  b2 = x - L;
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
	*xfirst = mode - (i0 - i1);
	*xlast  = mode + (i2 - i0);
	
	return sum;
}

