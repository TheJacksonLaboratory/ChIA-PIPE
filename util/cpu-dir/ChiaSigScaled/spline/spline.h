/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */


#ifndef _tk_spline_h
#define _tk_spline_h

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include "../CCCsig/bar.h"

// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files

// band matrix solver
class band_matrix {
private:
	std::vector< std::vector<double> > m_upper;  // upper band
	std::vector< std::vector<double> > m_lower;  // lower band
public:
	band_matrix() {};                             // constructor
	band_matrix(long dim, long n_u, long n_l);       // constructor
	~band_matrix() {};                            // destructor
	void resize(long dim, long n_u, long n_l);      // init with dim,n_u,n_l
	long dim() const;                             // matrix dimension
	long num_upper() const {
		return m_upper.size()-1;
	}
	long num_lower() const {
		return m_lower.size()-1;
	}
	// access operator
	double & operator () (long i, long j);            // write
	double   operator () (long i, long j) const;      // read
	// we can store an additional diogonal (in m_lower)
	double& saved_diag(long i);
	double  saved_diag(long i) const;
	void lu_decompose();
	std::vector<double> r_solve(const std::vector<double>& b) const;
	std::vector<double> l_solve(const std::vector<double>& b) const;
	std::vector<double> lu_solve(const std::vector<double>& b,
								 bool is_lu_decomposed=false);
	
};


//typedef uint32_t DeltaIndicatorCounter;

// spline interpolation
//typedef float precomputeType;
typedef double precomputeType;

class spline {
public:
	//WCH: for debugging, we use public for now
	//private:
	std::vector<double> m_x,m_y;           // x,y coordinates of points
	// interpolation parameters
	// f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
	std::vector<double> m_a,m_b,m_c;
	
	bool precomputed;
	uint32_t numMinus1;
	precomputeType* m_precompute;
	
public:
	spline();
	~spline();
	void set_points(const std::vector<double>& x,
					const std::vector<double>& y, bool cubic_spline=true);
	precomputeType operator() (uint32_t x) const;
#if 0
	inline precomputeType cachedAt(uint32_t x) const { return m_precompute[x]; }
#else
	inline double cachedAt(uint32_t x) const { return m_precompute[x]; }
#endif
	void precompute(uint32_t num, BIT_ARRAY* pDeltas);
	void setThreads(int n) { nThreads = n; }
	precomputeType at(uint32_t x) const;
	
private:
	int nThreads;
};

#endif /* _tk_spline_h */
