#include "spline.h"
#include <iostream>
#include <math.h>

#if 0
#define USE_SIMD_OPERATIONS 1
#ifdef USE_SIMD_OPERATIONS
#include "../CCCsig/vectorclass.h"
#endif
#include <string.h>
#endif

// ---------------------------------------------------------------------
// implementation part, which should be separated into a cpp file
// ---------------------------------------------------------------------


// band_matrix implementation
// -------------------------

band_matrix::band_matrix(long dim, long n_u, long n_l) {
	resize(dim, n_u, n_l);
}
void band_matrix::resize(long dim, long n_u, long n_l) {
	assert(dim>0);
	assert(n_u>=0);
	assert(n_l>=0);
	m_upper.resize(n_u+1);
	m_lower.resize(n_l+1);
	for(size_t i=0; i<m_upper.size(); i++) {
		m_upper[i].resize(dim);
	}
	for(size_t i=0; i<m_lower.size(); i++) {
		m_lower[i].resize(dim);
	}
}
long band_matrix::dim() const {
	if(m_upper.size()>0) {
		return m_upper[0].size();
	} else {
		return 0;
	}
}


// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & band_matrix::operator () (long i, long j) {
	long k=j-i;       // what band is the entry
	assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
	assert( (-num_lower()<=k) && (k<=num_upper()) );
	// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
	if(k>=0)   return m_upper[k][i];
	else	    return m_lower[-k][i];
}
double band_matrix::operator () (long i, long j) const {
	long k=j-i;       // what band is the entry
	assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
	assert( (-num_lower()<=k) && (k<=num_upper()) );
	// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
	if(k>=0)   return m_upper[k][i];
	else	    return m_lower[-k][i];
}
// second diag (used in LU decomposition), saved in m_lower
double band_matrix::saved_diag(long i) const {
	assert( (i>=0) && (i<dim()) );
	return m_lower[0][i];
}
double & band_matrix::saved_diag(long i) {
	assert( (i>=0) && (i<dim()) );
	return m_lower[0][i];
}

// LR-Decomposition of a band matrix
void band_matrix::lu_decompose() {
	long i_max,j_max;
	long j_min;
	double x;
	
	// preconditioning
	// normalize column i so that a_ii=1
	for(long i=0; i<this->dim(); i++) {
		assert(this->operator()(i,i)!=0.0);
		this->saved_diag(i)=1.0/this->operator()(i,i);
		j_min=std::max(0L,i-this->num_lower());
		j_max=std::min(this->dim()-1,i+this->num_upper());
		for(long j=j_min; j<=j_max; j++) {
			this->operator()(i,j) *= this->saved_diag(i);
		}
		this->operator()(i,i)=1.0;          // prevents rounding errors
	}
	
	// Gauss LR-Decomposition
	for(int k=0; k<this->dim(); k++) {
		i_max=std::min(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
		for(int i=k+1; i<=i_max; i++) {
			assert(this->operator()(k,k)!=0.0);
			x=-this->operator()(i,k)/this->operator()(k,k);
			this->operator()(i,k)=-x;                         // assembly part of L
			j_max=std::min(this->dim()-1,k+this->num_upper());
			for(int j=k+1; j<=j_max; j++) {
				// assembly part of R
				this->operator()(i,j)=this->operator()(i,j)+x*this->operator()(k,j);
			}
		}
	}
}
// solves Ly=b
std::vector<double> band_matrix::l_solve(const std::vector<double>& b) const {
	assert( this->dim()==(int)b.size() );
	std::vector<double> x(this->dim());
	long j_start;
	double sum;
	for(long i=0; i<this->dim(); i++) {
		sum=0;
		j_start=std::max(0L,i-this->num_lower());
		for(long j=j_start; j<i; j++) sum += this->operator()(i,j)*x[j];
		x[i]=(b[i]*this->saved_diag(i)) - sum;
	}
	return x;
}
// solves Rx=y
std::vector<double> band_matrix::r_solve(const std::vector<double>& b) const {
	assert( this->dim()==(int)b.size() );
	std::vector<double> x(this->dim());
	long j_stop;
	double sum;
	for(long i=this->dim()-1; i>=0; i--) {
		sum=0;
		j_stop=std::min(this->dim()-1,i+this->num_upper());
		for(long j=i+1; j<=j_stop; j++) sum += this->operator()(i,j)*x[j];
		x[i]=( b[i] - sum ) / this->operator()(i,i);
	}
	return x;
}

std::vector<double> band_matrix::lu_solve(const std::vector<double>& b,
										  bool is_lu_decomposed) {
	assert( this->dim()==(int)b.size() );
	std::vector<double>  x,y;
	if(is_lu_decomposed==false) {
		this->lu_decompose();
	}
	y=this->l_solve(b);
	x=this->r_solve(y);
	return x;
}





// spline implementation
// -----------------------

spline::spline()
: precomputed(false),numMinus1(0),m_precompute(NULL),nThreads(1)
{
	
}

spline::~spline()
{
	if (m_precompute) {
		free(m_precompute);
		m_precompute = NULL;
	}
}

void spline::set_points(const std::vector<double>& x,
						const std::vector<double>& y, bool cubic_spline) {
	assert(x.size()==y.size());
	m_x=x;
	m_y=y;
	size_t n=x.size();
	// TODO sort x and y, rather than returning an error
	for(int i=0; i<n-1; i++) {
		assert(m_x[i]<m_x[i+1]);
	}
	
	if(cubic_spline==true) { // cubic spline interpolation
		// setting up the matrix and right hand side of the equation system
		// for the parameters b[]
		band_matrix A(n,1,1);
		std::vector<double>  rhs(n);
		for(int i=1; i<n-1; i++) {
			A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
			A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
			A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
			rhs[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		}
		// boundary conditions, zero curvature b[0]=b[n-1]=0
		A(0,0)=2.0;
		A(0,1)=0.0;
		rhs[0]=0.0;
		A(n-1,n-1)=2.0;
		A(n-1,n-2)=0.0;
		rhs[n-1]=0.0;
		
		// solve the equation system to obtain the parameters b[]
		m_b=A.lu_solve(rhs);
		
		// calculate parameters a[] and c[] based on b[]
		m_a.resize(n);
		m_c.resize(n);
		for(int i=0; i<n-1; i++) {
			m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
			m_c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
			- 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
		}
	} else { // linear interpolation
		m_a.resize(n);
		m_b.resize(n);
		m_c.resize(n);
		for(int i=0; i<n-1; i++) {
			m_a[i]=0.0;
			m_b[i]=0.0;
			m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
		}
	}
	
	// for the right boundary we define
	// f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
	double h=x[n-1]-x[n-2];
	// m_b[n-1] is determined by the boundary condition
	m_a[n-1]=0.0;
	m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
	
	assert(x.size()<=0xFFFFFF);
	numMinus1 = (uint32_t)x.size();
	numMinus1--;
}

// TODO: OPTIMIZE multithreading
void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);

struct spline_precompute_left_fast_worker_t {
	spline_precompute_left_fast_worker_t(spline* s, precomputeType* vs, BIT_ARRAY* pds)
	: pSpline(s), values(vs), pDeltas(pds)
	{}
	
	spline* pSpline;
	precomputeType* values;
	BIT_ARRAY* pDeltas;
};

static void spline_precompute_left_fast_worker(void *data, int idx, int tid)
{
	spline_precompute_left_fast_worker_t *w = (spline_precompute_left_fast_worker_t*)data;
	if (barget(w->pDeltas, idx)) {
		double h=idx-w->pSpline->m_x[0];
		w->values[idx] = ((w->pSpline->m_b[0])*h + w->pSpline->m_c[0])*h + w->pSpline->m_y[0];
	}
}

struct spline_precompute_flanked_fast_worker_t {
	spline_precompute_flanked_fast_worker_t(spline *s, precomputeType *vs, std::vector<uint32_t>& cs, std::vector<uint32_t>& cn, BIT_ARRAY* pds, uint32_t ft)
	: pSpline(s), values(vs), chunkStarts(cs), chunkNums(cn), pDeltas(pds), flipTotal(ft)
	{}
	
	spline* pSpline;
	precomputeType* values;
	std::vector<uint32_t>& chunkStarts;
	std::vector<uint32_t>& chunkNums;
	BIT_ARRAY* pDeltas;
	uint32_t flipTotal;
};

static void spline_precompute_flanked_fast_worker(void *data, int idx, int tid)
{
	spline_precompute_flanked_fast_worker_t *w = (spline_precompute_flanked_fast_worker_t*)data;
	
	idx = w->flipTotal - idx;
	uint32_t i=0;
	uint32_t chunkStart=w->chunkStarts[idx];
	uint32_t num=w->chunkNums[idx];
	double h=0.;
	for(; i<num; ++i,h+=1.) {
		if (barget(w->pDeltas, chunkStart+i))
			w->values[chunkStart+i] = ((w->pSpline->m_a[idx]*h + w->pSpline->m_b[idx])*h + w->pSpline->m_c[idx])*h + w->pSpline->m_y[idx];
	}
}

struct spline_precompute_right_fast_worker_t {
	spline_precompute_right_fast_worker_t(spline *s, precomputeType *vs, uint32_t cs, uint32_t ci, BIT_ARRAY* pds)
	: pSpline(s), values(vs), chunkStart(cs), chunkIdx(ci), pDeltas(pds)
	{}
	
	spline* pSpline;
	precomputeType* values;
	uint32_t chunkStart;
	uint32_t chunkIdx;
	BIT_ARRAY* pDeltas;
};

static void spline_precompute_right_fast_worker(void *data, int idx, int tid)
{
	spline_precompute_right_fast_worker_t *w = (spline_precompute_right_fast_worker_t*)data;
	if (barget(w->pDeltas, w->chunkStart+idx)) {
		double h=idx;
		w->values[w->chunkStart+idx] = ((w->pSpline->m_b[w->chunkIdx])*h + w->pSpline->m_c[w->chunkIdx])*h + w->pSpline->m_y[w->chunkIdx];
	}
}

void spline::precompute(uint32_t n, BIT_ARRAY* pDeltas)
{
	// pre-compute for all possible values
	if (m_precompute) {
		free(m_precompute);
		m_precompute = NULL;
	}
	m_precompute = (precomputeType*)calloc(n, sizeof(precomputeType));
	
	// TODO: we can start from minDelta!!!
	// TODO: OPTIMIZE: we can interatively build up the precompute value
	//       NO need to perform a lower-bound check for every element!
	
	if (nThreads>1) {
		//precomputeLeft();
		spline_precompute_left_fast_worker_t wl(this, m_precompute, pDeltas);
		kt_for(nThreads, spline_precompute_left_fast_worker, &wl, (int)m_x[0]);
		
		//precomputeFlanked();
		std::vector<uint32_t> chunkStarts; chunkStarts.resize(numMinus1);
		std::vector<uint32_t> chunkNums; chunkNums.resize(numMinus1);
		//for(uint32_t idx=0; idx<=numMinus1; ++idx) {
		for(uint32_t idx=0; idx<numMinus1; ++idx) {
			chunkStarts[idx] = m_x[idx];
			chunkNums[idx] = m_x[idx+1]-m_x[idx];
		}
		spline_precompute_flanked_fast_worker_t wf(this, m_precompute, chunkStarts, chunkNums, pDeltas, numMinus1-1);
		kt_for(nThreads, spline_precompute_flanked_fast_worker, &wf, (int)numMinus1);
		
		//precomputeRight();
		uint32_t chunkStart = m_x[numMinus1];
		uint32_t num = n - m_x[numMinus1];
		spline_precompute_right_fast_worker_t wr(this, m_precompute, chunkStart, numMinus1, pDeltas);
		kt_for(nThreads, spline_precompute_right_fast_worker, &wr, (int)num);
		
	} else {
		//precomputeLeft();
		uint32_t i=0;
		uint32_t chunkStart=0;
		uint32_t num=m_x[0];
		double h;
		for(h=-m_x[0]; i<num; ++i,h+=1) {
			m_precompute[chunkStart+i] = ((m_b[0])*h + m_c[0])*h + m_y[0];
		}
		//precomputeFlanked();
		for(uint32_t idx=0; idx<numMinus1; ++idx) {
			i = 0;
			chunkStart = m_x[idx];
			num=m_x[idx+1]-m_x[idx];
			
			for(h=0; i<num; ++i,h+=1) {
				m_precompute[chunkStart+i] = ((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
			}
		}
		//precomputeRight();
		i = 0;
		chunkStart = m_x[numMinus1];
		num = n - m_x[numMinus1];
		for(h=0; i<num; ++i,h+=1) {
			m_precompute[chunkStart+i] = ((m_b[numMinus1])*h + m_c[numMinus1])*h + m_y[numMinus1];
		}
	}
#if 0
	// checking
	for(uint32_t i=0; i<n; ++i) {
		// pre-computation trigger:
		// non-zero counter w/o indicator, zero counter with indicator, both indicator and counter
		if (pDeltas[i]>0) {
			if (m_precompute[i]!=at(i)) {
				fprintf(stderr, "precompute(i:%u) new:%g orig:%g\n", i, m_precompute[i], at(i));
			}
		}
	}
#endif
	
	precomputed = true;
}

// TODO: OPTIMIZE: we will never have cases where x needs extrapolation!!!
precomputeType spline::at(uint32_t x) const
{
	//size_t n=m_x.size();
	// find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
	std::vector<double>::const_iterator it;
	it=std::lower_bound(m_x.begin(),m_x.end(),(double)x);
	//int idx=std::max( int(it-m_x.begin())-1, 0);
	long idx = it-m_x.begin(); if (0!=idx) idx--;
	
	double h=x-m_x[idx];
	double interpol;
	if(x<m_x[0]) {
		// extrapolation to the left
		interpol=((m_b[0])*h + m_c[0])*h + m_y[0];
	} else if(x>m_x[numMinus1]) {
		// extrapolation to the right
		interpol=((m_b[numMinus1])*h + m_c[numMinus1])*h + m_y[numMinus1];
	} else {
		// interpolation
		interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
	}
	return interpol;
}

precomputeType spline::operator() (uint32_t x) const {
	assert(precomputed);
	return m_precompute[x];
}

