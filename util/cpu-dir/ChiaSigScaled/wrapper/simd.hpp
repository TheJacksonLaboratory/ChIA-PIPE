#ifndef __SIMD_HPP_
#define __SIMD_HPP_

template <class T> struct simd_traits
{
	typedef T type;
	static const size_t size = 1;
};

#ifdef USE_SSE
template <> struct simd_traits<float>
{
	typedef vector4f type;
	static const size_t size = 4;
};

template <> struct simd_traits<double>
{
	typedef vector2d type;
	static const size_t size = 2;
};
#elif USE_AVX
template <> struct simd_traits<float>
{
	typedef vector8f type;
	static const size_t size = 8;
};

template <> struct simd_traits<double>
{
	typedef vector4d type;
	static const size_t size = 4;
};
#endif

template <class X> class simd_vector
{
public:
	
	typedef simd_traits<X>::value_type value_type;
	
	value_type operator[](size_t index) const
	{
		size_t size = simd_traits<X>::size;
		value_type v[size];
		(*this)().store_u(v);
		return v[index];
	}
};

// Common implementation for types that support vectorization
template <class T, class V> struct simd_functions_invoker
{
	inline static V
	set1(const T& a) { return V(a); }
	
	inline static V
	load_a(const T* src) { V res; res.load_a(src); return res; }
	
	inline static V
	load_u(const T* src) { V res; res.load_u(src); return res; }
	
	inline static void
	store_a(T* dst, const V& src) { src.store_a(dst); }
	
	inline static void
	store_u(T* dst, const V& src) { src.store_u(dst); }
};

// Specialization for types that don't support vectorization
template <class T> struct simd_functions_invoker<T,T>
{
	inline static T
	set1(const T& a) { return T(a); }
	
	inline static T
	load_a(const T* src) { return *src; }
	
	inline static T
	load_u(const T* src) { return *src; }
	
	inline static void
	store_a(T* dst, const T& src) { *dst = src; }
	
	inline static void
	store_u(T* dst, const T& src) { *dst = src; }
};

template <class T> inline typename simd_traits<T>::type set1(const T& a)
{
	return simd_functions_invoker<T,typename simd_traits<T>::type>::set1(a);
}

template <class T> inline typename simd_traits<T>::type load_a(const T* src)
{
	return simd_functions_invoker<T,typename simd_traits<T>::type>::load_a(src);
}

template <class T> inline typename simd_traits<T>::type load_u(const T* src)
{
	return simd_functions_invoker<T,typename simd_traits<T>::type>::load_u(src);
}

template <class T> inline void store_a(T* dst, const typename simd_traits<T>::type& src)
{
	simd_functions_invoker<T,typename simd_traits<T>::type>::store_a(dst,src);
}

template <class T> inline void store_u(T* dst, const typename simd_traits<T>::type& src)
{
	simd_functions_invoker<T,typename simd_traits<T>::type>::store_u(dst,src);
}

inline float hadd(const vector4f& rhs)
{
#if SSE_INSTR_SET >= 3 // SSE3
	__m128 tmp0 = _mm_hadd_ps(rhs,rhs);
	__m128 tmp1 = _mm_hadd_ps(tmp0,tmp0);
#else
	__m128 tmp0 = _mm_add_ps(rhs,_mm_movehl_ps(rhs,rhs));
	__m128 tmp1 = _mm_add_ss(tmp0,_mm_shuffle_ps(tmp0,tmp0,1));
#endif
	return _mm_cvtss_f32(tmp1);
}

#endif	//__SIMD_HPP_
