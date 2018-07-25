// absorbed.

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
