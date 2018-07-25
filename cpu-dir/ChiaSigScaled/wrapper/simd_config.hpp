#ifdef __SIMD_CONFIG_HPP__
#define __SIMD_CONFIG_HPP__

#if (defined(_M_AMD64) || defined(_M_X64) || defined(__amd64)) && ! defined(__x86_64__)
    #define __x86_64__ 1
#endif

// Find sse instruction set from compiler macros if SSE_INSTR_SET not defined
// Note: Not all compilers define these macros automatically
#ifndef SSE_INSTR_SET
    #if defined ( __AVX2__ )
        #define SSE_INSTR_SET 8
    #elif defined ( __AVX__ )
        #define SSE_INSTR_SET 7
    #elif defined ( __SSE4_2__ )
        #define SSE_INSTR_SET 6
    #elif defined ( __SSE4_1__ )
        #define SSE_INSTR_SET 5
    #elif defined ( __SSSE3__ )
        #define SSE_INSTR_SET 4
    #elif defined ( __SSE3__ )
        #define SSE_INSTR_SET 3
    #elif defined ( __SSE2__ ) || defined ( __x86_64__ )
        #define SSE_INSTR_SET 2
    #elif defined ( __SSE__ )
        #define SSE_INSTR_SET 1
    #elif defined ( _M_IX86_FP )           // Defined in MS compiler on 32bits system. 1: SSE, 2: SSE2
        #define SSE_INSTR_SET _M_IX86_FP
    #else
        #define SSE_INSTR_SET 0
    #endif // instruction set defines
#endif // SSE_INSTR_SET

// Include the appropriate header file for intrinsic functions
#if SSE_INSTR_SET > 7                  // AVX2 and later
    #ifdef __GNUC__
        #include <x86intrin.h>         // x86intrin.h includes header files for whatever instruction
                                       // sets are specified on the compiler command line, such as:
                                       // xopintrin.h, fma4intrin.h
    #else
        #include <immintrin.h>         // MS version of immintrin.h covers AVX, AVX2 and FMA3
    #endif // __GNUC__
#elif SSE_INSTR_SET == 7
    #include <immintrin.h>             // AVX
#elif SSE_INSTR_SET == 6
    #include <nmmintrin.h>             // SSE4.2
#elif SSE_INSTR_SET == 5
    #include <smmintrin.h>             // SSE4.1
#elif SSE_INSTR_SET == 4
    #include <tmmintrin.h>             // SSSE3
#elif SSE_INSTR_SET == 3
    #include <pmmintrin.h>             // SSE3
#elif SSE_INSTR_SET == 2
    #include <emmintrin.h>             // SSE2
#elif SSE_INSTR_SET == 1
    #include <xmmintrin.h>             // SSE
#endif

#endif //__SIMD_CONFIG_HPP__
