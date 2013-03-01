// Copyright (C) 2008-2013 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2013 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



#define arma_hot
#define arma_cold
#define arma_pure
#define arma_const
#define arma_aligned
#define arma_warn_unused
#define arma_deprecated
#define arma_malloc
#define arma_inline            inline
#define arma_noinline
#define arma_ignore(variable)  ((void)(variable))


#if defined(ARMA_BLAS_UNDERSCORE)
  #define arma_fortran2_noprefix(function) function##_
  #define arma_fortran2_prefix(function)   wrapper_##function##_
#else
  #define arma_fortran2_noprefix(function) function
  #define arma_fortran2_prefix(function)   wrapper_##function
#endif

#if defined(ARMA_USE_WRAPPER)
  #define arma_fortran(function) arma_fortran2_prefix(function)
  #define arma_atlas(function)   wrapper_##function
#else
  #define arma_fortran(function) arma_fortran2_noprefix(function)
  #define arma_atlas(function)   function
#endif

#define arma_fortran_prefix(function)   arma_fortran2_prefix(function)
#define arma_fortran_noprefix(function) arma_fortran2_noprefix(function)


#define ARMA_INCFILE_WRAP(x) <x>


#if (__cplusplus >= 201103L)
  #if !defined(ARMA_USE_CXX11)
    #define ARMA_USE_CXX11
  #endif
#endif


#if defined(ARMA_USE_CXX11)
  #if !defined(ARMA_USE_U64S64)
    #define ARMA_USE_U64S64
  #endif
#endif


#if defined(ARMA_64BIT_WORD)
  #if !defined(ARMA_USE_U64S64)
    #define ARMA_USE_U64S64
  #endif
#endif


#if defined(__INTEL_COMPILER)
  
  #if (__INTEL_COMPILER < 1000)
    #error "*** Need a newer compiler ***"
  #endif
  
  #define ARMA_GOOD_COMPILER
  #undef  ARMA_HAVE_STD_TR1
  
  #if (__INTEL_COMPILER <= 1110)
    #undef ARMA_HAVE_STD_ISFINITE
  #endif
  
#elif defined(__GNUG__)
  
  #if (__GNUC__ < 4)
    #error "*** Need a newer compiler ***"
  #endif
  
  #define ARMA_GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
  
  #define ARMA_GOOD_COMPILER
  #undef  ARMA_HAVE_STD_TR1
  
  #undef  arma_pure
  #undef  arma_const
  #undef  arma_aligned
  #undef  arma_warn_unused
  #undef  arma_deprecated
  #undef  arma_malloc
  #undef  arma_inline
  #undef  arma_noinline
  
  #define arma_pure               __attribute__((__pure__))
  #define arma_const              __attribute__((__const__))
  #define arma_aligned            __attribute__((__aligned__))
  #define arma_warn_unused        __attribute__((__warn_unused_result__))
  #define arma_deprecated         __attribute__((__deprecated__))
  #define arma_malloc             __attribute__((__malloc__))
  #define arma_inline      inline __attribute__((__always_inline__))
  #define arma_noinline           __attribute__((__noinline__))
  
  #if (ARMA_GCC_VERSION >= 40300)
    #undef  arma_hot
    #undef  arma_cold
    
    #define arma_hot  __attribute__((__hot__))
    #define arma_cold __attribute__((__cold__))
  #endif
  
  #if (ARMA_GCC_VERSION >= 40200)
    #if defined(_GLIBCXX_USE_C99_MATH_TR1) && defined(_GLIBCXX_USE_C99_COMPLEX_TR1)
      #define ARMA_HAVE_STD_TR1
    #endif
  #endif
  
  #if defined(__GXX_EXPERIMENTAL_CXX0X__)
    #undef ARMA_HAVE_STD_TR1
    
    #if !defined(ARMA_USE_CXX11)
      #define ARMA_USE_CXX11
    #endif
  #endif
  
  #if defined(__clang__)
    #undef ARMA_HAVE_STD_TR1
    //#undef ARMA_GOOD_COMPILER
  #endif
  
  #if ( (ARMA_GCC_VERSION >= 40700) && (ARMA_GCC_VERSION <= 40701) )
    #define ARMA_GCC47_BUG
    
    #warning "*** Detected GCC 4.7.0 / 4.7.1, which has a regression (bug)"
    #warning "*** See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53549   "
    #warning "*** A partial workaround for the bug has been activated,    " 
    #warning "*** which reduces some functionality in fixed-size matrices "
  #endif
  
  #undef ARMA_GCC_VERSION
  
#endif


#if defined(__APPLE__)
  #define ARMA_BLAS_SDOT_BUG
#endif


#if defined(_MSC_VER)
  
  #if (_MSC_VER < 1500)
    #error "*** Need a newer compiler ***"
  #endif
  
  #undef ARMA_GOOD_COMPILER
  #undef ARMA_HAVE_STD_ISFINITE
  #undef ARMA_HAVE_STD_SNPRINTF
  #undef ARMA_HAVE_LOG1P
  #undef ARMA_HAVE_STD_ISINF
  #undef ARMA_HAVE_STD_ISNAN
  #undef ARMA_HAVE_STD_TR1
  
  #undef  arma_inline
  #define arma_inline inline __forceinline
  
  // #if (_MSC_VER >= 1400)
  //   #undef  arma_aligned
  //   #define arma_aligned __declspec(align(16))
  // #endif
  
#endif


#if defined(__CUDACC__)
  #undef ARMA_HAVE_STD_ISFINITE
  #undef ARMA_HAVE_STD_SNPRINTF
  #undef ARMA_HAVE_LOG1P
  #undef ARMA_HAVE_STD_ISINF
  #undef ARMA_HAVE_STD_ISNAN
  #undef ARMA_HAVE_STD_TR1
#endif


#if defined(__SUNPRO_CC)
  #undef ARMA_HAVE_STD_ISFINITE
  #undef ARMA_HAVE_STD_SNPRINTF
  #undef ARMA_HAVE_LOG1P
  #undef ARMA_HAVE_STD_ISINF
  #undef ARMA_HAVE_STD_ISNAN
  #undef ARMA_HAVE_STD_TR1
#endif
