// Copyright (C) 2008-2013 Conrad Sanderson
// Copyright (C) 2008-2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



#define arma_hot
#define arma_cold
#define arma_pure
#define arma_const
#define arma_aligned
#define arma_align_mem
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


// posix_memalign() is part of IEEE standard 1003.1
// http://pubs.opengroup.org/onlinepubs/009696899/functions/posix_memalign.html
// http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/unistd.h.html
// http://sourceforge.net/p/predef/wiki/Standards/
#if ( defined(_POSIX_ADVISORY_INFO) && (_POSIX_ADVISORY_INFO >= 200112L) )
  #define ARMA_HAVE_POSIX_MEMALIGN
#endif


#if defined(__APPLE__)
  #define ARMA_BLAS_SDOT_BUG
  #undef  ARMA_HAVE_POSIX_MEMALIGN
#endif


#if defined(__MINGW32__)
  #undef ARMA_HAVE_POSIX_MEMALIGN
#endif


#if (__cplusplus >= 201103L)
  #if !defined(ARMA_USE_CXX11)
    #define ARMA_USE_CXX11
  #endif
#endif


#if defined(ARMA_USE_CXX11)
  
  // TODO: enable in version 4.0
  // #if defined(ARMA_USE_WRAPPER)
  //   #if !defined(ARMA_USE_CXX11_RNG)
  //     #define ARMA_USE_CXX11_RNG
  //   #endif
  // #endif
  
  #if !defined(ARMA_USE_U64S64)
    #define ARMA_USE_U64S64
  #endif
#endif


#if defined(ARMA_64BIT_WORD)
  #if !defined(ARMA_USE_U64S64)
    #define ARMA_USE_U64S64
  #endif
#endif


#if defined (__GNUG__)
  #define ARMA_FNSIG  __PRETTY_FUNCTION__
#elif defined (_MSC_VER)
  #define ARMA_FNSIG  __FUNCSIG__ 
#elif defined(__INTEL_COMPILER)
  #define ARMA_FNSIG  __FUNCTION__
#elif defined(ARMA_USE_CXX11)
  #define ARMA_FNSIG  __func__
#else 
  #define ARMA_FNSIG  "(unknown)"
#endif


#if defined(__INTEL_COMPILER)
  
  #if (__INTEL_COMPILER < 1000)
    #error "*** Need a newer compiler ***"
  #endif
  
  #if (__INTEL_COMPILER <= 1110)
    #undef ARMA_HAVE_STD_ISFINITE
  #endif
  
  #undef  ARMA_HAVE_STD_TR1
  
  #define ARMA_GOOD_COMPILER
  #define ARMA_HAVE_ICC_ASSUME_ALIGNED
  
  #if defined(__GNUG__)
    
    // #undef  arma_aligned
    // #define arma_aligned __attribute__((aligned))
    
    #undef  arma_align_mem
    #define arma_align_mem __attribute__((aligned(16)))
    
    #define ARMA_HAVE_ALIGNED_ATTRIBUTE
    
  #elif defined(_MSC_VER)
    
    // #if (_MANAGED == 1) || (_M_CEE == 1)
    //   
    //   // don't do any alignment when compiling in "managed code" mode 
    //   
    //   #undef  arma_aligned
    //   #define arma_aligned
    //   
    //   #undef  arma_align_mem
    //   #define arma_align_mem
    //   
    // #elif (_MSC_VER >= 1700)
    //   
    //   #undef  arma_align_mem
    //   #define arma_align_mem __declspec(align(16))
    //   
    //   #define ARMA_HAVE_ALIGNED_ATTRIBUTE
    //   
    // #endif
    
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
  #undef  arma_align_mem
  #undef  arma_warn_unused
  #undef  arma_deprecated
  #undef  arma_malloc
  #undef  arma_inline
  #undef  arma_noinline
  
  #define arma_pure               __attribute__((__pure__))
  #define arma_const              __attribute__((__const__))
  #define arma_aligned            __attribute__((__aligned__))
  #define arma_align_mem          __attribute__((__aligned__(16)))
  #define arma_warn_unused        __attribute__((__warn_unused_result__))
  #define arma_deprecated         __attribute__((__deprecated__))
  #define arma_malloc             __attribute__((__malloc__))
  #define arma_inline      inline __attribute__((__always_inline__))
  #define arma_noinline           __attribute__((__noinline__))
  
  #define ARMA_HAVE_ALIGNED_ATTRIBUTE
  
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
    #error "gcc versions 4.7.0 and 4.7.1 are unsupported; use 4.7.2 or later"
    // due to http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53549
  #endif
  
  #if (ARMA_GCC_VERSION >= 40700)
    #define ARMA_HAVE_GCC_ASSUME_ALIGNED
    // TODO: future versions of clang may also have __builtin_assume_aligned
  #endif
  
  #undef ARMA_GCC_VERSION
  
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
  
  // #undef  arma_inline
  // #define arma_inline inline __forceinline
  
  #pragma warning(push)
  
  #pragma warning(disable: 4127)  // conditional expression is constant
  #pragma warning(disable: 4510)  // default constructor could not be generated
  #pragma warning(disable: 4511)  // copy constructor can't be generated
  #pragma warning(disable: 4512)  // assignment operator can't be generated
  #pragma warning(disable: 4513)  // destructor can't be generated
  #pragma warning(disable: 4514)  // unreferenced inline function has been removed
  #pragma warning(disable: 4522)  // multiple assignment operators specified
  #pragma warning(disable: 4623)  // default constructor can't be generated
  #pragma warning(disable: 4624)  // destructor can't be generated
  #pragma warning(disable: 4625)  // copy constructor can't be generated
  #pragma warning(disable: 4626)  // assignment operator can't be generated
  #pragma warning(disable: 4710)  // function not inlined
  #pragma warning(disable: 4711)  // call was inlined
  #pragma warning(disable: 4714)  // __forceinline can't be inlined
  
  // #if (_MANAGED == 1) || (_M_CEE == 1)
  //   
  //   // don't do any alignment when compiling in "managed code" mode 
  //   
  //   #undef  arma_aligned
  //   #define arma_aligned
  //   
  //   #undef  arma_align_mem
  //   #define arma_align_mem
  //  
  // #elif (_MSC_VER >= 1700)
  //   
  //   #undef  arma_align_mem
  //   #define arma_align_mem __declspec(align(16))
  //   
  //   #define ARMA_HAVE_ALIGNED_ATTRIBUTE
  //   
  //   // disable warnings: "structure was padded due to __declspec(align(16))"
  //   #pragma warning(disable: 4324)
  //   
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
