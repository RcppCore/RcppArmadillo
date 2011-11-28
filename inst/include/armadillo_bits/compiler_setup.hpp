// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



#define arma_hot
#define arma_cold
#define arma_pure
#define arma_const
#define arma_inline  inline
#define arma_aligned
#define arma_warn_unused
#define arma_deprecated
#define arma_ignore(variable)  ((void)(variable))


#if defined(ARMA_BLAS_UNDERSCORE)
  #define arma_fortran2_noprefix(function) function##_
  #define arma_fortran2_prefix(function)   wrapper_##function##_
#else
  #define arma_fortran2_prefix(function)   wrapper_##function
  #define arma_fortran2_noprefix(function) function
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
  #undef  arma_inline
  #undef  arma_aligned
  #undef  arma_warn_unused
  #undef  arma_deprecated
  
  #define arma_pure               __attribute__((pure))
  #define arma_const              __attribute__((const))
  #define arma_inline      inline __attribute__((always_inline))
  #define arma_aligned            __attribute__((aligned))
  #define arma_warn_unused        __attribute__((warn_unused_result))
  #define arma_deprecated         __attribute__((deprecated))
  
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
  #endif
  
  #if (ARMA_GCC_VERSION >= 40300)
    #undef  arma_hot
    #undef  arma_cold
    
    #define arma_hot  __attribute__((hot))
    #define arma_cold __attribute__((cold))
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
