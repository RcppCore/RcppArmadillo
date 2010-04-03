// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
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

#if defined(__GNUG__)

  #if (__GNUC__ < 4)
    #error "*** Need a newer compiler ***"
  #endif
  
  #if (__GNUC_MINOR__ >= 3)
    #undef  arma_hot
    #undef  arma_cold
    
    #define arma_hot  __attribute__((hot))
    #define arma_cold __attribute__((cold))
  #endif

  #undef  arma_pure
  #undef  arma_const
  #undef  arma_inline
  #undef  arma_aligned
  #undef  arma_warn_unused

  #define arma_pure               __attribute__((pure))
  #define arma_const              __attribute__((const))
  #define arma_inline      inline __attribute__((always_inline))
  #define arma_aligned            __attribute__((aligned))
  #define arma_warn_unused        __attribute__((warn_unused_result))
  
  #define ARMA_GOOD_COMPILER
  
#elif defined(__INTEL_COMPILER)
  
  #if (__INTEL_COMPILER < 1000)
    #error "*** Need a newer compiler ***"
  #endif
  
  #define ARMA_GOOD_COMPILER
  
#elif defined(_MSC_VER)
  
  #pragma message ("*** WARNING: This compiler may have an incomplete implementation of the C++ standard ***")
  #undef ARMA_GOOD_COMPILER

#endif


#if defined(__CUDACC__)
  #undef ARMA_HAVE_STD_ISFINITE
  #undef ARMA_HAVE_STD_ISINF
  #undef ARMA_HAVE_STD_ISNAN
#endif


#if defined(__INTEL_COMPILER)
  #if (__INTEL_COMPILER <= 1110)
    #undef ARMA_HAVE_STD_ISFINITE
  #endif
#endif
