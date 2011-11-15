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


//! \addtogroup arma_config
//! @{



struct arma_config
  {
  #if defined(ARMA_MAT_PREALLOC)
    static const uword mat_prealloc = (sword(ARMA_MAT_PREALLOC) > 0) ? uword(ARMA_MAT_PREALLOC) : 1;
  #else
    static const uword mat_prealloc = 16;
  #endif
  
  #if defined(ARMA_USE_ATLAS)
    static const bool atlas = true;
  #else
    static const bool atlas = false;
  #endif
  
  
  #if defined(ARMA_USE_LAPACK)
    static const bool lapack = true;
  #else
    static const bool lapack = false;
  #endif
  
  
  #if defined(ARMA_USE_BLAS)
    static const bool blas = true;
  #else
    static const bool blas = false;
  #endif


  #if defined(ARMA_USE_BOOST)
    static const bool boost = true;
  #else
    static const bool boost = false;
  #endif
  

  #if defined(ARMA_USE_BOOST_DATE)
    static const bool boost_date = true;
  #else
    static const bool boost_date = false;
  #endif


  #if !defined(ARMA_NO_DEBUG) && !defined(NDEBUG)
    static const bool debug = true;
  #else
    static const bool debug = false;
  #endif
  
  
  #if defined(ARMA_EXTRA_DEBUG)
    static const bool extra_debug = true;
  #else
    static const bool extra_debug = false;
  #endif
  
  
  #if defined(ARMA_GOOD_COMPILER)
    static const bool good_comp = true;
  #else
    static const bool good_comp = false;
  #endif
  
  
  #if (  \
         defined(ARMA_EXTRA_MAT_PROTO)   || defined(ARMA_EXTRA_MAT_MEAT)   \
      || defined(ARMA_EXTRA_COL_PROTO)   || defined(ARMA_EXTRA_COL_MEAT)   \
      || defined(ARMA_EXTRA_ROW_PROTO)   || defined(ARMA_EXTRA_ROW_MEAT)   \
      || defined(ARMA_EXTRA_CUBE_PROTO)  || defined(ARMA_EXTRA_CUBE_MEAT)  \
      || defined(ARMA_EXTRA_FIELD_PROTO) || defined(ARMA_EXTRA_FIELD_MEAT) \
      )
    static const bool extra_code = true;
  #else
    static const bool extra_code = false;
  #endif
  };



//! @}
