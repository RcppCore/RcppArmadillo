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


//! \addtogroup version
//! @{



struct arma_version
  {
  static const unsigned int major = 0;
  static const unsigned int minor = 9;
  static const unsigned int patch = 8;
  
  static
  inline
  std::string
  as_string()
    {
    std::stringstream ss;
    ss << arma_version::major << '.' << arma_version::minor << '.' << arma_version::patch;
    
    return ss.str();
    }
  };



struct arma_config
  {
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
  };



//! @}
