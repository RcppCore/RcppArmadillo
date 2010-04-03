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



//! \addtogroup cmath_wrap
//! @{



template<typename eT>
arma_inline
bool
arma_isfinite(eT val)
  {
  return true;
  }



template<>
arma_inline
bool
arma_isfinite(float x)
  {
  #if defined(ARMA_HAVE_STD_ISFINITE)
    {
    return (std::isfinite(x) != 0);
    }
  #else
    {
    const bool x_is_inf = ( (x == x) && ((x - x) != float(0)) );
    const bool x_is_nan = (x != x);

    return ( (x_is_inf == false) && (x_is_nan == false) );
    }
  #endif
  }



template<>
arma_inline
bool
arma_isfinite(double x)
  {
  #if defined(ARMA_HAVE_STD_ISFINITE)
    {
    return (std::isfinite(x) != 0);
    }
  #else
    {
    const bool x_is_inf = ( (x == x) && ((x - x) != double(0)) );
    const bool x_is_nan = (x != x);

    return ( (x_is_inf == false) && (x_is_nan == false) );
    }
  #endif
  }



template<typename T>
arma_inline
bool
arma_isfinite(const std::complex<T>& x)
  {
  if( (arma_isfinite(x.real()) == false) || (arma_isfinite(x.imag()) == false) )
    {
    return false;
    }
  else
    {
    return true;
    }
  }



//! @}
