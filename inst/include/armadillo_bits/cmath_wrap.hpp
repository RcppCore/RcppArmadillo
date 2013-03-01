// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup cmath_wrap
//! @{



//
// wrappers for isfinite
//



template<typename eT>
arma_inline
bool
arma_isfinite(eT val)
  {
  arma_ignore(val);
    
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



//
// wrappers for trigonometric functions
//



// Wherever possible, try to use TR1 versions of the functions below,
// otherwise fall back to Boost Math.
//
// complex acos
// complex asin
// complex atan
//
// real    acosh
// real    asinh
// real    atanh
//
// complex acosh
// complex asinh
// complex atanh
// 
// 
// If TR1 not present and Boost math not present,
// we have our own rudimentary versions of:
// 
// real    acosh
// real    asinh
// real    atanh



#if defined(ARMA_USE_BOOST)
  #define arma_boost_wrap(trig_fn, val) ( (boost::math::trig_fn)(val) )
#else
  #define arma_boost_wrap(trig_fn, val) ( arma_stop( #trig_fn "(): need Boost libraries" ), val )
#endif


template<typename T>
arma_inline
std::complex<T>
arma_acos(const std::complex<T>& x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::acos(x);
    }
  #else
    {
    return arma_boost_wrap(acos, x);
    }
  #endif
  }



template<typename T>
arma_inline
std::complex<T>
arma_asin(const std::complex<T>& x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::asin(x);
    }
  #else
    {
    return arma_boost_wrap(asin, x);
    }
  #endif
  }



template<typename T>
arma_inline
std::complex<T>
arma_atan(const std::complex<T>& x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::atan(x);
    }
  #else
    {
    return arma_boost_wrap(atan, x);
    }
  #endif
  }



template<typename eT>
arma_inline
eT
arma_acosh(const eT x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::acosh(x);
    }
  #elif defined(ARMA_USE_BOOST)
    {
    return boost::math::acosh(x);
    }
  #else
    {
    if(x >= eT(1))
      {
      // http://functions.wolfram.com/ElementaryFunctions/ArcCosh/02/
      return std::log( x + std::sqrt(x*x - eT(1)) );
      }
    else
      {
      if(std::numeric_limits<eT>::has_quiet_NaN == true)
        {
        return -(std::numeric_limits<eT>::quiet_NaN());
        }
      else
        {
        return eT(0);
        }
      }
    }
  #endif
  }



template<typename eT>
arma_inline
eT
arma_asinh(const eT x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::asinh(x);
    }
  #elif defined(ARMA_USE_BOOST)
    {
    return boost::math::asinh(x);
    }
  #else
    {
    // http://functions.wolfram.com/ElementaryFunctions/ArcSinh/02/
    return std::log( x + std::sqrt(x*x + eT(1)) );
    }
  #endif
  }



template<typename eT>
arma_inline
eT
arma_atanh(const eT x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::atanh(x);
    }
  #elif defined(ARMA_USE_BOOST)
    {
    return boost::math::atanh(x);
    }
  #else
    {
    if( (x >= eT(-1)) && (x <= eT(+1)) )
      {
      // http://functions.wolfram.com/ElementaryFunctions/ArcTanh/02/
      return std::log( ( eT(1)+x ) / ( eT(1)-x ) ) / eT(2);
      }
    else
      {
      if(std::numeric_limits<eT>::has_quiet_NaN == true)
        {
        return -(std::numeric_limits<eT>::quiet_NaN());
        }
      else
        {
        return eT(0);
        }
      }
    }
  #endif
  }



template<typename T>
arma_inline
std::complex<T>
arma_acosh(const std::complex<T>& x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::acosh(x);
    }
  #else
    {
    return arma_boost_wrap(acosh, x);
    }
  #endif
  }



template<typename T>
arma_inline
std::complex<T>
arma_asinh(const std::complex<T>& x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::asinh(x);
    }
  #else
    {
    return arma_boost_wrap(asinh, x);
    }
  #endif
  }



template<typename T>
arma_inline
std::complex<T>
arma_atanh(const std::complex<T>& x)
  {
  #if defined(ARMA_HAVE_STD_TR1)
    {
    return std::tr1::atanh(x);
    }
  #else
    {
    return arma_boost_wrap(atanh, x);
    }
  #endif
  }



#undef arma_boost_wrap



//! @}
