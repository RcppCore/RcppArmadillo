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


//! \addtogroup eop_aux
//! @{



template<typename eT>
struct eop_aux_rand
  {
  arma_inline
  operator eT ()
    {
    return eT(std::rand()) / eT(RAND_MAX);
    }
  };


  
template<typename T>
struct eop_aux_rand< std::complex<T> >
  {
  arma_inline
  operator std::complex<T> ()
    {
    return std::complex<T>( T(eop_aux_rand<T>()), T(eop_aux_rand<T>()) );
    }
  };



template<typename eT>
struct eop_aux_randn
  {
  // TODO: implement a more efficient method for randn()
  //
  // possible option:
  // Marsaglia and Tsang Ziggurat technique to transform from a uniform to a normal distribution.
  // G. Marsaglia, W.W. Tsang.
  // "Ziggurat method for generating random variables",
  // J. Statistical Software, vol 5, 2000.
  // http://www.jstatsoft.org/v05/i08/
  
  inline
  operator eT () const
    {
    const u32 N  = 12;  // N must be >= 12 and an even number
    const u32 N2 = N/2;
    
    eT acc = eT(0);
    
    for(u32 i=0; i<N2; ++i)
      {
      const eT tmp1 = eT(std::rand()) / eT(RAND_MAX);
      const eT tmp2 = eT(std::rand()) / eT(RAND_MAX);
      acc += tmp1+tmp2;
      }
    
    return acc - eT(N2);
    }
  
  };



template<typename T>
struct eop_aux_randn< std::complex<T> >
  {
  arma_inline
  operator std::complex<T> () const
    {
    return std::complex<T>( T(eop_aux_randn<T>()), T(eop_aux_randn<T>()) );
    }

  };



class eop_aux
  {
  public:
  
  #if defined(ARMA_USE_BOOST)
    #define arma_boost_wrap(trig_fn, val) ( (boost::math::trig_fn)(val) )
  #else
    #define arma_boost_wrap(trig_fn, val) ( arma_stop( #trig_fn "(): need Boost libraries" ), val )
  #endif
  
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result acos (const eT& x) { return std::acos(double(x)); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result asin (const eT& x) { return std::asin(double(x)); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result atan (const eT& x) { return std::atan(double(x)); }
  
  template<typename eT> arma_inline static typename arma_float_only<eT>::result acos (const eT& x) { return std::acos(x); }
  template<typename eT> arma_inline static typename arma_float_only<eT>::result asin (const eT& x) { return std::asin(x); }
  template<typename eT> arma_inline static typename arma_float_only<eT>::result atan (const eT& x) { return std::atan(x); }
  
  template<typename  T> arma_inline static std::complex<T> acos (const std::complex<T>& x) { return arma_boost_wrap(acos,  x); }
  template<typename  T> arma_inline static std::complex<T> asin (const std::complex<T>& x) { return arma_boost_wrap(asin,  x); }
  template<typename  T> arma_inline static std::complex<T> atan (const std::complex<T>& x) { return arma_boost_wrap(atan,  x); }

  template<typename eT> arma_inline static typename arma_integral_only<eT>::result acosh (const eT& x) { return arma_boost_wrap(acosh, double(x)); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result asinh (const eT& x) { return arma_boost_wrap(asinh, double(x)); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result atanh (const eT& x) { return arma_boost_wrap(atanh, double(x)); }

  template<typename eT> arma_inline static typename arma_float_only<eT>::result acosh (const eT& x) { return arma_boost_wrap(acosh, x); }
  template<typename eT> arma_inline static typename arma_float_only<eT>::result asinh (const eT& x) { return arma_boost_wrap(asinh, x); }
  template<typename eT> arma_inline static typename arma_float_only<eT>::result atanh (const eT& x) { return arma_boost_wrap(atanh, x); }
  
  template<typename  T> arma_inline static std::complex<T> acosh (const std::complex<T>& x) { return arma_boost_wrap(acosh, x); }
  template<typename  T> arma_inline static std::complex<T> asinh (const std::complex<T>& x) { return arma_boost_wrap(asinh, x); }
  template<typename  T> arma_inline static std::complex<T> atanh (const std::complex<T>& x) { return arma_boost_wrap(atanh, x); }
  
  #undef arma_boost_wrap
  
  template<typename eT> arma_inline static eT              conj(const eT&              x) { return x;            }
  template<typename  T> arma_inline static std::complex<T> conj(const std::complex<T>& x) { return std::conj(x); }
  
  
  template<typename T1, typename T2> arma_inline static
  T1    pow(const T1    base, const T2 exponent) { return std::pow(base, exponent); }
  
  template<typename T2> arma_inline static
  char  pow(const char  base, const T2 exponent) { typedef char  out_t; return out_t( std::pow(double(base), double(exponent) ) ); }
  
  template<typename T2> arma_inline static
  short pow(const short base, const T2 exponent) { typedef short out_t; return out_t( std::pow(double(base), double(exponent) ) ); }
  
  template<typename T2> arma_inline static
  int   pow(const int   base, const T2 exponent) { typedef int   out_t; return out_t( std::pow(double(base), double(exponent) ) ); }
  
  template<typename T2> arma_inline static
  long  pow(const long  base, const T2 exponent) { typedef long  out_t; return out_t( std::pow(double(base), double(exponent) ) ); }
  
  
  template<typename T2> arma_inline static
  unsigned char  pow(const unsigned char  base, const T2 exponent) { typedef unsigned char  out_t; return out_t( std::pow(double(base), double(exponent) ) ); }
  
  template<typename T2> arma_inline static
  unsigned short pow(const unsigned short base, const T2 exponent) { typedef unsigned short out_t; return out_t( std::pow(double(base), double(exponent) ) ); }
  
  template<typename T2> arma_inline static
  unsigned int   pow(const unsigned int   base, const T2 exponent) { typedef unsigned int   out_t; return out_t( std::pow(double(base), double(exponent) ) ); }
  
  template<typename T2> arma_inline static
  unsigned long  pow(const unsigned long  base, const T2 exponent) { typedef unsigned long  out_t; return out_t( std::pow(double(base), double(exponent) ) ); }

  

  template<typename T1> arma_inline static
  T1 pow_int(const    T1    base, const int exponent) { return std::pow(base, exponent); }
  
  arma_inline static
  char  pow_int(const char  base, const int exponent) { typedef char  out_t; return out_t( std::pow(double(base), exponent) ); }
  
  arma_inline static
  short pow_int(const short base, const int exponent) { typedef short out_t; return out_t( std::pow(double(base), exponent) ); }
  
  arma_inline static
  int   pow_int(const int   base, const int exponent) { typedef int   out_t; return out_t( std::pow(double(base), exponent) ); }
  
  arma_inline static
  long  pow_int(const long  base, const int exponent) { typedef long  out_t; return out_t( std::pow(double(base), exponent) ); }
  
  
  arma_inline static
  unsigned char  pow_int(const unsigned char  base, const int exponent) { typedef unsigned char  out_t; return out_t( std::pow(double(base), exponent) ); }

  arma_inline static
  unsigned short pow_int(const unsigned short base, const int exponent) { typedef unsigned short out_t; return out_t( std::pow(double(base), exponent) ); }
  
  arma_inline static
  unsigned int   pow_int(const unsigned int   base, const int exponent) { typedef unsigned int   out_t; return out_t( std::pow(double(base), exponent) ); }
  
  arma_inline static
  unsigned long  pow_int(const unsigned long  base, const int exponent) { typedef unsigned long  out_t; return out_t( std::pow(double(base), exponent) ); }
  
  
  
  
  template<typename eT>
  inline
  static
  eT
  trunc_exp(const eT x)
    {
    if(std::numeric_limits<eT>::is_iec559 && (x >= Math<eT>::log_max() ))
      {
      return std::numeric_limits<eT>::max();
      }
    else
      {
      return std::exp(x);
      }
    }
  
  
  
  template<typename T>
  inline
  static
  std::complex<T>
  trunc_exp(const std::complex<T>& x)
    {
    return std::exp(x);
    }
  
  
  
  template<typename eT>
  inline 
  static
  eT
  trunc_log(const eT x)
    {
    if(std::numeric_limits<eT>::is_iec559)
      {
      if(x == std::numeric_limits<eT>::infinity())
        {
        return Math<eT>::log_max();
        }
      else
      if(x <= eT(0))
        {
        return Math<eT>::log_min();
        }
      else
        {
        return std::log(x);
        }
      }
    else
      {
      return std::log(x);
      }
    }
  
  
  
  template<typename T>
  inline 
  static
  std::complex<T>
  trunc_log(const std::complex<T>& x)
    {
    return std::log(x);
    }
  
  
  
  template<typename eT>
  arma_inline
  static
  typename arma_integral_only<eT>::result
  direct_eps(const eT& x)
    {
    return eT(0);
    }
  
  
  
  template<typename eT>
  inline
  static
  typename arma_float_only<eT>::result
  direct_eps(const eT& x)
    {
    //arma_extra_debug_sigprint();
    
    // acording to IEEE Standard for Floating-Point Arithmetic (IEEE 754)
    // the mantissa length for double is 53 bits = std::numeric_limits<double>::digits
    // the mantissa length for float  is 24 bits = std::numeric_limits<float >::digits
    
    //return std::pow( std::numeric_limits<eT>::radix, (std::floor(std::log10(std::abs(x))/std::log10(std::numeric_limits<eT>::radix))-(std::numeric_limits<eT>::digits-1)) );
    
    const eT radix_eT     = eT(std::numeric_limits<eT>::radix);
    const eT digits_m1_eT = eT(std::numeric_limits<eT>::digits - 1);
    
    // return std::pow( radix_eT, eT(std::floor(std::log10(std::abs(x))/std::log10(radix_eT)) - digits_m1_eT) );
    return eop_aux::pow( radix_eT, eT(std::floor(std::log10(std::abs(x))/std::log10(radix_eT)) - digits_m1_eT) );
    }
  
  
  
  template<typename T>
  inline
  static
  typename arma_float_only<T>::result
  direct_eps(const std::complex<T>& x)
    {
    //arma_extra_debug_sigprint();
    
    //return std::pow( std::numeric_limits<T>::radix, (std::floor(std::log10(std::abs(x))/std::log10(std::numeric_limits<T>::radix))-(std::numeric_limits<T>::digits-1)) );
    
    const T radix_T     = T(std::numeric_limits<T>::radix);
    const T digits_m1_T = T(std::numeric_limits<T>::digits - 1);
    
    return std::pow( radix_T, T(std::floor(std::log10(std::abs(x))/std::log10(radix_T)) - digits_m1_T) );
    }
  
  
  
  //! work around a bug in GCC 4.4
  template<typename eT> arma_inline static
  typename arma_unsigned_integral_only<eT>::result arma_abs(const eT& x)              { return x;           }
  
  template<typename eT> arma_inline static
  typename arma_signed_integral_only<eT>::result   arma_abs(const eT& x)              { return std::abs(x); }
  
  template<typename eT> arma_inline static
  typename arma_float_only<eT>::result             arma_abs(const eT& x)              { return std::abs(x); }
  
  template<typename T> arma_inline static
  typename arma_float_only<T>::result              arma_abs(const std::complex<T>& x) { return std::abs(x); }
  
  
  
  template<typename eT, typename eop_type>
  arma_inline
  static
  eT
  generate()
    {
         if(is_same_type<eop_type, eop_rand          >::value == true) { return eT(eop_aux_rand<eT>());  }
    else if(is_same_type<eop_type, eop_randn         >::value == true) { return eT(eop_aux_randn<eT>()); }
    else if(is_same_type<eop_type, eop_zeros         >::value == true) { return eT(0);                   }
    else if(is_same_type<eop_type, eop_ones_full     >::value == true) { return eT(1);                   }
    else if(is_same_type<eop_type, eop_cube_rand     >::value == true) { return eT(eop_aux_rand<eT>());  }
    else if(is_same_type<eop_type, eop_cube_randn    >::value == true) { return eT(eop_aux_randn<eT>()); }
    else if(is_same_type<eop_type, eop_cube_zeros    >::value == true) { return eT(0);                   }
    else if(is_same_type<eop_type, eop_cube_ones_full>::value == true) { return eT(1);                   }
    else                                                               { return eT(0);                   }
    }
  
  };



//! @}

