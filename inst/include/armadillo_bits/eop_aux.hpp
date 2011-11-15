// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
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
struct eop_aux_randu
  {
  arma_inline
  operator eT ()
    {
    return eT(std::rand()) / eT(RAND_MAX);
    }
  };


  
template<typename T>
struct eop_aux_randu< std::complex<T> >
  {
  arma_inline
  operator std::complex<T> ()
    {
    return std::complex<T>( T(eop_aux_randu<T>()), T(eop_aux_randu<T>()) );
    }
  };



template<typename eT>
struct eop_aux_randn
  {
  // // rudimentary method, based on the central limit theorem
  // // http://en.wikipedia.org/wiki/Central_limit_theorem
  // inline
  // operator eT () const
  //   {
  //   const uword N  = 12;  // N must be >= 12 and an even number
  //   const uword N2 = N/2;
  //   
  //   eT acc = eT(0);
  //   
  //   for(uword i=0; i<N2; ++i)
  //     {
  //     const eT tmp1 = eT(std::rand()) / eT(RAND_MAX);
  //     const eT tmp2 = eT(std::rand()) / eT(RAND_MAX);
  //     acc += tmp1+tmp2;
  //     }
  //   
  //   return acc - eT(N2);
  //   }
  
  
  // polar form of the Box-Muller transformation
  // http://en.wikipedia.org/wiki/Box-Muller_transformation
  // http://en.wikipedia.org/wiki/Marsaglia_polar_method
  inline
  operator eT () const
    {
    // make sure we are internally using at least floats
    typedef typename promote_type<eT,float>::result eTp;
    
    eTp tmp1;
    eTp tmp2;
    eTp w;
    
    do
      {
      tmp1 = eTp(2) * eTp(std::rand()) / eTp(RAND_MAX) - eTp(1);
      tmp2 = eTp(2) * eTp(std::rand()) / eTp(RAND_MAX) - eTp(1);
      w = tmp1*tmp1 + tmp2*tmp2;
      }
    while ( w >= eTp(1) );
    
    return eT( tmp1 * std::sqrt( (eTp(-2) * std::log(w)) / w) );
    }
  
  
  // other methods:
  // http://en.wikipedia.org/wiki/Ziggurat_algorithm
  //
  // Marsaglia and Tsang Ziggurat technique to transform from a uniform to a normal distribution.
  // G. Marsaglia, W.W. Tsang.
  // "Ziggurat method for generating random variables",
  // J. Statistical Software, vol 5, 2000.
  // http://www.jstatsoft.org/v05/i08/
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


//! use of the SFINAE approach to work around compiler limitations
//! http://en.wikipedia.org/wiki/SFINAE

class eop_aux
  {
  public:
  
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result    acos  (const eT x) { return eT( std::acos(double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result    asin  (const eT x) { return eT( std::asin(double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result    atan  (const eT x) { return eT( std::atan(double(x)) ); }
  
  template<typename eT> arma_inline static typename arma_float_only<eT>::result       acos  (const eT x) { return std::acos(x); }
  template<typename eT> arma_inline static typename arma_float_only<eT>::result       asin  (const eT x) { return std::asin(x); }
  template<typename eT> arma_inline static typename arma_float_only<eT>::result       atan  (const eT x) { return std::atan(x); }
  
  template<typename eT> arma_inline static typename arma_cx_only<eT>::result          acos  (const eT x) { return arma_acos(x); }
  template<typename eT> arma_inline static typename arma_cx_only<eT>::result          asin  (const eT x) { return arma_asin(x); }
  template<typename eT> arma_inline static typename arma_cx_only<eT>::result          atan  (const eT x) { return arma_atan(x); }
  
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result    acosh (const eT x) { return eT( arma_acosh(double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result    asinh (const eT x) { return eT( arma_asinh(double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result    atanh (const eT x) { return eT( arma_atanh(double(x)) ); }
  
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result acosh (const eT x) { return arma_acosh(x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result asinh (const eT x) { return arma_asinh(x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result atanh (const eT x) { return arma_atanh(x); }
  
  template<typename eT> arma_inline static typename arma_not_cx<eT>::result conj(const eT               x) { return x;            }
  template<typename  T> arma_inline static          std::complex<T>         conj(const std::complex<T>& x) { return std::conj(x); }
  
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result sqrt  (const eT x) { return eT( std::sqrt (double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result log10 (const eT x) { return eT( std::log10(double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result log   (const eT x) { return eT( std::log  (double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result exp   (const eT x) { return eT( std::exp  (double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result cos   (const eT x) { return eT( std::cos  (double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result sin   (const eT x) { return eT( std::sin  (double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result tan   (const eT x) { return eT( std::tan  (double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result cosh  (const eT x) { return eT( std::cosh (double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result sinh  (const eT x) { return eT( std::sinh (double(x)) ); }
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result tanh  (const eT x) { return eT( std::tanh (double(x)) ); }
  
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result sqrt  (const eT x) { return std::sqrt (x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result log10 (const eT x) { return std::log10(x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result log   (const eT x) { return std::log  (x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result exp   (const eT x) { return std::exp  (x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result cos   (const eT x) { return std::cos  (x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result sin   (const eT x) { return std::sin  (x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result tan   (const eT x) { return std::tan  (x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result cosh  (const eT x) { return std::cosh (x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result sinh  (const eT x) { return std::sinh (x); }
  template<typename eT> arma_inline static typename arma_float_or_cx_only<eT>::result tanh  (const eT x) { return std::tanh (x); }
  
  template<typename eT> arma_inline static typename arma_unsigned_integral_only<eT>::result neg (const eT x) { return  x; }
  template<typename eT> arma_inline static typename arma_signed_only<eT>::result            neg (const eT x) { return -x; }
  
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result floor(const eT  x) { return x;                                                }
  template<typename eT> arma_inline static typename arma_float_only<eT>::result    floor(const eT  x) { return std::floor(x);                                    }
  template<typename eT> arma_inline static typename arma_cx_only<eT>::result       floor(const eT& x) { return eT( std::floor(x.real()), std::floor(x.imag()) ); }
  
  template<typename eT> arma_inline static typename arma_integral_only<eT>::result  ceil(const eT  x) { return x;                                                }
  template<typename eT> arma_inline static typename arma_float_only<eT>::result     ceil(const eT  x) { return std::ceil(x);                                     }
  template<typename eT> arma_inline static typename arma_cx_only<eT>::result        ceil(const eT& x) { return eT( std::ceil(x.real()), std::ceil(x.imag()) );   }
  
  template<typename eT>
  arma_inline
  static
  typename arma_integral_only<eT>::result
  log2 (const eT x)
    {
    return eT( std::log(double(x))/ double(0.69314718055994530942) );
    }
  
  
  template<typename eT>
  arma_inline
  static
  typename arma_float_or_cx_only<eT>::result
  log2 (const eT x)
    {
    typedef typename get_pod_type<eT>::result T;
    return std::log(x) / T(0.69314718055994530942);
    }
  
  
  template<typename eT>
  arma_inline
  static
  typename arma_integral_only<eT>::result
  exp10 (const eT x)
    {
    return eT( std::pow(double(10), double(x)) );
    }
  
  
  template<typename eT>
  arma_inline
  static
  typename
  arma_float_or_cx_only<eT>::result
  exp10 (const eT x)
    {
    typedef typename get_pod_type<eT>::result T;
    return std::pow( T(10), x);
    }
  
  
  template<typename eT>
  arma_inline
  static
  typename arma_integral_only<eT>::result
  exp2 (const eT x)
    {
    return eT( std::pow(double(2), double(x)) );
    }
  
  
  template<typename eT>
  arma_inline
  static
  typename arma_float_or_cx_only<eT>::result
  exp2 (const eT x)
    {
    typedef typename get_pod_type<eT>::result T;
    return std::pow( T(2), x);
    }
  
  
  template<typename T1, typename T2>
  arma_inline
  static
  typename arma_float_or_cx_only<T1>::result
  pow(const T1 base, const T2 exponent)
    {
    return std::pow(base, exponent);
    }
  
  
  
  template<typename T1, typename T2>
  arma_inline
  static
  typename arma_integral_only<T1>::result
  pow(const T1 base, const T2 exponent)
    {
    return T1( std::pow( double(base), double(exponent) ) );
    }
  
  
  
  template<typename eT>
  arma_inline
  static
  typename arma_integral_only<eT>::result
  direct_eps(const eT)
    {
    return eT(0);
    }
  
  
  
  template<typename eT>
  inline
  static
  typename arma_float_only<eT>::result
  direct_eps(const eT x)
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
  direct_eps(const std::complex<T> x)
    {
    //arma_extra_debug_sigprint();
    
    //return std::pow( std::numeric_limits<T>::radix, (std::floor(std::log10(std::abs(x))/std::log10(std::numeric_limits<T>::radix))-(std::numeric_limits<T>::digits-1)) );
    
    const T radix_T     = T(std::numeric_limits<T>::radix);
    const T digits_m1_T = T(std::numeric_limits<T>::digits - 1);
    
    return std::pow( radix_T, T(std::floor(std::log10(std::abs(x))/std::log10(radix_T)) - digits_m1_T) );
    }
  
  
  
  //! work around a bug in GCC 4.4
  template<typename eT> arma_inline static
  typename arma_unsigned_integral_only<eT>::result arma_abs(const eT x)              { return x;           }
  
  template<typename eT> arma_inline static
  typename arma_signed_integral_only<eT>::result   arma_abs(const eT x)              { return std::abs(x); }
  
  template<typename eT> arma_inline static
  typename arma_float_only<eT>::result             arma_abs(const eT x)              { return std::abs(x); }
  
  template<typename T> arma_inline static
  typename arma_float_only<T>::result              arma_abs(const std::complex<T> x) { return std::abs(x); }
  
  };



//! @}

