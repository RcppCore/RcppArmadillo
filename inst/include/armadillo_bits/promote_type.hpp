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


//! \addtogroup promote_type
//! @{



template<typename T1, typename T2>
struct promote_type
  {
  inline static void check()
    {
    arma_static_assert<false> ERROR___UNSUPPORTED_MIXTURE_OF_TYPES;
    ERROR___UNSUPPORTED_MIXTURE_OF_TYPES = ERROR___UNSUPPORTED_MIXTURE_OF_TYPES;
    }
  
  typedef T1 result;
  };



struct promote_type_ok
  {
  arma_inline static void check() {}
  };


template<typename T> struct promote_type<T,               T> : public promote_type_ok { typedef T               result; };
template<typename T> struct promote_type<std::complex<T>, T> : public promote_type_ok { typedef std::complex<T> result; };

template<> struct promote_type<std::complex<double>, std::complex<float> > : public promote_type_ok { typedef std::complex<double> result; };
template<> struct promote_type<std::complex<double>, float>                : public promote_type_ok { typedef std::complex<double> result; };
template<> struct promote_type<std::complex<float>,  double>               : public promote_type_ok { typedef std::complex<double> result; };

template<typename T> struct promote_type<std::complex<T>, s32> : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<std::complex<T>, u32> : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<std::complex<T>, s16> : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<std::complex<T>, u16> : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<std::complex<T>, s8>  : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<std::complex<T>, u8>  : public promote_type_ok { typedef std::complex<T> result; };


template<> struct promote_type<double, float> : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, s32  > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, u32  > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, s16  > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, u16  > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, s8   > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, u8   > : public promote_type_ok { typedef double result; };

template<> struct promote_type<float, s32> : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, u32> : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, s16> : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, u16> : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, s8 > : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, u8 > : public promote_type_ok { typedef float result; };

template<> struct promote_type<s32, u32> : public promote_type_ok { typedef s32 result; };  // float ?  
template<> struct promote_type<s32, s16> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<s32, u16> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<s32, s8 > : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<s32, u8 > : public promote_type_ok { typedef s32 result; };

template<> struct promote_type<u32, s16> : public promote_type_ok { typedef s32 result; };  // float ?
template<> struct promote_type<u32, u16> : public promote_type_ok { typedef u32 result; };
template<> struct promote_type<u32, s8 > : public promote_type_ok { typedef s32 result; };  // float ?
template<> struct promote_type<u32, u8 > : public promote_type_ok { typedef u32 result; };

template<> struct promote_type<s16, u16> : public promote_type_ok { typedef s16 result; };  // s32 ?
template<> struct promote_type<s16, s8 > : public promote_type_ok { typedef s16 result; };
template<> struct promote_type<s16, u8 > : public promote_type_ok { typedef s16 result; };

template<> struct promote_type<u16, s8> : public promote_type_ok { typedef s16 result; };  // s32 ?
template<> struct promote_type<u16, u8> : public promote_type_ok { typedef u16 result; };

template<> struct promote_type<s8, u8> : public promote_type_ok { typedef s8 result; };  // s16 ?




//
// mirrored versions

template<typename T> struct promote_type<T, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };

template<> struct promote_type<std::complex<float>, std::complex<double> > : public promote_type_ok { typedef std::complex<double> result; };
template<> struct promote_type<float,               std::complex<double> > : public promote_type_ok { typedef std::complex<double> result; };
template<> struct promote_type<double,              std::complex<float>  > : public promote_type_ok { typedef std::complex<double> result; };

template<typename T> struct promote_type<s32, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<u32, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<s16, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<u16, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<s8,  std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };
template<typename T> struct promote_type<u8,  std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };


template<> struct promote_type<float, double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<s32  , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<u32  , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<s16  , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<u16  , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<s8   , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<u8   , double> : public promote_type_ok { typedef double result; };

template<> struct promote_type<s32, float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<u32, float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<s16, float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<u16, float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<s8 , float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<u8 , float> : public promote_type_ok { typedef float result; };

template<> struct promote_type<u32, s32> : public promote_type_ok { typedef s32 result; };  // float ?  
template<> struct promote_type<s16, s32> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<u16, s32> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<s8 , s32> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<u8 , s32> : public promote_type_ok { typedef s32 result; };

template<> struct promote_type<s16, u32> : public promote_type_ok { typedef s32 result; };  // float ?
template<> struct promote_type<u16, u32> : public promote_type_ok { typedef u32 result; };
template<> struct promote_type<s8 , u32> : public promote_type_ok { typedef s32 result; };  // float ?
template<> struct promote_type<u8 , u32> : public promote_type_ok { typedef u32 result; };

template<> struct promote_type<u16, s16> : public promote_type_ok { typedef s16 result; };  // s32 ?
template<> struct promote_type<s8 , s16> : public promote_type_ok { typedef s16 result; };
template<> struct promote_type<u8 , s16> : public promote_type_ok { typedef s16 result; };

template<> struct promote_type<s8, u16> : public promote_type_ok { typedef s16 result; };  // s32 ?
template<> struct promote_type<u8, u16> : public promote_type_ok { typedef u16 result; };

template<> struct promote_type<u8, s8> : public promote_type_ok { typedef s8 result; };  // s16 ?

  
  

//! @}
