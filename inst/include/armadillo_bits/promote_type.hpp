// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
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
struct is_promotable
  {
  static const bool value = false;
  typedef T1 result;
  };


struct is_promotable_ok
  {
  static const bool value = true;
  };


template<typename T> struct is_promotable<T,               T> : public is_promotable_ok { typedef T               result; };
template<typename T> struct is_promotable<std::complex<T>, T> : public is_promotable_ok { typedef std::complex<T> result; };

template<> struct is_promotable<std::complex<double>, std::complex<float> > : public is_promotable_ok { typedef std::complex<double> result; };
template<> struct is_promotable<std::complex<double>, float>                : public is_promotable_ok { typedef std::complex<double> result; };
template<> struct is_promotable<std::complex<float>,  double>               : public is_promotable_ok { typedef std::complex<double> result; };

// template<typename T> struct is_promotable<std::complex<T>, u64> : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<std::complex<T>, s32> : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<std::complex<T>, u32> : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<std::complex<T>, s16> : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<std::complex<T>, u16> : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<std::complex<T>, s8>  : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<std::complex<T>, u8>  : public is_promotable_ok { typedef std::complex<T> result; };


template<> struct is_promotable<double, float> : public is_promotable_ok { typedef double result; };
// template<> struct is_promotable<double, u64  > : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<double, s32  > : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<double, u32  > : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<double, s16  > : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<double, u16  > : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<double, s8   > : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<double, u8   > : public is_promotable_ok { typedef double result; };

// template<> struct is_promotable<float, u64> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<float, s32> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<float, u32> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<float, s16> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<float, u16> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<float, s8 > : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<float, u8 > : public is_promotable_ok { typedef float result; };

// template<> struct is_promotable<u64, u32> : public is_promotable_ok { typedef u64 result; };
// template<> struct is_promotable<u64, u16> : public is_promotable_ok { typedef u64 result; };
// template<> struct is_promotable<u64, u8 > : public is_promotable_ok { typedef u64 result; };

template<> struct is_promotable<s32, u32> : public is_promotable_ok { typedef s32 result; };  // float ?  
template<> struct is_promotable<s32, s16> : public is_promotable_ok { typedef s32 result; };
template<> struct is_promotable<s32, u16> : public is_promotable_ok { typedef s32 result; };
template<> struct is_promotable<s32, s8 > : public is_promotable_ok { typedef s32 result; };
template<> struct is_promotable<s32, u8 > : public is_promotable_ok { typedef s32 result; };

template<> struct is_promotable<u32, s16> : public is_promotable_ok { typedef s32 result; };  // float ?
template<> struct is_promotable<u32, u16> : public is_promotable_ok { typedef u32 result; };
template<> struct is_promotable<u32, s8 > : public is_promotable_ok { typedef s32 result; };  // float ?
template<> struct is_promotable<u32, u8 > : public is_promotable_ok { typedef u32 result; };

template<> struct is_promotable<s16, u16> : public is_promotable_ok { typedef s16 result; };  // s32 ?
template<> struct is_promotable<s16, s8 > : public is_promotable_ok { typedef s16 result; };
template<> struct is_promotable<s16, u8 > : public is_promotable_ok { typedef s16 result; };

template<> struct is_promotable<u16, s8> : public is_promotable_ok { typedef s16 result; };  // s32 ?
template<> struct is_promotable<u16, u8> : public is_promotable_ok { typedef u16 result; };

template<> struct is_promotable<s8, u8> : public is_promotable_ok { typedef s8 result; };  // s16 ?




//
// mirrored versions

template<typename T> struct is_promotable<T, std::complex<T> > : public is_promotable_ok { typedef std::complex<T> result; };

template<> struct is_promotable<std::complex<float>, std::complex<double> > : public is_promotable_ok { typedef std::complex<double> result; };
template<> struct is_promotable<float,               std::complex<double> > : public is_promotable_ok { typedef std::complex<double> result; };
template<> struct is_promotable<double,              std::complex<float>  > : public is_promotable_ok { typedef std::complex<double> result; };

// template<typename T> struct is_promotable<u64, std::complex<T> > : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<s32, std::complex<T> > : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<u32, std::complex<T> > : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<s16, std::complex<T> > : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<u16, std::complex<T> > : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<s8,  std::complex<T> > : public is_promotable_ok { typedef std::complex<T> result; };
template<typename T> struct is_promotable<u8,  std::complex<T> > : public is_promotable_ok { typedef std::complex<T> result; };


template<> struct is_promotable<float, double> : public is_promotable_ok { typedef double result; };
// template<> struct is_promotable<u64  , double> : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<s32  , double> : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<u32  , double> : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<s16  , double> : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<u16  , double> : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<s8   , double> : public is_promotable_ok { typedef double result; };
template<> struct is_promotable<u8   , double> : public is_promotable_ok { typedef double result; };

// template<> struct is_promotable<u64, float> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<s32, float> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<u32, float> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<s16, float> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<u16, float> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<s8 , float> : public is_promotable_ok { typedef float result; };
template<> struct is_promotable<u8 , float> : public is_promotable_ok { typedef float result; };

// template<> struct is_promotable<u32, u64> : public is_promotable_ok { typedef u64 result; };
// template<> struct is_promotable<u16, u64> : public is_promotable_ok { typedef u64 result; };
// template<> struct is_promotable<u8,  u64> : public is_promotable_ok { typedef u64 result; };

template<> struct is_promotable<u32, s32> : public is_promotable_ok { typedef s32 result; };  // float ?  
template<> struct is_promotable<s16, s32> : public is_promotable_ok { typedef s32 result; };
template<> struct is_promotable<u16, s32> : public is_promotable_ok { typedef s32 result; };
template<> struct is_promotable<s8 , s32> : public is_promotable_ok { typedef s32 result; };
template<> struct is_promotable<u8 , s32> : public is_promotable_ok { typedef s32 result; };

template<> struct is_promotable<s16, u32> : public is_promotable_ok { typedef s32 result; };  // float ?
template<> struct is_promotable<u16, u32> : public is_promotable_ok { typedef u32 result; };
template<> struct is_promotable<s8 , u32> : public is_promotable_ok { typedef s32 result; };  // float ?
template<> struct is_promotable<u8 , u32> : public is_promotable_ok { typedef u32 result; };

template<> struct is_promotable<u16, s16> : public is_promotable_ok { typedef s16 result; };  // s32 ?
template<> struct is_promotable<s8 , s16> : public is_promotable_ok { typedef s16 result; };
template<> struct is_promotable<u8 , s16> : public is_promotable_ok { typedef s16 result; };

template<> struct is_promotable<s8, u16> : public is_promotable_ok { typedef s16 result; };  // s32 ?
template<> struct is_promotable<u8, u16> : public is_promotable_ok { typedef u16 result; };

template<> struct is_promotable<u8, s8> : public is_promotable_ok { typedef s8 result; };  // s16 ?





template<typename T1, typename T2>
struct promote_type
  {
  inline static void check()
    {
    arma_static_assert< is_promotable<T1,T2>::value > ERROR___UNSUPPORTED_MIXTURE_OF_TYPES;
    ERROR___UNSUPPORTED_MIXTURE_OF_TYPES = ERROR___UNSUPPORTED_MIXTURE_OF_TYPES;
    }
  
  typedef typename is_promotable<T1,T2>::result result;
  };



template<typename T1, typename T2>
struct eT_promoter
  {
  typedef typename promote_type<typename T1::elem_type, typename T2::elem_type>::result eT;
  };



//! @}
