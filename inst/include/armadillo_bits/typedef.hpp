// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup typedef
//! @{


#if UCHAR_MAX >= 0xff
  //! unsigned 8 bit type
  typedef unsigned char u8;
  typedef          char s8;
#else
  #error "don't know how to typedef 'u8' on this system"
#endif

// NOTE:
// "signed char" is not the same as "char". 
// http://www.embedded.com/columns/programmingpointers/206107018
// http://en.wikipedia.org/wiki/C_variable_types_and_declarations


#if USHRT_MAX >= 0xffff
  //! unsigned 16 bit type  
  typedef unsigned short u16;
  typedef          short s16;
#else
  #error "don't know how to typedef 'u16' on this system"
#endif


#if   UINT_MAX  >= 0xffffffff
  typedef unsigned int  u32;
  typedef          int  s32;
#elif ULONG_MAX >= 0xffffffff
  typedef unsigned long u32;
  typedef          long s32;
#else
  #error "don't know how to typedef 'u32' on this system"
#endif



// //
// // only supported by C++0x, via #include <cstdint>, or by C99, via #include <stdint.h>
// 
// //! unsigned 8 bit type
// typedef uint8_t u8;
// 
// //! unsigned 16 bit type  
// typedef uint16_t u16;
// 
// //! unsigned 32 bit type
// typedef uint32_t u32;
//
// //! signed 32 bit type
// typedef int32_t s32;



#if   defined(ARMA_BLAS_LONG_LONG)
  typedef long long blas_int;
#elif defined(ARMA_BLAS_LONG)
  typedef long      blas_int;
#else
  typedef int       blas_int;
#endif



typedef std::complex<float>  cx_float;
typedef std::complex<double> cx_double;

typedef Mat<unsigned char>  uchar_mat;
typedef Col<unsigned char>  uchar_vec;
typedef Col<unsigned char>  uchar_colvec;
typedef Row<unsigned char>  uchar_rowvec;
typedef Cube<unsigned char> uchar_cube;

typedef Mat<u32>  umat;
typedef Col<u32>  uvec;
typedef Col<u32>  ucolvec;
typedef Row<u32>  urowvec;
typedef Cube<u32> ucube;

typedef Mat<s32>  imat;
typedef Col<s32>  ivec;
typedef Col<s32>  icolvec;
typedef Row<s32>  irowvec;
typedef Cube<s32> icube;

typedef Mat<float>  fmat;
typedef Col<float>  fvec;
typedef Col<float>  fcolvec;
typedef Row<float>  frowvec;
typedef Cube<float> fcube;

typedef Mat<double>  mat;
typedef Col<double>  vec;
typedef Col<double>  colvec;
typedef Row<double>  rowvec;
typedef Cube<double> cube;

typedef Mat<cx_float>  cx_fmat;
typedef Col<cx_float>  cx_fvec;
typedef Col<cx_float>  cx_fcolvec;
typedef Row<cx_float>  cx_frowvec;
typedef Cube<cx_float> cx_fcube;

typedef Mat<cx_double>  cx_mat;
typedef Col<cx_double>  cx_vec;
typedef Col<cx_double>  cx_colvec;
typedef Row<cx_double>  cx_rowvec;
typedef Cube<cx_double> cx_cube;



namespace junk
  {
  struct arma_elem_size_test
    {
  
    arma_static_assert<sizeof(u8) == 1> ERROR___TYPE_U8_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(s8) == 1> ERROR___TYPE_S8_HAS_UNSUPPORTED_SIZE;
    
    arma_static_assert<sizeof(u16) == 2> ERROR___TYPE_U16_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(s16) == 2> ERROR___TYPE_S16_HAS_UNSUPPORTED_SIZE;
    
    arma_static_assert<sizeof(u32) == 4> ERROR___TYPE_U32_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(s32) == 4> ERROR___TYPE_S32_HAS_UNSUPPORTED_SIZE;
    
    arma_static_assert<sizeof(float)  == 4> ERROR___TYPE_FLOAT_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(double) == 8> ERROR___TYPE_DOUBLE_HAS_UNSUPPORTED_SIZE;
    
    arma_static_assert<sizeof(std::complex<float>)  == 8>  ERROR___TYPE_COMPLEX_FLOAT_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(std::complex<double>) == 16> ERROR___TYPE_COMPLEX_DOUBLE_HAS_UNSUPPORTED_SIZE;
  
    };
  }


//! @}
