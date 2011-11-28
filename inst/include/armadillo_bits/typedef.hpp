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


#if defined(ARMA_64BIT_WORD)
  #if    ULONG_MAX >= 0xffffffffffffffff
    typedef unsigned long      u64;
    typedef          long      s64;
  #else
    #if ULLONG_MAX >= 0xffffffffffffffff
      typedef unsigned long long u64;
      typedef          long long s64;
    #else
      #error "don't know how to typedef 'u64' on this system"
    #endif
  #endif
#endif



// // only supported by C++11, via #include <cstdint>, or by C99, via #include <stdint.h>
// 
// typedef  uint8_t u8;
// typedef   int8_t s8;
// 
// typedef uint16_t u16;
// typedef  int16_t s16;
// 
// typedef uint32_t u32;
// typedef  int32_t s32;
// 
// typedef uint64_t u64;
// typedef  int64_t s64;



#if !defined(ARMA_64BIT_WORD)
  typedef u32 uword;
  typedef s32 sword;

  typedef u16 uhword;
  typedef s16 shword;
  
  #define ARMA_MAX_UWORD  0xffffffff
  #define ARMA_MAX_UHWORD 0xffff
#else
  typedef u64 uword;
  typedef s64 sword;
  
  typedef u32 uhword;
  typedef s32 shword;

  #define ARMA_MAX_UWORD  0xffffffffffffffff
  #define ARMA_MAX_UHWORD 0xffffffff
#endif



typedef std::complex<float>  cx_float;
typedef std::complex<double> cx_double;

typedef Mat <unsigned char> uchar_mat;
typedef Col <unsigned char> uchar_vec;
typedef Col <unsigned char> uchar_colvec;
typedef Row <unsigned char> uchar_rowvec;
typedef Cube<unsigned char> uchar_cube;

typedef Mat <u32> u32_mat;
typedef Col <u32> u32_vec;
typedef Col <u32> u32_colvec;
typedef Row <u32> u32_rowvec;
typedef Cube<u32> u32_cube;

typedef Mat <s32> s32_mat;
typedef Col <s32> s32_vec;
typedef Col <s32> s32_colvec;
typedef Row <s32> s32_rowvec;
typedef Cube<s32> s32_cube;

typedef Mat <uword> umat;
typedef Col <uword> uvec;
typedef Col <uword> ucolvec;
typedef Row <uword> urowvec;
typedef Cube<uword> ucube;

typedef Mat <sword> imat;
typedef Col <sword> ivec;
typedef Col <sword> icolvec;
typedef Row <sword> irowvec;
typedef Cube<sword> icube;

typedef Mat <float> fmat;
typedef Col <float> fvec;
typedef Col <float> fcolvec;
typedef Row <float> frowvec;
typedef Cube<float> fcube;

typedef Mat <double> mat;
typedef Col <double> vec;
typedef Col <double> colvec;
typedef Row <double> rowvec;
typedef Cube<double> cube;

typedef Mat <cx_float> cx_fmat;
typedef Col <cx_float> cx_fvec;
typedef Col <cx_float> cx_fcolvec;
typedef Row <cx_float> cx_frowvec;
typedef Cube<cx_float> cx_fcube;

typedef Mat <cx_double> cx_mat;
typedef Col <cx_double> cx_vec;
typedef Col <cx_double> cx_colvec;
typedef Row <cx_double> cx_rowvec;
typedef Cube<cx_double> cx_cube;



typedef void* void_ptr;



namespace junk
  {
  struct arma_elem_size_test
    {
    
    arma_static_check( (sizeof(u8) != 1), ERROR___TYPE_U8_HAS_UNSUPPORTED_SIZE );
    arma_static_check( (sizeof(s8) != 1), ERROR___TYPE_S8_HAS_UNSUPPORTED_SIZE );
    
    arma_static_check( (sizeof(u16) != 2), ERROR___TYPE_U16_HAS_UNSUPPORTED_SIZE );
    arma_static_check( (sizeof(s16) != 2), ERROR___TYPE_S16_HAS_UNSUPPORTED_SIZE );
    
    arma_static_check( (sizeof(u32) != 4), ERROR___TYPE_U32_HAS_UNSUPPORTED_SIZE );
    arma_static_check( (sizeof(s32) != 4), ERROR___TYPE_S32_HAS_UNSUPPORTED_SIZE );
    
    #if defined(ARMA_64BIT_WORD)
    arma_static_check( (sizeof(u64) != 8), ERROR___TYPE_U64_HAS_UNSUPPORTED_SIZE );
    arma_static_check( (sizeof(s64) != 8), ERROR___TYPE_S64_HAS_UNSUPPORTED_SIZE );
    #endif
    
    arma_static_check( (sizeof(float)  != 4), ERROR___TYPE_FLOAT_HAS_UNSUPPORTED_SIZE );
    arma_static_check( (sizeof(double) != 8), ERROR___TYPE_DOUBLE_HAS_UNSUPPORTED_SIZE );
    
    arma_static_check( (sizeof(std::complex<float>)  != 8),  ERROR___TYPE_COMPLEX_FLOAT_HAS_UNSUPPORTED_SIZE );
    arma_static_check( (sizeof(std::complex<double>) != 16), ERROR___TYPE_COMPLEX_DOUBLE_HAS_UNSUPPORTED_SIZE );
    
    };
  }


//! @}
