// Copyright (C) 2008-2013 Conrad Sanderson
// Copyright (C) 2008-2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup typedef
//! @{


#if   UCHAR_MAX >= 0xff
  typedef unsigned char    u8;
  typedef          char    s8;
#elif defined(UINT8_MAX)
  typedef          uint8_t u8;
  typedef           int8_t s8;
#else
  #error "don't know how to typedef 'u8' on this system"
#endif

// NOTE:
// "signed char" is not the same as "char". 
// http://www.embedded.com/columns/programmingpointers/206107018
// http://en.wikipedia.org/wiki/C_variable_types_and_declarations


#if   USHRT_MAX >= 0xffff
  typedef unsigned short    u16;
  typedef          short    s16;
#elif defined(UINT16_MAX)
  typedef          uint16_t u16;
  typedef           int16_t s16;
#else
  #error "don't know how to typedef 'u16' on this system"
#endif


#if   UINT_MAX  >= 0xffffffff
  typedef unsigned int      u32;
  typedef          int      s32;
#elif defined(UINT32_MAX)
  typedef          uint32_t u32;
  typedef           int32_t s32;
#else
  #error "don't know how to typedef 'u32' on this system"
#endif


#if defined(ARMA_USE_U64S64)
  #if   ULLONG_MAX >= 0xffffffffffffffff
    typedef unsigned long long u64;
    typedef          long long s64;
  #elif ULONG_MAX  >= 0xffffffffffffffff
    typedef unsigned long      u64;
    typedef          long      s64;
    #define ARMA_U64_IS_LONG
  #elif defined(UINT64_MAX)
    typedef          uint64_t  u64;
    typedef           int64_t  s64;
  #else
      #error "don't know how to typedef 'u64' on this system; please disable ARMA_64BIT_WORD and/or ARMA_USE_U64S64"
  #endif
#endif


#if !defined(ARMA_USE_U64S64) || (defined(ARMA_USE_U64S64) && !defined(ARMA_U64_IS_LONG))
  #define ARMA_ALLOW_LONG
#endif


typedef unsigned long ulng_t;
typedef          long slng_t;


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

#if defined(ARMA_USE_U64S64)
  typedef Mat <u64> u64_mat;
  typedef Col <u64> u64_vec;
  typedef Col <u64> u64_colvec;
  typedef Row <u64> u64_rowvec;
  typedef Cube<u64> u64_cube;

  typedef Mat <s64> s64_mat;
  typedef Col <s64> s64_vec;
  typedef Col <s64> s64_colvec;
  typedef Row <s64> s64_rowvec;
  typedef Cube<s64> s64_cube;
#endif

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



typedef SpMat <uword> sp_umat;
typedef SpCol <uword> sp_uvec;
typedef SpCol <uword> sp_ucolvec;
typedef SpRow <uword> sp_urowvec;

typedef SpMat <sword> sp_imat;
typedef SpCol <sword> sp_ivec;
typedef SpCol <sword> sp_icolvec;
typedef SpRow <sword> sp_irowvec;

typedef SpMat <float> sp_fmat;
typedef SpCol <float> sp_fvec;
typedef SpCol <float> sp_fcolvec;
typedef SpRow <float> sp_frowvec;

typedef SpMat <double> sp_mat;
typedef SpCol <double> sp_vec;
typedef SpCol <double> sp_colvec;
typedef SpRow <double> sp_rowvec;

typedef SpMat <cx_float> sp_cx_fmat;
typedef SpCol <cx_float> sp_cx_fvec;
typedef SpCol <cx_float> sp_cx_fcolvec;
typedef SpRow <cx_float> sp_cx_frowvec;

typedef SpMat <cx_double> sp_cx_mat;
typedef SpCol <cx_double> sp_cx_vec;
typedef SpCol <cx_double> sp_cx_colvec;
typedef SpRow <cx_double> sp_cx_rowvec;



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
    
    #if defined(ARMA_USE_U64S64)
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
