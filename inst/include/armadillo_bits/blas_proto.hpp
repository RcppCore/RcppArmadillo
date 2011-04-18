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



#ifdef ARMA_USE_BLAS


#if !defined(ARMA_BLAS_CAPITALS)

  #define arma_sdot  sdot
  #define arma_ddot  ddot

  #define arma_sgemv sgemv
  #define arma_dgemv dgemv
  #define arma_cgemv cgemv
  #define arma_zgemv zgemv

  #define arma_sgemm sgemm
  #define arma_dgemm dgemm
  #define arma_cgemm cgemm
  #define arma_zgemm zgemm

#else

  #define arma_sdot  SDOT
  #define arma_ddot  DDOT

  #define arma_sgemv SGEMV
  #define arma_dgemv DGEMV
  #define arma_cgemv CGEMV
  #define arma_zgemv ZGEMV

  #define arma_sgemm SGEMM
  #define arma_dgemm DGEMM
  #define arma_cgemm CGEMM
  #define arma_zgemm ZGEMM

#endif



//! \namespace blas namespace for BLAS functions
namespace blas
  {
  extern "C"
    {
    float  arma_fortran(arma_sdot)(blas_int* n, const float*  x, blas_int* incx, const float*  y, blas_int* incy);
    double arma_fortran(arma_ddot)(blas_int* n, const double* x, blas_int* incx, const double* y, blas_int* incy);
    
    void arma_fortran(arma_sgemv)(const char* transA, const blas_int* m, const blas_int* n, const float*  alpha, const float*  A, const blas_int* ldA, const float*  x, const blas_int* incx, const float*  beta, float*  y, const blas_int* incy);
    void arma_fortran(arma_dgemv)(const char* transA, const blas_int* m, const blas_int* n, const double* alpha, const double* A, const blas_int* ldA, const double* x, const blas_int* incx, const double* beta, double* y, const blas_int* incy);
    void arma_fortran(arma_cgemv)(const char* transA, const blas_int* m, const blas_int* n, const void*   alpha, const void*   A, const blas_int* ldA, const void*   x, const blas_int* incx, const void*   beta, void*   y, const blas_int* incy);
    void arma_fortran(arma_zgemv)(const char* transA, const blas_int* m, const blas_int* n, const void*   alpha, const void*   A, const blas_int* ldA, const void*   x, const blas_int* incx, const void*   beta, void*   y, const blas_int* incy);
    
    void arma_fortran(arma_sgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const float*  alpha, const float*  A, const blas_int* ldA, const float*  B, const blas_int* ldB, const float*  beta, float*  C, const blas_int* ldC);
    void arma_fortran(arma_dgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const double* alpha, const double* A, const blas_int* ldA, const double* B, const blas_int* ldB, const double* beta, double* C, const blas_int* ldC);
    void arma_fortran(arma_cgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const void*   alpha, const void*   A, const blas_int* ldA, const void*   B, const blas_int* ldB, const void*   beta, void*   C, const blas_int* ldC);
    void arma_fortran(arma_zgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const void*   alpha, const void*   A, const blas_int* ldA, const void*   B, const blas_int* ldB, const void*   beta, void*   C, const blas_int* ldC);
    
    // void   arma_fortran(arma_dswap)(const blas_int* n, double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void   arma_fortran(arma_dscal)(const blas_int* n, const double* alpha, double* x, const blas_int* incx);
    // void   arma_fortran(arma_dcopy)(const blas_int* n, const double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void   arma_fortran(arma_daxpy)(const blas_int* n, const double* alpha, const double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void   arma_fortran(arma_dger )(const blas_int* m, const blas_int* n, const double* alpha, const double* x, const blas_int* incx, const double* y, const blas_int* incy, double* A, const blas_int* ldA);
    }
  
  
  
  template<typename eT>
  inline
  eT
  dot(const u32 n_elem, const eT* x, const eT* y)
    {
    arma_ignore(n_elem);
    arma_ignore(x);
    arma_ignore(y);
    
    return eT(0);
    }
  
  
  
  template<>
  inline
  float
  dot(const u32 n_elem, const float* x, const float* y)
    {
    blas_int n   = blas_int(n_elem);
    blas_int inc = blas_int(1);
    
    return arma_fortran(arma_sdot)(&n, x, &inc, y, &inc);
    }
  
  
  
  template<>
  inline
  double
  dot(const u32 n_elem, const double* x, const double* y)
    {
    blas_int n   = blas_int(n_elem);
    blas_int inc = blas_int(1);
    
    return arma_fortran(arma_ddot)(&n, x, &inc, y, &inc);
    }
  
  
  
  template<typename eT>
  inline
  void
  gemv(const char* transA, const blas_int* m, const blas_int* n, const eT* alpha, const eT* A, const blas_int* ldA, const eT* x, const blas_int* incx, const eT* beta, eT* y, const blas_int* incy)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(arma_sgemv)(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgemv)(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cgemv)(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zgemv)(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    
    }
  
  
  
  template<typename eT>
  inline
  void
  gemm(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const eT* alpha, const eT* A, const blas_int* ldA, const eT* B, const blas_int* ldB, const eT* beta, eT* C, const blas_int* ldC)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(arma_sgemm)(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgemm)(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cgemm)(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zgemm)(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    
    }
  
  }


#endif
