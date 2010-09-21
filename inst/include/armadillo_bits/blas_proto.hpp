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


#ifdef ARMA_USE_BLAS

//! \namespace blas namespace for BLAS functions
namespace blas
  {
  extern "C"
    {
    float  sdot_(const blas_int* n, const float*  x, const blas_int* incx, const float*  y, const blas_int* incy);
    double ddot_(const blas_int* n, const double* x, const blas_int* incx, const double* y, const blas_int* incy);
    
    void sgemv_(const char* transA, const blas_int* m, const blas_int* n, const float*  alpha, const float*  A, const blas_int* ldA, const float*  x, const blas_int* incx, const float*  beta, float*  y, const blas_int* incy);
    void dgemv_(const char* transA, const blas_int* m, const blas_int* n, const double* alpha, const double* A, const blas_int* ldA, const double* x, const blas_int* incx, const double* beta, double* y, const blas_int* incy);
    void cgemv_(const char* transA, const blas_int* m, const blas_int* n, const void*   alpha, const void*   A, const blas_int* ldA, const void*   x, const blas_int* incx, const void*   beta, void*   y, const blas_int* incy);
    void zgemv_(const char* transA, const blas_int* m, const blas_int* n, const void*   alpha, const void*   A, const blas_int* ldA, const void*   x, const blas_int* incx, const void*   beta, void*   y, const blas_int* incy);
    
    void sgemm_(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const float*  alpha, const float*  A, const blas_int* ldA, const float*  B, const blas_int* ldB, const float*  beta, float*  C, const blas_int* ldC);
    void dgemm_(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const double* alpha, const double* A, const blas_int* ldA, const double* B, const blas_int* ldB, const double* beta, double* C, const blas_int* ldC);
    void cgemm_(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const void*   alpha, const void*   A, const blas_int* ldA, const void*   B, const blas_int* ldB, const void*   beta, void*   C, const blas_int* ldC);
    void zgemm_(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const void*   alpha, const void*   A, const blas_int* ldA, const void*   B, const blas_int* ldB, const void*   beta, void*   C, const blas_int* ldC);
    
    // void   dswap_(const blas_int* n, double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void   dscal_(const blas_int* n, const double* alpha, double* x, const blas_int* incx);
    // void   dcopy_(const blas_int* n, const double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void   daxpy_(const blas_int* n, const double* alpha, const double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void    dger_(const blas_int* m, const blas_int* n, const double* alpha, const double* x, const blas_int* incx, const double* y, const blas_int* incy, double* A, const blas_int* ldA);
    }
  
  
  
  template<typename eT>
  arma_inline
  eT
  dot_(const blas_int* n, const eT* x, const eT* y)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    const blas_int inc = 1;
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      return eT( sdot_(n, (const T*)x, &inc, (const T*)y, &inc) );
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      return eT( ddot_(n, (const T*)x, &inc, (const T*)y, &inc) );
      }
    else
      {
      return eT(0);  // prevent compiler warnings
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  gemv_(const char* transA, const blas_int* m, const blas_int* n, const eT* alpha, const eT* A, const blas_int* ldA, const eT* x, const blas_int* incx, const eT* beta, eT* y, const blas_int* incy)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgemv_(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgemv_(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      cgemv_(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      zgemv_(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    
    }
  
  
  
  template<typename eT>
  inline
  void
  gemm_(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const eT* alpha, const eT* A, const blas_int* ldA, const eT* B, const blas_int* ldB, const eT* beta, eT* C, const blas_int* ldC)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgemm_(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgemm_(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      cgemm_(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      zgemm_(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    
    }
  
  }

#endif
