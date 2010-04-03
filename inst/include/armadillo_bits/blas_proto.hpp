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
    void sgemv_(const char* transA, const int* m, const int* n, const float*  alpha, const float*  A, const int* ldA, const float*  x, const int* incx, const float*  beta, float*  y, const int* incy);
    void dgemv_(const char* transA, const int* m, const int* n, const double* alpha, const double* A, const int* ldA, const double* x, const int* incx, const double* beta, double* y, const int* incy);
    void cgemv_(const char* transA, const int* m, const int* n, const void*   alpha, const void*   A, const int* ldA, const void*   x, const int* incx, const void*   beta, void*   y, const int* incy);
    void zgemv_(const char* transA, const int* m, const int* n, const void*   alpha, const void*   A, const int* ldA, const void*   x, const int* incx, const void*   beta, void*   y, const int* incy);
    
    void sgemm_(const char* transA, const char* transB, const int* m, const int* n, const int* k, const float*  alpha, const float*  A, const int* ldA, const float*  B, const int* ldB, const float*  beta, float*  C, const int* ldC);
    void dgemm_(const char* transA, const char* transB, const int* m, const int* n, const int* k, const double* alpha, const double* A, const int* ldA, const double* B, const int* ldB, const double* beta, double* C, const int* ldC);
    void cgemm_(const char* transA, const char* transB, const int* m, const int* n, const int* k, const void*   alpha, const void*   A, const int* ldA, const void*   B, const int* ldB, const void*   beta, void*   C, const int* ldC);
    void zgemm_(const char* transA, const char* transB, const int* m, const int* n, const int* k, const void*   alpha, const void*   A, const int* ldA, const void*   B, const int* ldB, const void*   beta, void*   C, const int* ldC);

//     float  sdot_(const int* n, const float*  x, const int* incx, const float*  y, const int* incy);
//     double ddot_(const int* n, const double* x, const int* incx, const double* y, const int* incy);

//     void   dswap_(const int* n, double* x, const int* incx, double* y, const int* incy);
//     void   dscal_(const int* n, const double* alpha, double* x, const int* incx);
//     void   dcopy_(const int* n, const double* x, const int* incx, double* y, const int* incy);
//     void   daxpy_(const int* n, const double* alpha, const double* x, const int* incx, double* y, const int* incy);
//     void    dger_(const int* m, const int* n, const double* alpha, const double* x, const int* incx, const double* y, const int* incy, double* A, const int* ldA);
    }
  
  
  
//   template<typename eT>
//   inline
//   eT
//   dot_(const int* n, const eT* x, const int* incx, const eT* y, const int* incy)
//     {
//     arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
//     
//     if(is_float<eT>::value == true)
//       {
//       typedef float T;
//       return sdot_(n, (const T*)x, incx, (const T*)y, incy);
//       }
//     else
//     if(is_double<eT>::value == true)
//       {
//       typedef double T;
//       return ddot_(n, (const T*)x, incx, (const T*)y, incy);
//       }
//     
//     return eT();  // prevent compiler warnings
//     }
  
  
  
  template<typename eT>
  inline
  void
  gemv_(const char* transA, const int* m, const int* n, const eT* alpha, const eT* A, const int* ldA, const eT* x, const int* incx, const eT* beta, eT* y, const int* incy)
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
  gemm_(const char* transA, const char* transB, const int* m, const int* n, const int* k, const eT* alpha, const eT* A, const int* ldA, const eT* B, const int* ldB, const eT* beta, eT* C, const int* ldC)
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
