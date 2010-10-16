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



#ifdef ARMA_USE_BLAS


//! \namespace blas namespace for BLAS functions
namespace blas
  {
  extern "C"
    {
    float  arma_fortran(sdot)(const blas_int* n, const float*  x, const blas_int* incx, const float*  y, const blas_int* incy);
    double arma_fortran(ddot)(const blas_int* n, const double* x, const blas_int* incx, const double* y, const blas_int* incy);
    
    void arma_fortran(sgemv)(const char* transA, const blas_int* m, const blas_int* n, const float*  alpha, const float*  A, const blas_int* ldA, const float*  x, const blas_int* incx, const float*  beta, float*  y, const blas_int* incy);
    void arma_fortran(dgemv)(const char* transA, const blas_int* m, const blas_int* n, const double* alpha, const double* A, const blas_int* ldA, const double* x, const blas_int* incx, const double* beta, double* y, const blas_int* incy);
    void arma_fortran(cgemv)(const char* transA, const blas_int* m, const blas_int* n, const void*   alpha, const void*   A, const blas_int* ldA, const void*   x, const blas_int* incx, const void*   beta, void*   y, const blas_int* incy);
    void arma_fortran(zgemv)(const char* transA, const blas_int* m, const blas_int* n, const void*   alpha, const void*   A, const blas_int* ldA, const void*   x, const blas_int* incx, const void*   beta, void*   y, const blas_int* incy);
    
    void arma_fortran(sgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const float*  alpha, const float*  A, const blas_int* ldA, const float*  B, const blas_int* ldB, const float*  beta, float*  C, const blas_int* ldC);
    void arma_fortran(dgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const double* alpha, const double* A, const blas_int* ldA, const double* B, const blas_int* ldB, const double* beta, double* C, const blas_int* ldC);
    void arma_fortran(cgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const void*   alpha, const void*   A, const blas_int* ldA, const void*   B, const blas_int* ldB, const void*   beta, void*   C, const blas_int* ldC);
    void arma_fortran(zgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const void*   alpha, const void*   A, const blas_int* ldA, const void*   B, const blas_int* ldB, const void*   beta, void*   C, const blas_int* ldC);
    
    // void   arma_fortran(dswap)(const blas_int* n, double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void   arma_fortran(dscal)(const blas_int* n, const double* alpha, double* x, const blas_int* incx);
    // void   arma_fortran(dcopy)(const blas_int* n, const double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void   arma_fortran(daxpy)(const blas_int* n, const double* alpha, const double* x, const blas_int* incx, double* y, const blas_int* incy);
    // void    arma_fortran(dger)(const blas_int* m, const blas_int* n, const double* alpha, const double* x, const blas_int* incx, const double* y, const blas_int* incy, double* A, const blas_int* ldA);
    }
  
  
  
  template<typename eT>
  arma_inline
  eT
  dot(const blas_int* n, const eT* x, const eT* y)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    const blas_int inc = 1;
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      return eT( arma_fortran(sdot)(n, (const T*)x, &inc, (const T*)y, &inc) );
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      return eT( arma_fortran(ddot)(n, (const T*)x, &inc, (const T*)y, &inc) );
      }
    else
      {
      return eT(0);  // prevent compiler warnings
      }
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
      arma_fortran(sgemv)(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgemv)(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(cgemv)(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(zgemv)(transA, m, n, (const T*)alpha, (const T*)A, ldA, (const T*)x, incx, (const T*)beta, (T*)y, incy);
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
      arma_fortran(sgemm)(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgemm)(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(cgemm)(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(zgemm)(transA, transB, m, n, k, (const T*)alpha, (const T*)A, ldA, (const T*)B, ldB, (const T*)beta, (T*)C, ldC);
      }
    
    }
  
  }


#endif
