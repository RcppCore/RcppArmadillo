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


//! \namespace blas namespace for BLAS functions
namespace blas
  {
  
  
  template<typename eT>
  inline
  eT
  dot(const uword n_elem, const eT* x, const eT* y)
    {
    arma_ignore(n_elem);
    arma_ignore(x);
    arma_ignore(y);
    
    return eT(0);
    }
  
  
  
  template<>
  inline
  float
  dot(const uword n_elem, const float* x, const float* y)
    {
    blas_int n   = blas_int(n_elem);
    blas_int inc = blas_int(1);
    
    return arma_fortran(arma_sdot)(&n, x, &inc, y, &inc);
    }
  
  
  
  template<>
  inline
  double
  dot(const uword n_elem, const double* x, const double* y)
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
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
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
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
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
