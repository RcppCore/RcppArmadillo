// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifdef ARMA_USE_BLAS


//! \namespace blas namespace for BLAS functions
namespace blas
  {
  
  
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
  
  
  
  template<typename eT>
  inline
  eT
  dot(const uword n_elem, const eT* x, const eT* y)
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value == true)
      {
      #if defined(ARMA_BLAS_SDOT_BUG)
        {
        if(n_elem == 0)  { return eT(0); }
        
        const char trans   = 'T';
        
        const blas_int m   = blas_int(n_elem);
        const blas_int n   = 1;
        //const blas_int lda = (n_elem > 0) ? blas_int(n_elem) : blas_int(1);
        const blas_int inc = 1;
        
        const eT alpha     = eT(1);
        const eT beta      = eT(0);
        
        eT result[2];  // paranoia: using two elements instead of one
        
        //blas::gemv(&trans, &m, &n, &alpha, x, &lda, y, &inc, &beta, &result[0], &inc);
        blas::gemv(&trans, &m, &n, &alpha, x, &m, y, &inc, &beta, &result[0], &inc);
        
        return result[0];
        }
      #else
        {
        blas_int n   = blas_int(n_elem);
        blas_int inc = 1;
        
        typedef float T;
        return arma_fortran(arma_sdot)(&n, (const T*)x, &inc, (const T*)y, &inc);
        }
      #endif
      }
    else
    if(is_double<eT>::value == true)
      {
      blas_int n   = blas_int(n_elem);
      blas_int inc = 1;
      
      typedef double T;
      return arma_fortran(arma_ddot)(&n, (const T*)x, &inc, (const T*)y, &inc);
      }
    else
    if( (is_supported_complex_float<eT>::value == true) || (is_supported_complex_double<eT>::value == true) )
      {
      if(n_elem == 0)  { return eT(0); }
      
      // using gemv() workaround due to compatibility issues with cdotu() and zdotu()
      
      const char trans   = 'T';
      
      const blas_int m   = blas_int(n_elem);
      const blas_int n   = 1;
      //const blas_int lda = (n_elem > 0) ? blas_int(n_elem) : blas_int(1);
      const blas_int inc = 1;
      
      const eT alpha     = eT(1);
      const eT beta      = eT(0);
      
      eT result[2];  // paranoia: using two elements instead of one
      
      //blas::gemv(&trans, &m, &n, &alpha, x, &lda, y, &inc, &beta, &result[0], &inc);
      blas::gemv(&trans, &m, &n, &alpha, x, &m, y, &inc, &beta, &result[0], &inc);
      
      return result[0];
      }
    else
      {
      return eT(0);
      }
    }
  }


#endif
