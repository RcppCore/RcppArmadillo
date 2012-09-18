// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


#ifdef ARMA_USE_ATLAS


//! \namespace atlas namespace for ATLAS functions (imported from the global namespace)
namespace atlas
  {
  
  template<typename eT>
  inline static const eT& tmp_real(const eT& X)              { return X; }
  
  template<typename T>
  inline static const  T  tmp_real(const std::complex<T>& X) { return X.real(); }
  
  
  
  template<typename eT>
  arma_inline
  eT
  cblas_dot(const int N, const eT* X, const eT* Y)
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      return eT( arma_atlas(cblas_sdot)(N, (const T*)X, 1, (const T*)Y, 1) );
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      return eT( arma_atlas(cblas_ddot)(N, (const T*)X, 1, (const T*)Y, 1) );
      }
    else
      {
      return eT(0);
      }
    }
  
  
  
  template<typename eT>
  arma_inline
  eT
  cx_cblas_dot(const int N, const eT* X, const eT* Y)
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef typename std::complex<float> T;
      
      T out;    
      arma_atlas(cblas_cdotu_sub)(N, (const T*)X, 1, (const T*)Y, 1, &out);
      
      return eT(out);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef typename std::complex<double> T;
      
      T out;
      arma_atlas(cblas_zdotu_sub)(N, (const T*)X, 1, (const T*)Y, 1, &out);
      
      return eT(out);
      }
    else
      {
      return eT(0);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  cblas_gemv
    (
    const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
    const int M, const int N,
    const eT alpha,
    const eT *A, const int lda,
    const eT *X, const int incX,
    const eT beta,
    eT *Y, const int incY
    )
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_atlas(cblas_sgemv)(Order, TransA, M, N, (const T)tmp_real(alpha), (const T*)A, lda, (const T*)X, incX, (const T)tmp_real(beta), (T*)Y, incY);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_atlas(cblas_dgemv)(Order, TransA, M, N, (const T)tmp_real(alpha), (const T*)A, lda, (const T*)X, incX, (const T)tmp_real(beta), (T*)Y, incY);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_atlas(cblas_cgemv)(Order, TransA, M, N, (const T*)&alpha, (const T*)A, lda, (const T*)X, incX, (const T*)&beta, (T*)Y, incY);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_atlas(cblas_zgemv)(Order, TransA, M, N, (const T*)&alpha, (const T*)A, lda, (const T*)X, incX, (const T*)&beta, (T*)Y, incY);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  cblas_gemm
    (
    const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
    const int K, const eT alpha, const eT *A,
    const int lda, const eT *B, const int ldb,
    const eT beta, eT *C, const int ldc
    )
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_atlas(cblas_sgemm)(Order, TransA, TransB, M, N, K, (const T)tmp_real(alpha), (const T*)A, lda, (const T*)B, ldb, (const T)tmp_real(beta), (T*)C, ldc);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_atlas(cblas_dgemm)(Order, TransA, TransB, M, N, K, (const T)tmp_real(alpha), (const T*)A, lda, (const T*)B, ldb, (const T)tmp_real(beta), (T*)C, ldc);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_atlas(cblas_cgemm)(Order, TransA, TransB, M, N, K, (const T*)&alpha, (const T*)A, lda, (const T*)B, ldb, (const T*)&beta, (T*)C, ldc);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_atlas(cblas_zgemm)(Order, TransA, TransB, M, N, K, (const T*)&alpha, (const T*)A, lda, (const T*)B, ldb, (const T*)&beta, (T*)C, ldc);
      }
    }
  
  
  
  template<typename eT>
  inline
  int
  clapack_getrf
    (
    const enum CBLAS_ORDER Order, const int M, const int N,
    eT *A, const int lda, int *ipiv
    )
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      return arma_atlas(clapack_sgetrf)(Order, M, N, (T*)A, lda, ipiv);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      return arma_atlas(clapack_dgetrf)(Order, M, N, (T*)A, lda, ipiv);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      return arma_atlas(clapack_cgetrf)(Order, M, N, (T*)A, lda, ipiv);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      return arma_atlas(clapack_zgetrf)(Order, M, N, (T*)A, lda, ipiv);
      }
    else
      {
      return -1;
      }
    }
  
  
  
  template<typename eT>
  inline
  int
  clapack_getri
    (
    const enum CBLAS_ORDER Order, const int N, eT *A,
    const int lda, const int *ipiv
    )
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      return arma_atlas(clapack_sgetri)(Order, N, (T*)A, lda, ipiv);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      return arma_atlas(clapack_dgetri)(Order, N, (T*)A, lda, ipiv);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      return arma_atlas(clapack_cgetri)(Order, N, (T*)A, lda, ipiv);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      return arma_atlas(clapack_zgetri)(Order, N, (T*)A, lda, ipiv);
      }
    else
      {
      return -1;
      }
    }
  
  
  
  template<typename eT>
  inline
  int
  clapack_gesv
    (
    const enum CBLAS_ORDER Order,
    const int N, const int NRHS,
    eT* A, const int lda, int* ipiv,
    eT* B, const int ldb
    )
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      return arma_atlas(clapack_sgesv)(Order, N, NRHS, (T*)A, lda, ipiv, (T*)B, ldb);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      return arma_atlas(clapack_dgesv)(Order, N, NRHS, (T*)A, lda, ipiv, (T*)B, ldb);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      return arma_atlas(clapack_cgesv)(Order, N, NRHS, (T*)A, lda, ipiv, (T*)B, ldb);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      return arma_atlas(clapack_zgesv)(Order, N, NRHS, (T*)A, lda, ipiv, (T*)B, ldb);
      }
    else
      {
      return -1;
      }
    }



  }

#endif
