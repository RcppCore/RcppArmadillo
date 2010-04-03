// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Edmund Highcock (edmund dot highcock at merton dot ox dot ac dot uk)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


#ifdef ARMA_USE_LAPACK

//! \namespace lapack namespace for LAPACK functions
namespace lapack
  {
  //! \addtogroup LAPACK
  //! @{
  
  extern "C"
    {
    // LU factorisation
    void sgetrf_(int* m, int* n,  float* a, int* lda, int* ipiv, int* info);
    void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
    void cgetrf_(int* m, int* n,   void* a, int* lda, int* ipiv, int* info);
    void zgetrf_(int* m, int* n,   void* a, int* lda, int* ipiv, int* info);
    
    // matrix inversion
    void sgetri_(int* n,  float* a, int* lda, int* ipiv,  float* work, int* lwork, int* info);
    void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);
    void cgetri_(int* n,  void*  a, int* lda, int* ipiv,   void* work, int* lwork, int* info);
    void zgetri_(int* n,  void*  a, int* lda, int* ipiv,   void* work, int* lwork, int* info);
    
    // eigenvector decomposition of symmetric real matrices
    void ssyev_(char* jobz, char* uplo, int* n,  float* a, int* lda,  float* w,  float* work, int* lwork, int* info);
    void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);

    // eigenvector decomposition of hermitian matrices (complex)
    void cheev_(char* jobz, char* uplo, int* n,   void* a, int* lda,  float* w,   void* work, int* lwork,  float* rwork, int* info);
    void zheev_(char* jobz, char* uplo, int* n,   void* a, int* lda, double* w,   void* work, int* lwork, double* rwork, int* info);

    // eigenvector decomposition of general real matrices
    void sgeev_(char* jobvl, char* jobvr, int* n,  float* a, int* lda,  float* wr,  float* wi,  float* vl, int* ldvl,  float* vr, int* ldvr,  float* work, int* lwork, int* info);
    void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info);

    // eigenvector decomposition of general complex matrices
    void cgeev_(char* jobvr, char* jobvl, int* n, void* a, int* lda, void* w, void* vl, int* ldvl, void* vr, int* ldvr, void* work, int* lwork,  float* rwork, int* info);
    void zgeev_(char* jobvl, char* jobvr, int* n, void* a, int* lda, void* w, void* vl, int *ldvl, void* vr, int *ldvr, void* work, int* lwork, double* rwork, int* info);

    // Cholesky decomposition
    void spotrf_(char* uplo, int* n,  float* a, int* lda, int* info);
    void dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
    void cpotrf_(char* uplo, int* n,   void* a, int* lda, int* info);
    void zpotrf_(char* uplo, int* n,   void* a, int* lda, int* info);

    // QR decomposition
    void sgeqrf_(int* m, int* n,  float* a, int* lda,  float* tau,  float* work, int* lwork, int* info);
    void dgeqrf_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info);
    void cgeqrf_(int* m, int* n,   void* a, int* lda,   void* tau,   void* work, int* lwork, int* info);
    void zgeqrf_(int* m, int* n,   void* a, int* lda,   void* tau,   void* work, int* lwork, int* info);

    // Q matrix calculation from QR decomposition (real matrices)
    void sorgqr_(int* m, int* n, int* k,  float* a, int* lda,  float* tau,  float* work, int* lwork, int* info);
    void dorgqr_(int* m, int* n, int* k, double* a, int* lda, double* tau, double* work, int* lwork, int* info);
    
    // Q matrix calculation from QR decomposition (complex matrices)
    void cungqr_(int* m, int* n, int* k,   void* a, int* lda,   void* tau,   void* work, int* lwork, int* info);
    void zungqr_(int* m, int* n, int* k,   void* a, int* lda,   void* tau,   void* work, int* lwork, int* info);

    // SVD (real matrices)
    void sgesvd_(char* jobu, char* jobvt, int* m, int* n, float*  a, int* lda, float*  s, float*  u, int* ldu, float*  vt, int* ldvt, float*  work, int* lwork, int* info);
    void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info);
    
    // SVD (complex matrices)
    void cgesvd_(char* jobu, char* jobvt, int* m, int* n, void*   a, int* lda, float*  s, void*   u, int* ldu, void*   vt, int* ldvt, void*   work, int* lwork, float*  rwork, int* info);
    void zgesvd_(char* jobu, char* jobvt, int* m, int* n, void*   a, int* lda, double* s, void*   u, int* ldu, void*   vt, int* ldvt, void*   work, int* lwork, double* rwork, int* info);

    // solve system of linear equations, using LU decomposition
    void sgesv_(int* n, int* nrhs, float*  a, int* lda, int* ipiv, float*  b, int* ldb, int* info);
    void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
    void cgesv_(int* n, int* nrhs, void*   a, int* lda, int* ipiv, void*   b, int* ldb, int* info);
    void zgesv_(int* n, int* nrhs, void*   a, int* lda, int* ipiv, void*   b, int* ldb, int* info);

    // solve over/underdetermined system of linear equations
    void sgels_(char* trans, int* m, int* n, int* nrhs, float*  a, int* lda, float*  b, int* ldb, float*  work, int* lwork, int* info);
    void dgels_(char* trans, int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* work, int* lwork, int* info);
    void cgels_(char* trans, int* m, int* n, int* nrhs, void*   a, int *lda, void*   b, int* ldb, void*   work, int* lwork, int* info);
    void zgels_(char* trans, int* m, int* n, int* nrhs, void*   a, int *lda, void*   b, int* ldb, void*   work, int* lwork, int* info);

    // void dgeqp3_(int* m, int* n, double* a, int* lda, int* jpvt, double* tau, double* work, int* lwork, int* info);
    // void dormqr_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* lwork, int* info);
    // void  dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* info);
    // void dtrtrs_(char* uplo, char* trans, char* diag, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* info);
    // void  dgees_(char* jobvs, char* sort, int* select, int* n, double* a, int* lda, int* sdim, double* wr, double* wi, double* vs, int* ldvs, double* work, int* lwork, int* bwork, int* info);
    
    }

  template<typename eT>
  inline
  void
  getrf_(int* m, int* n, eT* a, int* lda, int* ipiv, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgetrf_(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgetrf_(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      cgetrf_(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      zgetrf_(m, n, (T*)a, lda, ipiv, info);
      }
    }
    
    
    
  template<typename eT>
  inline
  void
  getri_(int* n,  eT* a, int* lda, int* ipiv, eT* work, int* lwork, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgetri_(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgetri_(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      cgetri_(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      zgetri_(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  syev_(char* jobz, char* uplo, int* n, eT* a, int* lda, eT* w,  eT* work, int* lwork, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      ssyev_(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dsyev_(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  heev_
    (
    char* jobz, char* uplo, int* n,
    eT* a, int* lda, typename eT::value_type* w,
    eT* work, int* lwork, typename eT::value_type* rwork,
    int* info
    )
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef float T;
      typedef typename std::complex<T> cx_T;
      cheev_(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      typedef typename std::complex<T> cx_T;
      zheev_(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info);
      }
    }
  
	 
  template<typename eT>
  inline
  void
  geev_
    (
    char* jobvl, char* jobvr, int* n, 
    eT* a, int* lda, eT* wr, eT* wi, eT* vl, 
    int* ldvl, eT* vr, int* ldvr, 
    eT* work, int* lwork,
    int* info
    )
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();

    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgeev_(jobvl, jobvr, n,  (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgeev_(jobvl, jobvr, n,  (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info);
      }
    }


  template<typename eT>
  inline
  void
  cx_geev_
    (
    char* jobvl, char* jobvr, int* n, 
    eT* a, int* lda, eT* w, 
    eT* vl, int* ldvl, 
    eT* vr, int* ldvr, 
    eT* work, int* lwork, typename eT::value_type* rwork, 
    int* info
    )
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef float T;
      typedef typename std::complex<T> cx_T;
      cgeev_(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      typedef typename std::complex<T> cx_T;
      zgeev_(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info);
      }
    }
  


  
  template<typename eT>
  inline
  void
  potrf_(char* uplo, int* n, eT* a, int* lda, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      spotrf_(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dpotrf_(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      cpotrf_(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      zpotrf_(uplo, n, (T*)a, lda, info);
      }
    
    }
  
  
  
  template<typename eT>
  inline
  void
  geqrf_(int* m, int* n, eT* a, int* lda, eT* tau, eT* work, int* lwork, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgeqrf_(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgeqrf_(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      cgeqrf_(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      zgeqrf_(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    
    }
  
  
  
  template<typename eT>
  inline
  void
  orgqr_(int* m, int* n, int* k, eT* a, int* lda, eT* tau, eT* work, int* lwork, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sorgqr_(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dorgqr_(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    }


  
  template<typename eT>
  inline
  void
  ungqr_(int* m, int* n, int* k, eT* a, int* lda, eT* tau, eT* work, int* lwork, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef float T;
      cungqr_(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      zungqr_(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    }
  
  
  template<typename eT>
  inline
  void
  gesvd_
    (
    char* jobu, char* jobvt, int* m, int* n, eT* a, int* lda,
    eT* s, eT* u, int* ldu, eT* vt, int* ldvt,
    eT* work, int* lwork, int* info
    )
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgesvd_(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgesvd_(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info);
      }
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gesvd_
    (
    char* jobu, char* jobvt, int* m, int* n, std::complex<T>* a, int* lda,
    T* s, std::complex<T>* u, int* ldu, std::complex<T>* vt, int* ldvt, 
    std::complex<T>* work, int* lwork, T* rwork, int* info
    )
    {
    arma_type_check<is_supported_blas_type<T>::value == false>::apply();
    arma_type_check<is_supported_blas_type< std::complex<T> >::value == false>::apply();
    
    if(is_float<T>::value == true)
      {
      typedef float bT;
      cgesvd_
        (
        jobu, jobvt, m, n, (std::complex<bT>*)a, lda,
        (bT*)s, (std::complex<bT>*)u, ldu, (std::complex<bT>*)vt, ldvt,
        (std::complex<bT>*)work, lwork, (bT*)rwork, info
        );
      }
    else
    if(is_double<T>::value == true)
      {
      typedef double bT;
      zgesvd_
        (
        jobu, jobvt, m, n, (std::complex<bT>*)a, lda,
        (bT*)s, (std::complex<bT>*)u, ldu, (std::complex<bT>*)vt, ldvt,
        (std::complex<bT>*)work, lwork, (bT*)rwork, info
        );
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  gesv_(int* n, int* nrhs, eT* a, int* lda, int* ipiv, eT* b, int* ldb, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgesv_(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgesv_(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      cgesv_(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      zgesv_(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
//sgels_(char* trans, int* m, int* n, int* nrhs, float*  a, int* lda, float*  b, int* ldb, float*  work, int* lwork, int* info);
  gels_(char* trans, int* m, int* n, int* nrhs, eT* a, int* lda, eT* b, int* ldb, eT* work, int* lwork, int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      sgels_(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      dgels_(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      cgels_(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      zgels_(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    }
  
  
  
  //! @}
  }

#endif

