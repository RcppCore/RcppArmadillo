// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// Copyright (C) 2009      Edmund Highcock
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
    void arma_fortran(sgetrf)(blas_int* m, blas_int* n,  float* a, blas_int* lda, blas_int* ipiv, blas_int* info);
    void arma_fortran(dgetrf)(blas_int* m, blas_int* n, double* a, blas_int* lda, blas_int* ipiv, blas_int* info);
    void arma_fortran(cgetrf)(blas_int* m, blas_int* n,   void* a, blas_int* lda, blas_int* ipiv, blas_int* info);
    void arma_fortran(zgetrf)(blas_int* m, blas_int* n,   void* a, blas_int* lda, blas_int* ipiv, blas_int* info);
    
    // matrix inversion
    void arma_fortran(sgetri)(blas_int* n,  float* a, blas_int* lda, blas_int* ipiv,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(dgetri)(blas_int* n, double* a, blas_int* lda, blas_int* ipiv, double* work, blas_int* lwork, blas_int* info);
    void arma_fortran(cgetri)(blas_int* n,  void*  a, blas_int* lda, blas_int* ipiv,   void* work, blas_int* lwork, blas_int* info);
    void arma_fortran(zgetri)(blas_int* n,  void*  a, blas_int* lda, blas_int* ipiv,   void* work, blas_int* lwork, blas_int* info);
    
    // eigenvector decomposition of symmetric real matrices
    void arma_fortran(ssyev)(char* jobz, char* uplo, blas_int* n,  float* a, blas_int* lda,  float* w,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(dsyev)(char* jobz, char* uplo, blas_int* n, double* a, blas_int* lda, double* w, double* work, blas_int* lwork, blas_int* info);

    // eigenvector decomposition of hermitian matrices (complex)
    void arma_fortran(cheev)(char* jobz, char* uplo, blas_int* n,   void* a, blas_int* lda,  float* w,   void* work, blas_int* lwork,  float* rwork, blas_int* info);
    void arma_fortran(zheev)(char* jobz, char* uplo, blas_int* n,   void* a, blas_int* lda, double* w,   void* work, blas_int* lwork, double* rwork, blas_int* info);

    // eigenvector decomposition of general real matrices
    void arma_fortran(sgeev)(char* jobvl, char* jobvr, blas_int* n,  float* a, blas_int* lda,  float* wr,  float* wi,  float* vl, blas_int* ldvl,  float* vr, blas_int* ldvr,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(dgeev)(char* jobvl, char* jobvr, blas_int* n, double* a, blas_int* lda, double* wr, double* wi, double* vl, blas_int* ldvl, double* vr, blas_int* ldvr, double* work, blas_int* lwork, blas_int* info);

    // eigenvector decomposition of general complex matrices
    void arma_fortran(cgeev)(char* jobvr, char* jobvl, blas_int* n, void* a, blas_int* lda, void* w, void* vl, blas_int* ldvl, void* vr, blas_int* ldvr, void* work, blas_int* lwork,  float* rwork, blas_int* info);
    void arma_fortran(zgeev)(char* jobvl, char* jobvr, blas_int* n, void* a, blas_int* lda, void* w, void* vl, blas_int* ldvl, void* vr, blas_int* ldvr, void* work, blas_int* lwork, double* rwork, blas_int* info);

    // Cholesky decomposition
    void arma_fortran(spotrf)(char* uplo, blas_int* n,  float* a, blas_int* lda, blas_int* info);
    void arma_fortran(dpotrf)(char* uplo, blas_int* n, double* a, blas_int* lda, blas_int* info);
    void arma_fortran(cpotrf)(char* uplo, blas_int* n,   void* a, blas_int* lda, blas_int* info);
    void arma_fortran(zpotrf)(char* uplo, blas_int* n,   void* a, blas_int* lda, blas_int* info);

    // QR decomposition
    void arma_fortran(sgeqrf)(blas_int* m, blas_int* n,  float* a, blas_int* lda,  float* tau,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(dgeqrf)(blas_int* m, blas_int* n, double* a, blas_int* lda, double* tau, double* work, blas_int* lwork, blas_int* info);
    void arma_fortran(cgeqrf)(blas_int* m, blas_int* n,   void* a, blas_int* lda,   void* tau,   void* work, blas_int* lwork, blas_int* info);
    void arma_fortran(zgeqrf)(blas_int* m, blas_int* n,   void* a, blas_int* lda,   void* tau,   void* work, blas_int* lwork, blas_int* info);

    // Q matrix calculation from QR decomposition (real matrices)
    void arma_fortran(sorgqr)(blas_int* m, blas_int* n, blas_int* k,  float* a, blas_int* lda,  float* tau,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(dorgqr)(blas_int* m, blas_int* n, blas_int* k, double* a, blas_int* lda, double* tau, double* work, blas_int* lwork, blas_int* info);
    
    // Q matrix calculation from QR decomposition (complex matrices)
    void arma_fortran(cungqr)(blas_int* m, blas_int* n, blas_int* k,   void* a, blas_int* lda,   void* tau,   void* work, blas_int* lwork, blas_int* info);
    void arma_fortran(zungqr)(blas_int* m, blas_int* n, blas_int* k,   void* a, blas_int* lda,   void* tau,   void* work, blas_int* lwork, blas_int* info);

    // SVD (real matrices)
    void arma_fortran(sgesvd)(char* jobu, char* jobvt, blas_int* m, blas_int* n, float*  a, blas_int* lda, float*  s, float*  u, blas_int* ldu, float*  vt, blas_int* ldvt, float*  work, blas_int* lwork, blas_int* info);
    void arma_fortran(dgesvd)(char* jobu, char* jobvt, blas_int* m, blas_int* n, double* a, blas_int* lda, double* s, double* u, blas_int* ldu, double* vt, blas_int* ldvt, double* work, blas_int* lwork, blas_int* info);
    
    // SVD (complex matrices)
    void arma_fortran(cgesvd)(char* jobu, char* jobvt, blas_int* m, blas_int* n, void*   a, blas_int* lda, float*  s, void*   u, blas_int* ldu, void*   vt, blas_int* ldvt, void*   work, blas_int* lwork, float*  rwork, blas_int* info);
    void arma_fortran(zgesvd)(char* jobu, char* jobvt, blas_int* m, blas_int* n, void*   a, blas_int* lda, double* s, void*   u, blas_int* ldu, void*   vt, blas_int* ldvt, void*   work, blas_int* lwork, double* rwork, blas_int* info);

    // solve system of linear equations, using LU decomposition
    void arma_fortran(sgesv)(blas_int* n, blas_int* nrhs, float*  a, blas_int* lda, blas_int* ipiv, float*  b, blas_int* ldb, blas_int* info);
    void arma_fortran(dgesv)(blas_int* n, blas_int* nrhs, double* a, blas_int* lda, blas_int* ipiv, double* b, blas_int* ldb, blas_int* info);
    void arma_fortran(cgesv)(blas_int* n, blas_int* nrhs, void*   a, blas_int* lda, blas_int* ipiv, void*   b, blas_int* ldb, blas_int* info);
    void arma_fortran(zgesv)(blas_int* n, blas_int* nrhs, void*   a, blas_int* lda, blas_int* ipiv, void*   b, blas_int* ldb, blas_int* info);

    // solve over/underdetermined system of linear equations
    void arma_fortran(sgels)(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, float*  a, blas_int* lda, float*  b, blas_int* ldb, float*  work, blas_int* lwork, blas_int* info);
    void arma_fortran(dgels)(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, double* a, blas_int* lda, double* b, blas_int* ldb, double* work, blas_int* lwork, blas_int* info);
    void arma_fortran(cgels)(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, void*   a, blas_int* lda, void*   b, blas_int* ldb, void*   work, blas_int* lwork, blas_int* info);
    void arma_fortran(zgels)(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, void*   a, blas_int* lda, void*   b, blas_int* ldb, void*   work, blas_int* lwork, blas_int* info);

    // void arma_fortran(dgeqp3)(blas_int* m, blas_int* n, double* a, blas_int* lda, blas_int* jpvt, double* tau, double* work, blas_int* lwork, blas_int* info);
    // void arma_fortran(dormqr)(char* side, char* trans, blas_int* m, blas_int* n, blas_int* k, double* a, blas_int* lda, double* tau, double* c, blas_int* ldc, double* work, blas_int* lwork, blas_int* info);
    // void  arma_fortran(dposv)(char* uplo, blas_int* n, blas_int* nrhs, double* a, blas_int* lda, double* b, blas_int* ldb, blas_int* info);
    // void arma_fortran(dtrtrs)(char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, double* a, blas_int* lda, double* b, blas_int* ldb, blas_int* info);
    // void  arma_fortran(dgees)(char* jobvs, char* sort, blas_int* select, blas_int* n, double* a, blas_int* lda, blas_int* sdim, double* wr, double* wi, double* vs, blas_int* ldvs, double* work, blas_int* lwork, blas_int* bwork, blas_int* info);
    }
  
  
  
  template<typename eT>
  inline
  void
  getrf(blas_int* m, blas_int* n, eT* a, blas_int* lda, blas_int* ipiv, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sgetrf)(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgetrf)(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(cgetrf)(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(zgetrf)(m, n, (T*)a, lda, ipiv, info);
      }
    }
    
    
    
  template<typename eT>
  inline
  void
  getri(blas_int* n,  eT* a, blas_int* lda, blas_int* ipiv, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(cgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(zgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  syev(char* jobz, char* uplo, blas_int* n, eT* a, blas_int* lda, eT* w,  eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(ssyev)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dsyev)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  heev
    (
    char* jobz, char* uplo, blas_int* n,
    eT* a, blas_int* lda, typename eT::value_type* w,
    eT* work, blas_int* lwork, typename eT::value_type* rwork,
    blas_int* info
    )
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef float T;
      typedef typename std::complex<T> cx_T;
      arma_fortran(cheev)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      typedef typename std::complex<T> cx_T;
      arma_fortran(zheev)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info);
      }
    }
  
	 
  template<typename eT>
  inline
  void
  geev
    (
    char* jobvl, char* jobvr, blas_int* n, 
    eT* a, blas_int* lda, eT* wr, eT* wi, eT* vl, 
    blas_int* ldvl, eT* vr, blas_int* ldvr, 
    eT* work, blas_int* lwork,
    blas_int* info
    )
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();

    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sgeev)(jobvl, jobvr, n,  (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgeev)(jobvl, jobvr, n,  (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info);
      }
    }


  template<typename eT>
  inline
  void
  cx_geev
    (
    char* jobvl, char* jobvr, blas_int* n, 
    eT* a, blas_int* lda, eT* w, 
    eT* vl, blas_int* ldvl, 
    eT* vr, blas_int* ldvr, 
    eT* work, blas_int* lwork, typename eT::value_type* rwork, 
    blas_int* info
    )
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef float T;
      typedef typename std::complex<T> cx_T;
      arma_fortran(cgeev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      typedef typename std::complex<T> cx_T;
      arma_fortran(zgeev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info);
      }
    }
  


  
  template<typename eT>
  inline
  void
  potrf(char* uplo, blas_int* n, eT* a, blas_int* lda, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(spotrf)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dpotrf)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(cpotrf)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(zpotrf)(uplo, n, (T*)a, lda, info);
      }
    
    }
  
  
  
  template<typename eT>
  inline
  void
  geqrf(blas_int* m, blas_int* n, eT* a, blas_int* lda, eT* tau, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(cgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(zgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    
    }
  
  
  
  template<typename eT>
  inline
  void
  orgqr(blas_int* m, blas_int* n, blas_int* k, eT* a, blas_int* lda, eT* tau, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sorgqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dorgqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    }


  
  template<typename eT>
  inline
  void
  ungqr(blas_int* m, blas_int* n, blas_int* k, eT* a, blas_int* lda, eT* tau, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(cungqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(zungqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    }
  
  
  template<typename eT>
  inline
  void
  gesvd
    (
    char* jobu, char* jobvt, blas_int* m, blas_int* n, eT* a, blas_int* lda,
    eT* s, eT* u, blas_int* ldu, eT* vt, blas_int* ldvt,
    eT* work, blas_int* lwork, blas_int* info
    )
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sgesvd)(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgesvd)(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info);
      }
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gesvd
    (
    char* jobu, char* jobvt, blas_int* m, blas_int* n, std::complex<T>* a, blas_int* lda,
    T* s, std::complex<T>* u, blas_int* ldu, std::complex<T>* vt, blas_int* ldvt, 
    std::complex<T>* work, blas_int* lwork, T* rwork, blas_int* info
    )
    {
    arma_type_check<is_supported_blas_type<T>::value == false>::apply();
    arma_type_check<is_supported_blas_type< std::complex<T> >::value == false>::apply();
    
    if(is_float<T>::value == true)
      {
      typedef float bT;
      arma_fortran(cgesvd)
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
      arma_fortran(zgesvd)
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
  gesv(blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, blas_int* ipiv, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(cgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(zgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  gels(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, eT* b, blas_int* ldb, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(cgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(zgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    }
  
  
  
  //! @}
  }


#endif
