// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// Copyright (C)      2009 Edmund Highcock
// Copyright (C)      2011 James Sanders
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


#if !defined(ARMA_BLAS_CAPITALS)
  
  #define arma_sgetrf sgetrf
  #define arma_dgetrf dgetrf
  #define arma_cgetrf cgetrf
  #define arma_zgetrf zgetrf
  
  #define arma_sgetri sgetri
  #define arma_dgetri dgetri
  #define arma_cgetri cgetri
  #define arma_zgetri zgetri
  
  #define arma_strtri strtri
  #define arma_dtrtri dtrtri
  #define arma_ctrtri ctrtri
  #define arma_ztrtri ztrtri
  
  #define arma_ssyev  ssyev
  #define arma_dsyev  dsyev
  
  #define arma_cheev  cheev
  #define arma_zheev  zheev
  
  #define arma_sgeev  sgeev
  #define arma_dgeev  dgeev
  
  #define arma_cgeev  cgeev
  #define arma_zgeev  zgeev
  
  #define arma_spotrf spotrf
  #define arma_dpotrf dpotrf
  #define arma_cpotrf cpotrf
  #define arma_zpotrf zpotrf
  
  #define arma_spotri spotri
  #define arma_dpotri dpotri
  #define arma_cpotri cpotri
  #define arma_zpotri zpotri
  
  #define arma_sgeqrf sgeqrf
  #define arma_dgeqrf dgeqrf
  #define arma_cgeqrf cgeqrf
  #define arma_zgeqrf zgeqrf
  
  #define arma_sorgqr sorgqr
  #define arma_dorgqr dorgqr
  
  #define arma_cungqr cungqr
  #define arma_zungqr zungqr
  
  #define arma_sgesvd sgesvd
  #define arma_dgesvd dgesvd
  
  #define arma_cgesvd cgesvd
  #define arma_zgesvd zgesvd
  
  #define arma_sgesv  sgesv
  #define arma_dgesv  dgesv
  #define arma_cgesv  cgesv
  #define arma_zgesv  zgesv
  
  #define arma_sgels  sgels
  #define arma_dgels  dgels
  #define arma_cgels  cgels
  #define arma_zgels  zgels
  
  #define arma_strtrs strtrs
  #define arma_dtrtrs dtrtrs
  #define arma_ctrtrs ctrtrs
  #define arma_ztrtrs ztrtrs

  #define arma_sgees  sgees
  #define arma_dgees  dgees
  #define arma_cgees  cgees
  #define arma_zgees  zgees
  
  #define arma_strsyl strsyl
  #define arma_dtrsyl dtrsyl
  #define arma_ctrsyl ctrsyl
  #define arma_ztrsyl ztrsyl
  
  #define arma_ssytrf ssytrf
  #define arma_dsytrf dsytrf
  #define arma_csytrf csytrf
  #define arma_zsytrf zsytrf
  
  #define arma_ssytri ssytri
  #define arma_dsytri dsytri
  #define arma_csytri csytri
  #define arma_zsytri zsytri
  
#else
  
  #define arma_sgetrf SGETRF
  #define arma_dgetrf DGETRF
  #define arma_cgetrf CGETRF
  #define arma_zgetrf ZGETRF
  
  #define arma_sgetri SGETRI
  #define arma_dgetri DGETRI
  #define arma_cgetri CGETRI
  #define arma_zgetri ZGETRI
  
  #define arma_strtri STRTRI
  #define arma_dtrtri DTRTRI
  #define arma_ctrtri CTRTRI
  #define arma_ztrtri ZTRTRI
  
  #define arma_ssyev  SSYEV
  #define arma_dsyev  DSYEV
  
  #define arma_cheev  CHEEV
  #define arma_zheev  ZHEEV
  
  #define arma_sgeev  SGEEV
  #define arma_dgeev  DGEEV
  
  #define arma_cgeev  CGEEV
  #define arma_zgeev  ZGEEV
  
  #define arma_spotrf SPOTRF
  #define arma_dpotrf DPOTRF
  #define arma_cpotrf CPOTRF
  #define arma_zpotrf ZPOTRF
  
  #define arma_spotri SPOTRI
  #define arma_dpotri DPOTRI
  #define arma_cpotri CPOTRI
  #define arma_zpotri ZPOTRI
  
  #define arma_sgeqrf SGEQRF
  #define arma_dgeqrf DGEQRF
  #define arma_cgeqrf CGEQRF
  #define arma_zgeqrf ZGEQRF
  
  #define arma_sorgqr SORGQR
  #define arma_dorgqr DORGQR
  
  #define arma_cungqr CUNGQR
  #define arma_zungqr ZUNGQR
  
  #define arma_sgesvd SGESVD
  #define arma_dgesvd DGESVD
  
  #define arma_cgesvd CGESVD
  #define arma_zgesvd ZGESVD
  
  #define arma_sgesv  SGESV
  #define arma_dgesv  DGESV
  #define arma_cgesv  CGESV
  #define arma_zgesv  ZGESV
  
  #define arma_sgels  SGELS
  #define arma_dgels  DGELS
  #define arma_cgels  CGELS
  #define arma_zgels  ZGELS
  
  #define arma_strtrs STRTRS
  #define arma_dtrtrs DTRTRS
  #define arma_ctrtrs CTRTRS
  #define arma_ztrtrs ZTRTRS

  #define arma_sgees  SGEES
  #define arma_dgees  DGEES
  #define arma_cgees  CGEES
  #define arma_zgees  ZGEES

  #define arma_strsyl STRSYL
  #define arma_dtrsyl DTRSYL
  #define arma_ctrsyl CTRSYL
  #define arma_ztrsyl ZTRSYL
  
  #define arma_ssytrf SSYTRF
  #define arma_dsytrf DSYTRF
  #define arma_csytrf CSYTRF
  #define arma_zsytrf ZSYTRF
  
  #define arma_ssytri SSYTRI
  #define arma_dsytri DSYTRI
  #define arma_csytri CSYTRI
  #define arma_zsytri ZSYTRI
  
#endif



//! \namespace lapack namespace for LAPACK functions
namespace lapack
  {
  //! \addtogroup LAPACK
  //! @{
  
  extern "C"
    {
    // LU factorisation
    void arma_fortran(arma_sgetrf)(blas_int* m, blas_int* n,  float* a, blas_int* lda, blas_int* ipiv, blas_int* info);
    void arma_fortran(arma_dgetrf)(blas_int* m, blas_int* n, double* a, blas_int* lda, blas_int* ipiv, blas_int* info);
    void arma_fortran(arma_cgetrf)(blas_int* m, blas_int* n,   void* a, blas_int* lda, blas_int* ipiv, blas_int* info);
    void arma_fortran(arma_zgetrf)(blas_int* m, blas_int* n,   void* a, blas_int* lda, blas_int* ipiv, blas_int* info);
    
    // matrix inversion (using LU factorisation result)
    void arma_fortran(arma_sgetri)(blas_int* n,  float* a, blas_int* lda, blas_int* ipiv,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_dgetri)(blas_int* n, double* a, blas_int* lda, blas_int* ipiv, double* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_cgetri)(blas_int* n,  void*  a, blas_int* lda, blas_int* ipiv,   void* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_zgetri)(blas_int* n,  void*  a, blas_int* lda, blas_int* ipiv,   void* work, blas_int* lwork, blas_int* info);
    
    // matrix inversion (triangular matrices)
    void arma_fortran(arma_strtri)(char* uplo, char* diag, blas_int* n,  float* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_dtrtri)(char* uplo, char* diag, blas_int* n, double* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_ctrtri)(char* uplo, char* diag, blas_int* n,   void* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_ztrtri)(char* uplo, char* diag, blas_int* n,   void* a, blas_int* lda, blas_int* info);
    
    // eigenvector decomposition of symmetric real matrices
    void arma_fortran(arma_ssyev)(char* jobz, char* uplo, blas_int* n,  float* a, blas_int* lda,  float* w,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_dsyev)(char* jobz, char* uplo, blas_int* n, double* a, blas_int* lda, double* w, double* work, blas_int* lwork, blas_int* info);
    
    // eigenvector decomposition of hermitian matrices (complex)
    void arma_fortran(arma_cheev)(char* jobz, char* uplo, blas_int* n,   void* a, blas_int* lda,  float* w,   void* work, blas_int* lwork,  float* rwork, blas_int* info);
    void arma_fortran(arma_zheev)(char* jobz, char* uplo, blas_int* n,   void* a, blas_int* lda, double* w,   void* work, blas_int* lwork, double* rwork, blas_int* info);
    
    // eigenvector decomposition of general real matrices
    void arma_fortran(arma_sgeev)(char* jobvl, char* jobvr, blas_int* n,  float* a, blas_int* lda,  float* wr,  float* wi,  float* vl, blas_int* ldvl,  float* vr, blas_int* ldvr,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_dgeev)(char* jobvl, char* jobvr, blas_int* n, double* a, blas_int* lda, double* wr, double* wi, double* vl, blas_int* ldvl, double* vr, blas_int* ldvr, double* work, blas_int* lwork, blas_int* info);
    
    // eigenvector decomposition of general complex matrices
    void arma_fortran(arma_cgeev)(char* jobvr, char* jobvl, blas_int* n, void* a, blas_int* lda, void* w, void* vl, blas_int* ldvl, void* vr, blas_int* ldvr, void* work, blas_int* lwork,  float* rwork, blas_int* info);
    void arma_fortran(arma_zgeev)(char* jobvl, char* jobvr, blas_int* n, void* a, blas_int* lda, void* w, void* vl, blas_int* ldvl, void* vr, blas_int* ldvr, void* work, blas_int* lwork, double* rwork, blas_int* info);
    
    // Cholesky decomposition
    void arma_fortran(arma_spotrf)(char* uplo, blas_int* n,  float* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_dpotrf)(char* uplo, blas_int* n, double* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_cpotrf)(char* uplo, blas_int* n,   void* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_zpotrf)(char* uplo, blas_int* n,   void* a, blas_int* lda, blas_int* info);
    
    // matrix inversion (using Cholesky decomposition result)
    void arma_fortran(arma_spotri)(char* uplo, blas_int* n,  float* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_dpotri)(char* uplo, blas_int* n, double* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_cpotri)(char* uplo, blas_int* n,   void* a, blas_int* lda, blas_int* info);
    void arma_fortran(arma_zpotri)(char* uplo, blas_int* n,   void* a, blas_int* lda, blas_int* info);
    
    // QR decomposition
    void arma_fortran(arma_sgeqrf)(blas_int* m, blas_int* n,  float* a, blas_int* lda,  float* tau,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_dgeqrf)(blas_int* m, blas_int* n, double* a, blas_int* lda, double* tau, double* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_cgeqrf)(blas_int* m, blas_int* n,   void* a, blas_int* lda,   void* tau,   void* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_zgeqrf)(blas_int* m, blas_int* n,   void* a, blas_int* lda,   void* tau,   void* work, blas_int* lwork, blas_int* info);
    
    // Q matrix calculation from QR decomposition (real matrices)
    void arma_fortran(arma_sorgqr)(blas_int* m, blas_int* n, blas_int* k,  float* a, blas_int* lda,  float* tau,  float* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_dorgqr)(blas_int* m, blas_int* n, blas_int* k, double* a, blas_int* lda, double* tau, double* work, blas_int* lwork, blas_int* info);
    
    // Q matrix calculation from QR decomposition (complex matrices)
    void arma_fortran(arma_cungqr)(blas_int* m, blas_int* n, blas_int* k,   void* a, blas_int* lda,   void* tau,   void* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_zungqr)(blas_int* m, blas_int* n, blas_int* k,   void* a, blas_int* lda,   void* tau,   void* work, blas_int* lwork, blas_int* info);
    
    // SVD (real matrices)
    void arma_fortran(arma_sgesvd)(char* jobu, char* jobvt, blas_int* m, blas_int* n, float*  a, blas_int* lda, float*  s, float*  u, blas_int* ldu, float*  vt, blas_int* ldvt, float*  work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_dgesvd)(char* jobu, char* jobvt, blas_int* m, blas_int* n, double* a, blas_int* lda, double* s, double* u, blas_int* ldu, double* vt, blas_int* ldvt, double* work, blas_int* lwork, blas_int* info);
    
    // SVD (complex matrices)
    void arma_fortran(arma_cgesvd)(char* jobu, char* jobvt, blas_int* m, blas_int* n, void*   a, blas_int* lda, float*  s, void*   u, blas_int* ldu, void*   vt, blas_int* ldvt, void*   work, blas_int* lwork, float*  rwork, blas_int* info);
    void arma_fortran(arma_zgesvd)(char* jobu, char* jobvt, blas_int* m, blas_int* n, void*   a, blas_int* lda, double* s, void*   u, blas_int* ldu, void*   vt, blas_int* ldvt, void*   work, blas_int* lwork, double* rwork, blas_int* info);
    
    // solve system of linear equations, using LU decomposition
    void arma_fortran(arma_sgesv)(blas_int* n, blas_int* nrhs, float*  a, blas_int* lda, blas_int* ipiv, float*  b, blas_int* ldb, blas_int* info);
    void arma_fortran(arma_dgesv)(blas_int* n, blas_int* nrhs, double* a, blas_int* lda, blas_int* ipiv, double* b, blas_int* ldb, blas_int* info);
    void arma_fortran(arma_cgesv)(blas_int* n, blas_int* nrhs, void*   a, blas_int* lda, blas_int* ipiv, void*   b, blas_int* ldb, blas_int* info);
    void arma_fortran(arma_zgesv)(blas_int* n, blas_int* nrhs, void*   a, blas_int* lda, blas_int* ipiv, void*   b, blas_int* ldb, blas_int* info);
    
    // solve over/underdetermined system of linear equations
    void arma_fortran(arma_sgels)(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, float*  a, blas_int* lda, float*  b, blas_int* ldb, float*  work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_dgels)(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, double* a, blas_int* lda, double* b, blas_int* ldb, double* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_cgels)(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, void*   a, blas_int* lda, void*   b, blas_int* ldb, void*   work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_zgels)(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, void*   a, blas_int* lda, void*   b, blas_int* ldb, void*   work, blas_int* lwork, blas_int* info);
    
    // solve a triangular system of linear equations
    void arma_fortran(arma_strtrs)(char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, const float*  a, blas_int* lda, float*  b, blas_int* ldb, blas_int* info);
    void arma_fortran(arma_dtrtrs)(char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, const double* a, blas_int* lda, double* b, blas_int* ldb, blas_int* info);
    void arma_fortran(arma_ctrtrs)(char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, const void*   a, blas_int* lda, void*   b, blas_int* ldb, blas_int* info);
    void arma_fortran(arma_ztrtrs)(char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, const void*   a, blas_int* lda, void*   b, blas_int* ldb, blas_int* info);
    
    // Schur decomposition (real matrices)
    void arma_fortran(arma_sgees)(char* jobvs, char* sort, blas_int* select, blas_int* n, float*  a, blas_int* lda, blas_int* sdim, float*  wr, float*  wi, float*  vs, blas_int* ldvs, float*  work, blas_int* lwork, blas_int* bwork, blas_int* info);
    void arma_fortran(arma_dgees)(char* jobvs, char* sort, blas_int* select, blas_int* n, double* a, blas_int* lda, blas_int* sdim, double* wr, double* wi, double* vs, blas_int* ldvs, double* work, blas_int* lwork, blas_int* bwork, blas_int* info);
    
    // Schur decomposition (complex matrices)
    void arma_fortran(arma_cgees)(char* jobvs, char* sort, blas_int* select, blas_int* n, void* a, blas_int* lda, blas_int* sdim, void* w, void* vs, blas_int* ldvs, void* work, blas_int* lwork, float*  rwork, blas_int* bwork, blas_int* info);
    void arma_fortran(arma_zgees)(char* jobvs, char* sort, blas_int* select, blas_int* n, void* a, blas_int* lda, blas_int* sdim, void* w, void* vs, blas_int* ldvs, void* work, blas_int* lwork, double* rwork, blas_int* bwork, blas_int* info);
    
    // solve a Sylvester equation ax + xb = c, with a and b assumed to be in Schur form
    void arma_fortran(arma_strsyl)(char* transa, char* transb, blas_int* isgn, blas_int* m, blas_int* n, const float*  a, blas_int* lda, const float*  b, blas_int* ldb, float*  c, blas_int* ldc, float*  scale, blas_int* info);
    void arma_fortran(arma_dtrsyl)(char* transa, char* transb, blas_int* isgn, blas_int* m, blas_int* n, const double* a, blas_int* lda, const double* b, blas_int* ldb, double* c, blas_int* ldc, double* scale, blas_int* info);
    void arma_fortran(arma_ctrsyl)(char* transa, char* transb, blas_int* isgn, blas_int* m, blas_int* n, const void*   a, blas_int* lda, const void*   b, blas_int* ldb, void*   c, blas_int* ldc, float*  scale, blas_int* info);
    void arma_fortran(arma_ztrsyl)(char* transa, char* transb, blas_int* isgn, blas_int* m, blas_int* n, const void*   a, blas_int* lda, const void*   b, blas_int* ldb, void*   c, blas_int* ldc, double* scale, blas_int* info);
    
    void arma_fortran(arma_ssytrf)(char* uplo, blas_int* n, float*  a, blas_int* lda, blas_int* ipiv, float*  work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_dsytrf)(char* uplo, blas_int* n, double* a, blas_int* lda, blas_int* ipiv, double* work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_csytrf)(char* uplo, blas_int* n, void*   a, blas_int* lda, blas_int* ipiv, void*   work, blas_int* lwork, blas_int* info);
    void arma_fortran(arma_zsytrf)(char* uplo, blas_int* n, void*   a, blas_int* lda, blas_int* ipiv, void*   work, blas_int* lwork, blas_int* info);
    
    void arma_fortran(arma_ssytri)(char* uplo, blas_int* n, float*  a, blas_int* lda, blas_int* ipiv, float*  work, blas_int* info);
    void arma_fortran(arma_dsytri)(char* uplo, blas_int* n, double* a, blas_int* lda, blas_int* ipiv, double* work, blas_int* info);
    void arma_fortran(arma_csytri)(char* uplo, blas_int* n, void*   a, blas_int* lda, blas_int* ipiv, void*   work, blas_int* info);
    void arma_fortran(arma_zsytri)(char* uplo, blas_int* n, void*   a, blas_int* lda, blas_int* ipiv, void*   work, blas_int* info);
    
    // void arma_fortran(arma_dgeqp3)(blas_int* m, blas_int* n, double* a, blas_int* lda, blas_int* jpvt, double* tau, double* work, blas_int* lwork, blas_int* info);
    // void arma_fortran(arma_dormqr)(char* side, char* trans, blas_int* m, blas_int* n, blas_int* k, double* a, blas_int* lda, double* tau, double* c, blas_int* ldc, double* work, blas_int* lwork, blas_int* info);
    // void  arma_fortran(arma_dposv)(char* uplo, blas_int* n, blas_int* nrhs, double* a, blas_int* lda, double* b, blas_int* ldb, blas_int* info);
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
      arma_fortran(arma_sgetrf)(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgetrf)(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cgetrf)(m, n, (T*)a, lda, ipiv, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zgetrf)(m, n, (T*)a, lda, ipiv, info);
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
      arma_fortran(arma_sgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  trtri(char* uplo, char* diag, blas_int* n, eT* a, blas_int* lda, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(arma_strtri)(uplo, diag, n, (T*)a, lda, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dtrtri)(uplo, diag, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_ctrtri)(uplo, diag, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_ztrtri)(uplo, diag, n, (T*)a, lda, info);
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
      arma_fortran(arma_ssyev)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dsyev)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info);
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
      arma_fortran(arma_cheev)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      typedef typename std::complex<T> cx_T;
      arma_fortran(arma_zheev)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info);
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
      arma_fortran(arma_sgeev)(jobvl, jobvr, n,  (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgeev)(jobvl, jobvr, n,  (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info);
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
      arma_fortran(arma_cgeev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      typedef typename std::complex<T> cx_T;
      arma_fortran(arma_zgeev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info);
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
      arma_fortran(arma_spotrf)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dpotrf)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cpotrf)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zpotrf)(uplo, n, (T*)a, lda, info);
      }
    
    }
  
  
  
  template<typename eT>
  inline
  void
  potri(char* uplo, blas_int* n, eT* a, blas_int* lda, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(arma_spotri)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dpotri)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cpotri)(uplo, n, (T*)a, lda, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zpotri)(uplo, n, (T*)a, lda, info);
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
      arma_fortran(arma_sgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
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
      arma_fortran(arma_sorgqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dorgqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
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
      arma_fortran(arma_cungqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_zungqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info);
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
      arma_fortran(arma_sgesvd)(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgesvd)(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info);
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
      arma_fortran(arma_cgesvd)
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
      arma_fortran(arma_zgesvd)
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
      arma_fortran(arma_sgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info);
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
      arma_fortran(arma_sgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_cgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  trtrs(char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, const eT* a, blas_int* lda, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(arma_strtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dtrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_ctrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_ztrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  gees(char* jobvs, char* sort, blas_int* select, blas_int* n, eT* a, blas_int* lda, blas_int* sdim, eT* wr, eT* wi, eT* vs, blas_int* ldvs, eT* work, blas_int* lwork, blas_int* bwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(sgees)(jobvs, sort, select, n, (T*)a, lda, sdim, (T*)wr, (T*)wi, (T*)vs, ldvs, (T*)work, lwork, bwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dgees)(jobvs, sort, select, n, (T*)a, lda, sdim, (T*)wr, (T*)wi, (T*)vs, ldvs, (T*)work, lwork, bwork, info);
      }
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gees(char* jobvs, char* sort, blas_int* select, blas_int* n, std::complex<T>* a, blas_int* lda, blas_int* sdim, std::complex<T>* w, std::complex<T>* vs, blas_int* ldvs, std::complex<T>* work, blas_int* lwork, T* rwork, blas_int* bwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<T>::value == false>::apply();
    arma_type_check<is_supported_blas_type< std::complex<T> >::value == false>::apply();
    
    if(is_float<T>::value == true)
      {
      typedef float bT;
      typedef std::complex<bT> cT;
      arma_fortran(cgees)(jobvs, sort, select, n, (cT*)a, lda, sdim, (cT*)w, (cT*)vs, ldvs, (cT*)work, lwork, (bT*)rwork, bwork, info);
      }
    else
    if(is_double<T>::value == true)
      {
      typedef double bT;
      typedef std::complex<bT> cT;
      arma_fortran(zgees)(jobvs, sort, select, n, (cT*)a, lda, sdim, (cT*)w, (cT*)vs, ldvs, (cT*)work, lwork, (bT*)rwork, bwork, info);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  trsyl(char* transa, char* transb, blas_int* isgn, blas_int* m, blas_int* n, const eT* a, blas_int* lda, const eT* b, blas_int* ldb, eT* c, blas_int* ldc, eT* scale, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(strsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (T*)scale, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(dtrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (T*)scale, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(ctrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (float*)scale, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(ztrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (double*)scale, info);
      }
    }
  
  
  template<typename eT>
  inline
  void
  sytrf(char* uplo, blas_int* n, eT* a, blas_int* lda, blas_int* ipiv, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(arma_ssytrf)(uplo, n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dsytrf)(uplo, n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_csytrf)(uplo, n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zsytrf)(uplo, n, (T*)a, lda, ipiv, (T*)work, lwork, info);
      }
    }
  
  
  template<typename eT>
  inline
  void
  sytri(char* uplo, blas_int* n, eT* a, blas_int* lda, blas_int* ipiv, eT* work, blas_int* info)
    {
    arma_type_check<is_supported_blas_type<eT>::value == false>::apply();
    
    if(is_float<eT>::value == true)
      {
      typedef float T;
      arma_fortran(arma_ssytri)(uplo, n, (T*)a, lda, ipiv, (T*)work, info);
      }
    else
    if(is_double<eT>::value == true)
      {
      typedef double T;
      arma_fortran(arma_dsytri)(uplo, n, (T*)a, lda, ipiv, (T*)work, info);
      }
    else
    if(is_supported_complex_float<eT>::value == true)
      {
      typedef std::complex<float> T;
      arma_fortran(arma_csytri)(uplo, n, (T*)a, lda, ipiv, (T*)work, info);
      }
    else
    if(is_supported_complex_double<eT>::value == true)
      {
      typedef std::complex<double> T;
      arma_fortran(arma_zsytri)(uplo, n, (T*)a, lda, ipiv, (T*)work, info);
      }
    }
  
  
  //! @}
  }


#endif
