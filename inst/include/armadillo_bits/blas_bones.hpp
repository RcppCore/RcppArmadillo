// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



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


#endif
