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


//! \addtogroup auxlib
//! @{


//! wrapper for accessing external functions defined in ATLAS, LAPACK or BLAS libraries
class auxlib
  {
  public:
  
  //
  // inv
  
  template<typename eT, typename T1>
  inline static bool inv(Mat<eT>& out, const Base<eT,T1>& X);
  
  template<typename eT>
  inline static bool inv(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename eT>
  inline static bool inv_noalias_tinymat(Mat<eT>& out, const Mat<eT>& X, const u32 N);
  
  template<typename eT>
  inline static bool inv_inplace_tinymat(Mat<eT>& out, const u32 N);
  
  template<typename eT>
  inline static bool inv_inplace_lapack(Mat<eT>& out);
  
  
  //
  // det
  
  template<typename eT, typename T1>
  inline static eT det(const Base<eT,T1>& X);
  
  template<typename eT>
  inline static eT det_tinymat(const Mat<eT>& X, const u32 N);
  
  template<typename eT>
  inline static eT det_lapack(const Mat<eT>& X, const bool make_copy);
  
  
  //
  // log_det
  
  template<typename eT, typename T1>
  inline static void log_det(eT& out_val, typename get_pod_type<eT>::result& out_sign, const Base<eT,T1>& X);
  
  
  //
  // lu
  
  template<typename eT, typename T1>
  inline static void lu(Mat<eT>& L, Mat<eT>& U, podarray<blas_int>& ipiv, const Base<eT,T1>& X);
  
  template<typename eT, typename T1>
  inline static void lu(Mat<eT>& L, Mat<eT>& U, Mat<eT>& P, const Base<eT,T1>& X);
  
  template<typename eT, typename T1>
  inline static void lu(Mat<eT>& L, Mat<eT>& U, const Base<eT,T1>& X);
  
  
  //
  // eig
  
  template<typename eT, typename T1> 
  inline static bool eig_sym(Col<eT>& eigval, const Base<eT,T1>& X);
  
  template<typename T, typename T1> 
  inline static bool eig_sym(Col<T>& eigval, const Base<std::complex<T>,T1>& X);
  
  template<typename eT, typename T1>
  inline static bool eig_sym(Col<eT>& eigval, Mat<eT>& eigvec, const Base<eT,T1>& X);
  
  template<typename T, typename T1>
  inline static bool eig_sym(Col<T>& eigval, Mat< std::complex<T> >& eigvec, const Base<std::complex<T>,T1>& X);
  
  template<typename T, typename T1>
  inline static bool eig_gen(Col< std::complex<T> >& eigval, Mat<T>& l_eigvec, Mat<T>& r_eigvec, const Base<T,T1>& X, const char side);
  
  template<typename T, typename T1>
  inline static bool eig_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& l_eigvec, Mat< std::complex<T> >& r_eigvec, const Base< std::complex<T>, T1 >& X, const char side);
  
  
  //
  // chol
  
  template<typename eT, typename T1>
  inline static bool chol(Mat<eT>& out, const Base<eT,T1>& X);
  
  
  //
  // qr
  
  template<typename eT, typename T1>
  inline static bool qr(Mat<eT>& Q, Mat<eT>& R, const Base<eT,T1>& X);
  
  
  //
  // svd
  
  template<typename eT, typename T1>
  inline static bool svd(Col<eT>& S, const Base<eT,T1>& X, u32& n_rows, u32& n_cols);
  
  template<typename T, typename T1>
  inline static bool svd(Col<T>& S, const Base<std::complex<T>, T1>& X, u32& n_rows, u32& n_cols);
  
  template<typename eT, typename T1>
  inline static bool svd(Col<eT>& S, const Base<eT,T1>& X);
  
  template<typename T, typename T1>
  inline static bool svd(Col<T>& S, const Base<std::complex<T>, T1>& X);
  
  template<typename eT, typename T1>
  inline static bool svd(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, const Base<eT,T1>& X);
  
  template<typename T, typename T1>
  inline static bool svd(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, const Base< std::complex<T>, T1>& X);
  
  
  //
  // solve
  
  template<typename eT>
  inline static bool solve(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename eT>
  inline static bool solve_od(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename eT>
  inline static bool solve_ud(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  };


//! @}
