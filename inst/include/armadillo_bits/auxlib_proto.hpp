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


//! \addtogroup auxlib
//! @{


//! wrapper for accessing external functions defined in ATLAS, LAPACK or BLAS libraries
class auxlib
  {
  public:
  
  //
  // inv
  
  template<typename eT>
  inline static bool inv_noalias(Mat<eT>& out, const Mat<eT>& X);
  
  template<typename eT>
  inline static bool inv_inplace(Mat<eT>& X);
  
  
  //
  // det
  
  template<typename eT>
  inline static eT det(const Mat<eT>& X);
  
  template<typename eT>
  inline static void log_det(eT& out_val, typename get_pod_type<eT>::result& out_sign, const Mat<eT>& X);
  
  
  //
  // lu
  
  template<typename eT>
  inline static void lu(Mat<eT>& L, Mat<eT>& U, podarray<int>& ipiv, const Mat<eT>& X_orig);
  
  template<typename eT>
  inline static void lu(Mat<eT>& L, Mat<eT>& U, Mat<eT>& P, const Mat<eT>& X);
  
  template<typename eT>
  inline static void lu(Mat<eT>& L, Mat<eT>& U, const Mat<eT>& X);
  
  
  //
  // eig
  
  template<typename eT> 
  inline static void eig_sym(Col<eT>& eigval, const Mat<eT>& A);
  
  template<typename T> 
  inline static void eig_sym(Col<T>& eigval, const Mat< std::complex<T> >& A);

  template<typename eT>
  inline static void eig_sym(Col<eT>& eigval, Mat<eT>& eigvec, const Mat<eT>& A);
  
  template<typename T>
  inline static void eig_sym(Col<T>& eigval, Mat< std::complex<T> >& eigvec, const Mat< std::complex<T> >& A);

  template<typename eT>
  inline static void eig_gen(Col< std::complex<eT> >& eigval, Mat<eT>& l_eigvec, Mat<eT>& r_eigvec, const Mat<eT>& A, const char side);

  template<typename T>
  inline static void eig_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& l_eigvec, Mat< std::complex<T> >& r_eigvec, const Mat< std::complex<T> >& A, const char side);
  
  //
  // chol
  
  template<typename eT>
  inline static bool chol(Mat<eT>& out, const Mat<eT>& X);
  
  
  //
  // qr
  
  template<typename eT>
  inline static bool qr(Mat<eT>& Q, Mat<eT>& R, const Mat<eT>& X);
  
  
  //
  // svd
  
  template<typename eT>
  inline static bool svd(Col<eT>& S, const Mat<eT>& X);
  
  template<typename T>
  inline static bool svd(Col<T>& S, const Mat< std::complex<T> >& X);
  
  template<typename eT>
  inline static bool svd(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, const Mat<eT>& X);
  
  template<typename T>
  inline static bool svd(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, const Mat< std::complex<T> >& X);
  
  
  // solve
  
  template<typename eT>
  inline static bool solve(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename eT>
  inline static bool solve_od(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename eT>
  inline static bool solve_ud(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  };


//! @}
