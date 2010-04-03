// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_solve
//! @{



//! Solve a system of linear equations, i.e., A*X = B, where X is unknown.
//! For a square matrix A, this function is conceptually the same as X = inv(A)*B,
//! but is done more efficiently.
//! The number of rows in A and B must be the same.
//! B can be either a column vector or a matrix.
//! This function will also try to provide approximate solutions
//! to under-determined as well as over-determined systems (non-square A matrices).

template<typename eT, typename T1, typename T2>
inline
bool
solve(Mat<eT>& X, const Base<eT,T1>& A_in, const Base<eT,T2>& B_in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(A_in.get_ref());
  const unwrap<T2> tmp2(B_in.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in A and B must be the same" );
  
  bool status;
  
  if(A.n_rows == A.n_cols)
    {
    status = auxlib::solve(X, A, B);
    }
  else
  if(A.n_rows > A.n_cols)
    {
    arma_extra_debug_print("solve(): detected over-determined system");
    status = auxlib::solve_od(X, A, B);
    }
  else
    {
    arma_extra_debug_print("solve(): detected under-determined system");
    status = auxlib::solve_ud(X, A, B);
    }
  
  if(status == false)
    {
    X.reset();
    }
  
  return status;
  }



template<typename eT, typename T1, typename T2>
inline
Mat<eT>
solve(const Base<eT,T1>& A_in, const Base<eT,T2>& B_in)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> X;
  const bool status = solve(X, A_in, B_in);
  
  if(status == false)
    {
    arma_print("solve(): solution not found");
    
    X.reset();
    }
  
  return X;
  }



//! @}
