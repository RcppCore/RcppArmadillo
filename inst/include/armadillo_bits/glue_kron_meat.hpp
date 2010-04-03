// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup glue_kron
//! @{



//! \brief
//! both input matrices have the same element type
template<typename eT>
inline
void
glue_kron::direct_kron(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  const u32 A_rows = A.n_rows;
  const u32 A_cols = A.n_cols;
  const u32 B_rows = B.n_rows;
  const u32 B_cols = B.n_cols;
  
  out.set_size(A_rows*B_rows, A_cols*B_cols);
  
  for(u32 i = 0; i < A_rows; i++)
    {
    for(u32 j = 0; j < A_cols; j++)
      {
      out.submat(i*B_rows, j*B_cols, (i+1)*B_rows-1, (j+1)*B_cols-1) = A(i,j) * B; 
      }
    }  
  }



//! \brief
//! different types of input matrices
//! A -> complex, B -> basic element type
template<typename T>
inline
void
glue_kron::direct_kron(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A, const Mat<T>& B)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const u32 A_rows = A.n_rows;
  const u32 A_cols = A.n_cols;
  const u32 B_rows = B.n_rows;
  const u32 B_cols = B.n_cols;
  
  out.set_size(A_rows*B_rows, A_cols*B_cols);
  
  Mat<eT> tmp_B = conv_to< Mat<eT> >::from(B);
  
  for(u32 i = 0; i < A_rows; i++)
    {
    for(u32 j = 0; j < A_cols; j++)
      {
      out.submat(i*B_rows, j*B_cols, (i+1)*B_rows-1, (j+1)*B_cols-1) = A(i,j) * tmp_B; 
      }
    }  
  }



//! \brief
//! different types of input matrices
//! A -> basic element type, B -> complex
template<typename T>
inline
void
glue_kron::direct_kron(Mat< std::complex<T> >& out, const Mat<T>& A, const Mat< std::complex<T> >& B)
  {
  arma_extra_debug_sigprint();
  
  const u32 A_rows = A.n_rows;
  const u32 A_cols = A.n_cols;
  const u32 B_rows = B.n_rows;
  const u32 B_cols = B.n_cols;
  
  out.set_size(A_rows*B_rows, A_cols*B_cols);
  
  for(u32 i = 0; i < A_rows; i++)
    {
    for(u32 j = 0; j < A_cols; j++)
      {
      out.submat(i*B_rows, j*B_cols, (i+1)*B_rows-1, (j+1)*B_cols-1) = A(i,j) * B; 
      }
    }  
  }



//! \brief
//! apply Kronecker product for two objects with same element type
template<typename T1, typename T2>
inline
void
glue_kron::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_kron>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> A_tmp(X.A, out);
  const unwrap_check<T2> B_tmp(X.B, out);
  
  const Mat<eT>& A = A_tmp.M;
  const Mat<eT>& B = B_tmp.M;
  
  glue_kron::direct_kron(out, A, B); 
  }



//! @}
