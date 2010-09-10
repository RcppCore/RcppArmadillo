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


//! \addtogroup fn_accu
//! @{



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_unwrap(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  return arrayops::accumulate( A.memptr(), A.n_elem );
  }



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_proxy(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.get_ref());
  
  const u32 N = A.n_elem;
  
  eT val = eT(0);
  
  for(u32 i=0; i<N; ++i)
    {
    val += A[i];
    }
  
  return val;
  }



//! accumulate the elements of a matrix
template<typename T1>
arma_inline
arma_warn_unused
typename T1::elem_type
accu(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return (is_Mat<T1>::value == true) ? accu_unwrap(X) : accu_proxy(X);
  }



//! explicit handling of Hamming norm (also known as zero norm)
template<typename T1>
arma_inline
arma_warn_unused
u32
accu(const mtOp<u32,T1,op_rel_noteq>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  const u32 n_elem = A.n_elem;
  const eT  val    = X.aux;
  
  u32 n_nonzero = 0;
  for(u32 i=0; i<n_elem; ++i)
    {
    if(A[i] != val)
      {
      ++n_nonzero;
      }
    }
  
  return n_nonzero;
  }



//! accumulate the elements of a cube
template<typename T1>
arma_hot
arma_warn_unused
inline
typename T1::elem_type
accu(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const ProxyCube<T1> A(X.get_ref());
  
  const u32 n_elem = A.n_elem;
  
  eT val = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    val += A[i];
    }
  
  return val;
  }



//! accumulate the elements of a diagview
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();  
  
  const u32 n_elem = X.n_elem;
  eT val = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    val += X[i];
    }
  
  return val;
  }



//! accumulate the elements of a subview (submatrix)
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const subview<eT>& S)
  {
  arma_extra_debug_sigprint();  
  
  const u32 S_n_rows = S.n_rows;
  const u32 S_n_cols = S.n_cols;
  
  eT val = eT(0);
  
  for(u32 col=0; col<S_n_cols; ++col)
    {
    val += arrayops::accumulate( S.colptr(col), S_n_rows );
    }
  
  return val;
  }



//! accumulate the elements of a subview_row
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const subview_row<eT>& S)
  {
  arma_extra_debug_sigprint();  
  
  const Mat<eT>& X = S.m;
  
  const u32 row       = S.aux_row1;
  const u32 start_col = S.aux_col1;
  const u32 end_col   = S.aux_col2;
  
  eT val = eT(0);
  
  for(u32 col=start_col; col<=end_col; ++col)
    {
    val += X.at(row,col);
    }
  
  return val;
  }



//! accumulate the elements of a subview_col
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const subview_col<eT>& S)
  {
  arma_extra_debug_sigprint();
  
  return arrayops::accumulate( S.colptr(0), S.n_rows );
  }



//! @}
