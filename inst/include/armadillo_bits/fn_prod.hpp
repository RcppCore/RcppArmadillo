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


//! \addtogroup fn_prod
//! @{


//! \brief
//! Delayed product of elements of a matrix along a specified dimension (either rows or columns).
//! The result is stored in a dense matrix that has either one column or one row.
//! For dim = 0, find the sum of each column (i.e. traverse across rows)
//! For dim = 1, find the sum of each row (i.e. traverse across columns)
//! The default is dim = 0.
//! NOTE: this function works differently than in Matlab/Octave.

template<typename T1>
arma_inline
const Op<T1, op_prod>
prod(const Base<typename T1::elem_type,T1>& X, const u32 dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_prod>(X.get_ref(), dim, 0);
  }



//! \brief
//! Immediate 'product of all values' operation for a row vector
template<typename eT>
inline
eT
prod(const Row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (X.n_elem < 1), "prod(): given object has no elements" );
  
  const u32 n_elem = X.n_elem;
  const eT* X_mem  = X.memptr();
  
  eT val = X_mem[0];
  
  for(u32 i=1; i<n_elem; ++i)
    {
    val *= X_mem[i];
    }
  
  return val;
  }



//! \brief
//! Immediate 'product of all values' operation for a column vector
template<typename eT>
inline
eT
prod(const Col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (X.n_elem < 1), "prod(): given object has no elements" );
  
  const u32 n_elem = X.n_elem;
  const eT* X_mem  = X.memptr();
  
  eT val = X_mem[0];
  
  for(u32 i=1; i<n_elem; ++i)
    {
    val *= X_mem[i];
    }
  
  return val;
  }



//! \brief
//! Immediate 'product of all values' operation,
//! invoked, for example, by: prod(prod(A))

template<typename T1>
inline
typename T1::elem_type
prod(const Op<T1, op_prod>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("prod(): two consecutive prod() calls detected");
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem < 1), "prod(): given object has no elements" );
  
  const u32 n_elem = X.n_elem;
  const eT* X_mem  = X.memptr();
  
  eT val = X_mem[0];
  
  for(u32 i=1; i<n_elem; ++i)
    {
    val *= X_mem[i];
    }
  
  return val;
  }



template<typename T1>
inline
const Op<Op<T1, op_prod>, op_prod>
prod(const Op<T1, op_prod>& in, const u32 dim)
  {
  arma_extra_debug_sigprint();
  
  return Op<Op<T1, op_prod>, op_prod>(in, dim, 0);
  }



//! product of all values of a subview_row
template<typename eT>
inline
eT
prod(const subview_row<eT>& S)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (S.n_elem < 1), "prod(): given object has no elements" );
  
  const Mat<eT>& X = S.m;
  
  const u32 row       = S.aux_row1;
  const u32 start_col = S.aux_col1;
  const u32 end_col   = S.aux_col2;
  
  eT val = X.at(row,start_col);
  
  for(u32 col=start_col+1; col<=end_col; ++col)
    {
    val *= X.at(row,col);
    }
  
  return val;
  }



//! product of all values of a subview_col
template<typename eT>
inline
eT
prod(const subview_col<eT>& S)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (S.n_elem < 1), "prod(): given object has no elements" );
  
  const eT* S_colptr = S.colptr(0);
  const u32 n_rows   = S.n_rows;
  
  eT val = S_colptr[0];
  
  for(u32 row=1; row<n_rows; ++row)
    {
    val *= S_colptr[row];
    }
  
  return val;
  }



//! product of all values of a diagview
template<typename eT>
inline
eT
prod(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (X.n_elem < 1), "prod(): given object has no elements" );
  
  const u32 n_elem = X.n_elem;
  
  eT val = X[0];
  
  for(u32 i=1; i<n_elem; ++i)
    {
    val *= X[i];
    }
  
  return val;
  }



//! @}
