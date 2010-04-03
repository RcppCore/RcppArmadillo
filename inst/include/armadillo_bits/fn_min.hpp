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


//! \addtogroup fn_min
//! @{

//! \brief
//! Delayed 'minimum values' operation.
//! The dimension, along which the minima are found, is set via 'dim'.
//! For dim = 0, the maximum value of each column is found (i.e. searches by traversing across rows).
//! For dim = 1, the maximum value of each row is found (i.e. searches by traversing across columns).
//! The default is dim = 0.
//! NOTE: This function works differently than in Matlab/Octave.

template<typename T1>
arma_inline
const Op<T1, op_min>
min(const Base<typename T1::elem_type,T1>& X, const u32 dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_min>(X.get_ref(), dim, 0);
  }


//! Immediate 'find the minimum value in a row vector' operation
template<typename eT>
inline
arma_warn_unused
eT
min(const Row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "min(): given vector has no elements" );
  
  return op_min::direct_min(A.mem, A.n_elem);
  }



//! Immediate 'find the minimum value in a column vector'
template<typename eT>
inline
arma_warn_unused
eT
min(const Col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "min(): given vector has no elements" );
  
  return op_min::direct_min(A.mem, A.n_elem);
  }



//! \brief
//! Immediate 'find minimum value' operation,
//! invoked, for example, by: min(min(A))
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
min(const Op<T1, op_min>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("min(): two consecutive min() calls detected");
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.m);
  const Mat<eT>& X = tmp1.M;
  
  arma_debug_check( (X.n_elem == 0), "min(): given matrix has no elements" );
  
  return op_min::direct_min(X.mem, X.n_elem);
  }



template<typename T1>
inline
const Op< Op<T1, op_min>, op_min>
min(const Op<T1, op_min>& in, const u32 dim)
  {
  arma_extra_debug_sigprint();
  
  return Op< Op<T1, op_min>, op_min>(in, dim, 0);
  }



template<typename eT>
inline
arma_warn_unused
eT
min(const subview_row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "min(): given vector has no elements" );
  
  return op_min::direct_min(A);
  }



template<typename eT>
inline
arma_warn_unused
eT
min(const subview_col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "min(): given vector has no elements" );
  
  return op_min::direct_min(A);
  }



template<typename eT>
inline
arma_warn_unused
eT
min(const diagview<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "min(): given vector has no elements" );
  
  return op_min::direct_min(A);
  }



template<typename eT>
inline
arma_warn_unused
eT
min(const Op<subview<eT>, op_min>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("max(): two consecutive max() calls detected");
  
  const subview<eT>& X = in.m;
  
  arma_debug_check( (X.n_elem == 0), "max(): given matrix has no elements" );
  
  return op_min::direct_min(X);
  }



//! @}
