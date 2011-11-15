// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_sum
//! @{


//! \brief
//! Delayed sum of elements of a matrix along a specified dimension (either rows or columns).
//! The result is stored in a dense matrix that has either one column or one row.
//! For dim = 0, find the sum of each column.
//! For dim = 1, find the sum of each row.
//! The default is dim = 0.
//! NOTE: this function works differently than in Matlab/Octave.

template<typename T1>
arma_inline
const Op<T1, op_sum>
sum(const Base<typename T1::elem_type,T1>& X, const uword dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_sum>(X.get_ref(), dim, 0);
  }


//! \brief
//! Immediate 'sum all values' operation for a row vector
template<typename eT>
inline
arma_warn_unused
eT
sum(const Row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return accu(X);
  }



//! \brief
//! Immediate 'sum all values' operation for a column vector
template<typename eT>
inline
arma_warn_unused
eT
sum(const Col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return accu(X);
  }



//! \brief
//! Immediate 'sum all values' operation,
//! invoked, for example, by: sum(sum(A))

template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
sum(const Op<T1, op_sum>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("sum(): two consecutive sum() calls detected");
  
  return accu(in.m);
  }



template<typename T1>
arma_inline
const Op<Op<T1, op_sum>, op_sum>
sum(const Op<T1, op_sum>& in, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op<Op<T1, op_sum>, op_sum>(in, dim, 0);
  }



//! sum all values of a subview_row
template<typename eT>
inline
arma_warn_unused
eT
sum(const subview_row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return accu(X);
  }



//! sum all values of a subview_col
template<typename eT>
inline
arma_warn_unused
eT
sum(const subview_col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return accu(X);
  }



//! sum all values of a diagview
template<typename eT>
inline
arma_warn_unused
eT
sum(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return accu(X);
  }



template<typename eT, typename T1>
inline
arma_warn_unused
eT
sum(const subview_elem1<eT,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return accu(A);
  }



//! @}
