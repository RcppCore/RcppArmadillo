// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_max
//! @{


//! \brief
//! Delayed 'maximum values' operation.
//! The dimension, along which the maxima are found, is set via 'dim'.
//! For dim = 0, the maximum value of each column is found (i.e. searches by traversing across rows).
//! For dim = 1, the maximum value of each row is found (i.e. searches by traversing across columns).
//! The default is dim = 0.

template<typename T1>
arma_inline
const Op<T1, op_max>
max(const Base<typename T1::elem_type,T1>& X, const uword dim = 0)
  {
  arma_extra_debug_sigprint();

  return Op<T1, op_max>(X.get_ref(), dim, 0);
  }


//! Immediate 'find the maximum value in a row vector' operation
template<typename eT>
inline
arma_warn_unused
eT
max(const Row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_elem = A.n_elem;
  
  arma_debug_check( (A_n_elem == 0), "max(): given object has no elements" );
  
  return op_max::direct_max(A.mem, A_n_elem);
  }



//! Immediate 'find the maximum value in a column vector' operation
template<typename eT>
inline
arma_warn_unused
eT
max(const Col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_elem = A.n_elem;
  
  arma_debug_check( (A_n_elem == 0), "max(): given object has no elements" );
  
  return op_max::direct_max(A.mem, A_n_elem);
  }



//! \brief
//! Immediate 'find maximum value' operation,
//! invoked, for example, by: max(max(A))
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
max(const Op<T1, op_max>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("max(): two consecutive max() calls detected");
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.m);
  const Mat<eT>& X = tmp1.M;
  
  const uword X_n_elem = X.n_elem;
  
  arma_debug_check( (X_n_elem == 0), "max(): given object has no elements" );
  
  return op_max::direct_max(X.mem, X_n_elem);
  }



template<typename T1>
arma_inline
const Op< Op<T1, op_max>, op_max>
max(const Op<T1, op_max>& in, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op< Op<T1, op_max>, op_max>(in, dim, 0);
  }



template<typename eT>
inline
arma_warn_unused
eT
max(const subview_row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "max(): given object has no elements" );
  
  return op_max::direct_max(A);
  }



template<typename eT>
inline
arma_warn_unused
eT
max(const subview_col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "max(): given object has no elements" );
  
  return op_max::direct_max(A.colptr(0), A.n_rows);
  }



template<typename eT>
inline
arma_warn_unused
eT
max(const Op<subview<eT>, op_max>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("max(): two consecutive max() calls detected");
  
  const subview<eT>& X = in.m;
  
  arma_debug_check( (X.n_elem == 0), "max(): given object has no elements" );
  
  return op_max::direct_max(X);
  }



template<typename eT>
inline
arma_warn_unused
eT
max(const diagview<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "max(): given object has no elements" );
  
  return op_max::direct_max(A);
  }



template<typename eT, typename T1>
inline
arma_warn_unused
eT
max(const subview_elem1<eT,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> X(A);
  
  const uword X_n_elem = X.n_elem;
  
  arma_debug_check( (X_n_elem == 0), "max(): given object has no elements" );
  
  return op_max::direct_max(X.mem, X_n_elem);
  }



//! @}
