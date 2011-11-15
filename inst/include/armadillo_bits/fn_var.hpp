// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_var
//! @{



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_var>
var(const Base<typename T1::elem_type,T1>& X, const uword norm_type = 0, const uword dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_var>(X.get_ref(), norm_type, dim);
  }



//! Immediate 'find the variance of a row vector' operation
template<typename eT>
inline
arma_warn_unused
typename get_pod_type<eT>::result
var(const Row<eT>& A, const uword norm_type = 0)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_elem = A.n_elem;
  
  arma_debug_check( (A_n_elem == 0), "var(): given object has no elements" );
  
  return op_var::direct_var(A.mem, A_n_elem, norm_type);
  }



//! Immediate 'find the variance of a column vector' operation
template<typename eT>
inline
arma_warn_unused
typename get_pod_type<eT>::result
var(const Col<eT>& A, const uword norm_type = 0)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_elem = A.n_elem;
  
  arma_debug_check( (A_n_elem == 0), "var(): given object has no elements" );
  
  return op_var::direct_var(A.mem, A_n_elem, norm_type);
  }



template<typename eT>
inline
arma_warn_unused
typename get_pod_type<eT>::result
var(const subview_row<eT>& A, const uword norm_type = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "var(): given object has no elements" );
  
  return op_var::direct_var(A, norm_type);
  }



template<typename eT>
inline
arma_warn_unused
typename get_pod_type<eT>::result
var(const subview_col<eT>& A, const uword norm_type = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "var(): given object has no elements" );
  
  return op_var::direct_var(A.colptr(0), A.n_rows, norm_type);
  }



template<typename eT>
inline
arma_warn_unused
typename get_pod_type<eT>::result
var(const diagview<eT>& A, const uword norm_type = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "var(): given object has no elements" );
  
  return op_var::direct_var(A, norm_type);
  }



template<typename eT, typename T1>
inline
arma_warn_unused
typename get_pod_type<eT>::result
var(const subview_elem1<eT,T1>& A, const uword norm_type = 0)
  {
  arma_extra_debug_sigprint();
  
  const Col<eT> X(A);
  
  return var(X, norm_type);
  }



//! @}
