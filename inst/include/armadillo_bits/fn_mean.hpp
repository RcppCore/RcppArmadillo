// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_mean
//! @{



template<typename T1>
arma_inline
const Op<T1, op_mean>
mean(const Base<typename T1::elem_type,T1>& X, const u32 dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_mean>(X.get_ref(), dim, 0);
  }



//! Immediate 'find the mean value of a row vector' operation
template<typename eT>
inline
arma_warn_unused
eT
mean(const Row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "mean(): given vector has no elements" );
  
  return op_mean::direct_mean(A.mem, A.n_elem);
  }



//! Immediate 'find the mean value of a column vector' operation
template<typename eT>
inline
arma_warn_unused
eT
mean(const Col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "mean(): given vector has no elements" );
  
  return op_mean::direct_mean(A.mem, A.n_elem);
  }



//! \brief
//! Immediate 'find mean value' operation,
//! invoked, for example, by: mean(mean(A))
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
mean(const Op<T1, op_mean>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("mean(): two consecutive mean() calls detected");
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.m);
  const Mat<eT>& X = tmp1.M;
  
  arma_debug_check( (X.n_elem == 0), "mean(): given matrix has no elements" );
  
  return op_mean::direct_mean(X.mem, X.n_elem);
  }



template<typename T1>
arma_inline
const Op< Op<T1, op_mean>, op_mean>
mean(const Op<T1, op_mean>& in, const u32 dim)
  {
  arma_extra_debug_sigprint();
  
  return Op< Op<T1, op_mean>, op_mean>(in, dim, 0);
  }



template<typename eT>
inline
arma_warn_unused
eT
mean(const subview_row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "mean(): given vector has no elements" );
  
  return op_mean::direct_mean(A);
  }



template<typename eT>
inline
arma_warn_unused
eT
mean(const subview_col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "mean(): given vector has no elements" );
  
  return op_mean::direct_mean(A);
  }



template<typename eT>
inline
arma_warn_unused
eT
mean(const Op<subview<eT>, op_mean>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("mean(): two consecutive mean() calls detected");
  
  const subview<eT>& X = in.m;
  
  arma_debug_check( (X.n_elem == 0), "mean(): given matrix has no elements" );
  
  return op_mean::direct_mean(X);
  }



template<typename eT>
inline
arma_warn_unused
eT
mean(const diagview<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "mean(): given vector has no elements" );
  
  return op_mean::direct_mean(A);
  }



//! @}
