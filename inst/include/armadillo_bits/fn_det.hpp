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


//! \addtogroup fn_det
//! @{



//! determinant of mat
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
det(const Base<typename T1::elem_type,T1>& X, const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( !A.is_square(), "det(): matrix must be square" );
  
  return auxlib::det(A);
  }



//! determinant of diagmat
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
det(const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy<T1> A(X.m);
  
  arma_debug_check( (A.n_elem == 0), "det(): given object has no elements" );
  
  eT val = A[0];
  
  for(u32 i=1; i<A.n_elem; ++i)
    {
    val *= A[i];
    }
  
  return val;
  }



//! determinant of inv(A), without doing the inverse operation
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
det(const Op<T1,op_inv>& in, const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  isnt_fltpt<eT>::check();
  
  eT tmp = det(in.m);
  arma_warn( (tmp == eT(0)), "det(): warning: denominator is zero" );
  
  return eT(1) / tmp;
  }



//! determinant of trans(A)
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
det(const Op<T1,op_trans>& in, const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& X = tmp.M;

  return det(X);
  }



//! @}
