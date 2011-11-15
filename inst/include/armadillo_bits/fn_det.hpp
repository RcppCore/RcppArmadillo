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


//! \addtogroup fn_det
//! @{



//! determinant of mat
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
det
  (
  const Base<typename T1::elem_type,T1>& X,
  const bool slow = false,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return auxlib::det(X, slow);
  }



//! determinant of diagmat
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
det
  (
  const Op<T1, op_diagmat>& X,
  const bool slow = false
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(slow);
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy<T1> A(X.m);
  
  const uword A_n_elem = A.n_elem;
  
  eT val = eT(1);
  
  for(uword i=0; i<A_n_elem; ++i)
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
det
  (
  const Op<T1,op_inv>& in,
  const bool slow = false,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  eT tmp = det(in.m, slow);
  arma_warn( (tmp == eT(0)), "det(): warning: denominator is zero" );
  
  return eT(1) / tmp;
  }



//! determinant of trans(A)
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
det
  (
  const Op<T1,op_htrans>& in,
  const bool slow = false,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk1 = 0,
  const typename         arma_not_cx<typename T1::elem_type>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  return det(X, slow);
  }



//! @}
