// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_logmat
//! @{



template<typename T1>
arma_inline
typename enable_if2< (is_supported_blas_type<typename T1::elem_type>::value && is_cx<typename T1::elem_type>::no), const mtOp<std::complex<typename T1::elem_type>, T1, op_logmat> >::result
logmat(const Base<typename T1::elem_type,T1>& X, const uword n_iters = 100u)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<std::complex<typename T1::elem_type>, T1, op_logmat>(X.get_ref(), n_iters, uword(0));
  }



template<typename T1>
arma_inline
typename enable_if2< (is_supported_blas_type<typename T1::elem_type>::value && is_cx<typename T1::elem_type>::yes), const Op<T1, op_logmat_cx> >::result
logmat(const Base<typename T1::elem_type,T1>& X, const uword n_iters = 100u)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_logmat_cx>(X.get_ref(), n_iters, uword(0));
  }



template<typename T1>
inline
typename enable_if2< (is_supported_blas_type<typename T1::elem_type>::value && is_cx<typename T1::elem_type>::no), bool >::result
logmat(Mat< std::complex<typename T1::elem_type> >& Y, const Base<typename T1::elem_type,T1>& X, const uword n_iters = 100u)
  {
  arma_extra_debug_sigprint();
  
  return op_logmat::apply_direct(Y, X.get_ref(), n_iters);
  }



template<typename T1>
inline
typename enable_if2< (is_supported_blas_type<typename T1::elem_type>::value && is_cx<typename T1::elem_type>::yes), bool >::result
logmat(Mat<typename T1::elem_type>& Y, const Base<typename T1::elem_type,T1>& X, const uword n_iters = 100u)
  {
  arma_extra_debug_sigprint();
  
  return op_logmat_cx::apply_direct(Y, X.get_ref(), n_iters);
  }



//! @}
