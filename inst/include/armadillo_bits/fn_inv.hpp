// Copyright (C) 2008-2011 Conrad Sanderson
// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_inv
//! @{



//! delayed matrix inverse (general matrices)
template<typename T1>
arma_inline
const Op<T1, op_inv>
inv
  (
  const Base<typename T1::elem_type,T1>& X,
  const bool slow = false,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_inv>(X.get_ref(), ((slow == false) ? 0 : 1), 0);
  }



//! remove the inverse operation if applied twice consecutively
template<typename T1>
arma_inline
const T1&
inv(const Op<T1, op_inv>& X, const bool slow = false)
  {
  arma_extra_debug_sigprint();
  arma_ignore(slow);
  
  return X.m;
  }



//! delayed matrix inverse (triangular matrices)
template<typename T1>
arma_inline
const Op<T1, op_inv_tr>
inv
  (
  const Op<T1, op_trimat>& X,
  const bool slow = false,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(slow);
  arma_ignore(junk);
  
  return Op<T1, op_inv_tr>(X.m, X.aux_uword_a, 0);
  }



//! delayed matrix inverse (symmetric positive definite matrices)
template<typename T1>
arma_inline
const Op<T1, op_inv_sympd>
inv
  (
  const Op<T1, op_sympd>& X,
  const bool slow = false,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(slow);
  arma_ignore(junk);
  
  return Op<T1, op_inv_sympd>(X.m, 0, 0);
  }



template<typename T1>
inline
bool
inv
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& X,
  const bool slow = false,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  try
    {
    out = inv(X,slow);
    }
  catch(std::runtime_error&)
    {
    return false;
    }
  
  return true;
  }



//! @}
