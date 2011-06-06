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
  
  return Op<T1, op_inv_tr>(X.m, X.aux_u32_a, 0);
  }



//! delayed matrix inverse (symmetric matrices)
template<typename T1>
arma_inline
const Op<T1, op_inv_sym>
inv
  (
  const Op<T1, op_symmat>& X,
  const bool slow = false,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(slow);
  arma_ignore(junk);
  
  return Op<T1, op_inv_sym>(X.m, X.aux_u32_a, 0);
  }



//! @}
