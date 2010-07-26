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


//! \addtogroup fn_inv
//! @{



//! delayed matrix inverse
template<typename T1>
arma_inline
const Op<T1, op_inv>
inv(const Base<typename T1::elem_type,T1>& X, const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_inv>(X.get_ref());
  }



//! remove the inverse operation if applied twice consecutively
template<typename T1>
arma_inline
const T1&
inv(const Op<T1, op_inv>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.m;
  }



//! @}
