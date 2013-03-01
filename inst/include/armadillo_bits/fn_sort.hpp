// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_sort
//! @{


template<typename T1>
arma_inline
typename
enable_if2
  <
  ( (is_arma_type<T1>::value == true) && (resolves_to_vector<T1>::value == false) ),
  const Op<T1, op_sort>
  >::result
sort
  (
  const T1&   X,
  const uword sort_type = 0,
  const uword dim       = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_sort>(X, sort_type, dim);
  }



template<typename T1>
arma_inline
typename
enable_if2
  <
  ( (is_arma_type<T1>::value == true) && (resolves_to_vector<T1>::value == true) ),
  const Op<T1, op_sort>
  >::result
sort
  (
  const T1&   X,
  const uword sort_type = 0
  )
  {
  arma_extra_debug_sigprint();
  
  const uword dim = (T1::is_col) ? 0 : 1;
  
  return Op<T1, op_sort>(X, sort_type, dim);
  }



//! @}
