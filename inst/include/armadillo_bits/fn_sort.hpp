// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


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
