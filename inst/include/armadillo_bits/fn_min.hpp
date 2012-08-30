// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
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


//! \addtogroup fn_min
//! @{

//! \brief
//! Delayed 'minimum values' operation.
//! The dimension, along which the minima are found, is set via 'dim'.
//! For dim = 0, the minimum value of each column is found (i.e. searches by traversing across rows).
//! For dim = 1, the minimum value of each row is found (i.e. searches by traversing across columns).
//! The default is dim = 0.

template<typename T1>
arma_inline
const Op<T1, op_min>
min
  (
  const T1& X,
  const uword dim = 0,
  const typename enable_if< is_arma_type<T1>::value       == true  >::result* junk1 = 0,
  const typename enable_if< resolves_to_vector<T1>::value == false >::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return Op<T1, op_min>(X, dim, 0);
  }


template<typename T1>
arma_inline
const Op<T1, op_min>
min
  (
  const T1& X,
  const uword dim,
  const typename enable_if<resolves_to_vector<T1>::value == true>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_min>(X, dim, 0);
  }



template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
min
  (
  const T1& X,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename enable_if<resolves_to_vector<T1>::value == true>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return op_min::min(X);
  }



//! \brief
//! Immediate 'find minimum value' operation,
//! invoked, for example, by: min(min(A))
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
min(const Op<T1, op_min>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("min(): two consecutive min() calls detected");
  
  return op_min::min(in.m);
  }



template<typename T1>
arma_inline
const Op< Op<T1, op_min>, op_min>
min(const Op<T1, op_min>& in, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op< Op<T1, op_min>, op_min>(in, dim, 0);
  }



template<typename T>
arma_inline
arma_warn_unused
const typename arma_scalar_only<T>::result &
min(const T& x)
  {
  return x;
  }



template<typename T1>
inline
arma_warn_unused
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value == true) && (resolves_to_sparse_vector<T1>::value == true),
  typename T1::elem_type
  >::result
min(const T1& x)
  {
  arma_extra_debug_sigprint();
  
  return spop_min::vector_min(x);
  }



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value == true) && (resolves_to_sparse_vector<T1>::value == false),
  const SpOp<T1, spop_min>
  >::result
min(const T1& X, const uword dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_min>(X, dim, 0);
  }



template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
min(const SpOp<T1, spop_min>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("min(): two consecutive min() calls detected");
  
  return spop_min::vector_min(X.m);
  }



template<typename T1>
inline
const SpOp< SpOp<T1, spop_min>, spop_min>
min(const SpOp<T1, spop_min>& in, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return SpOp< SpOp<T1, spop_min>, spop_min>(in, dim, 0);
  }



//! @}
