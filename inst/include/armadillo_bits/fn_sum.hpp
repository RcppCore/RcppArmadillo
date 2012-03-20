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


//! \addtogroup fn_sum
//! @{


//! \brief
//! Delayed sum of elements of a matrix along a specified dimension (either rows or columns).
//! The result is stored in a dense matrix that has either one column or one row.
//! For dim = 0, find the sum of each column.
//! For dim = 1, find the sum of each row.
//! The default is dim = 0.
//! NOTE: this function works differently than in Matlab/Octave.

template<typename T1>
arma_inline
const Op<T1, op_sum>
sum
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword dim = 0,
  const typename enable_if<resolves_to_vector<T1>::value == false>::result* junk1 = 0,
  const typename enable_if<is_basevec<T1>::value == false>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return Op<T1, op_sum>(X.get_ref(), dim, 0);
  }



template<typename T1>
arma_inline
const Op<T1, op_sum>
sum
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword dim,
  const typename enable_if<resolves_to_vector<T1>::value == true>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_sum>(X.get_ref(), dim, 0);
  }



//! \brief
//! Immediate 'sum all values' operation,
//! invoked, for example, by: sum(sum(A))

template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
sum(const Op<T1, op_sum>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("sum(): two consecutive sum() calls detected");
  
  return accu(in.m);
  }



template<typename T1>
arma_inline
const Op<Op<T1, op_sum>, op_sum>
sum(const Op<T1, op_sum>& in, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op<Op<T1, op_sum>, op_sum>(in, dim, 0);
  }



//! \brief
//! Immediate 'sum all values' operation for expressions which resolve to a vector 
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
sum
  (
  const Base<typename T1::elem_type,T1>& X,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename enable_if<resolves_to_vector<T1>::value == true>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return accu(X.get_ref());
  }



template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
sum
  (
  const T1& X,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename enable_if<is_basevec<T1>::value == true>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return accu(X);
  }



template<typename T>
arma_inline
arma_warn_unused
const typename arma_scalar_only<T>::result &
sum(const T& x)
  {
  return x;
  }



//! @}
