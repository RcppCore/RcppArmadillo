// Copyright (C) 2009-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_mean
//! @{



template<typename T1>
arma_inline
const Op<T1, op_mean>
mean
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
  
  return Op<T1, op_mean>(X.get_ref(), dim, 0);
  }



template<typename T1>
arma_inline
const Op<T1, op_mean>
mean
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword dim,
  const typename enable_if<resolves_to_vector<T1>::value == true>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_mean>(X.get_ref(), dim, 0);
  }



//! \brief
//! Immediate 'find mean value' operation,
//! invoked, for example, by: mean(mean(A))
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
mean(const Op<T1, op_mean>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("mean(): two consecutive mean() calls detected");
  
  typedef typename T1::elem_type eT;
  
  return op_mean::mean_all(in.m);
  }



template<typename T1>
arma_inline
const Op< Op<T1, op_mean>, op_mean>
mean(const Op<T1, op_mean>& in, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op< Op<T1, op_mean>, op_mean>(in, dim, 0);
  }



template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
mean
  (
  const Base<typename T1::elem_type,T1>& X,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename enable_if<resolves_to_vector<T1>::value == true>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return op_mean::mean_all(X);
  }



template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
mean
  (
  const T1& X,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename enable_if<is_basevec<T1>::value == true>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return op_mean::mean_all(X);
  }



template<typename T>
arma_inline
arma_warn_unused
const typename arma_scalar_only<T>::result &
mean(const T& x)
  {
  return x;
  }



//! @}
