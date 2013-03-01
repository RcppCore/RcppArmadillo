// Copyright (C) 2009-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup operator_relational
//! @{


// <  : lt
// >  : gt
// <= : lteq
// >= : gteq
// == : eq
// != : noteq



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_complex<typename T1::elem_type>::value == false) && (is_complex<typename T2::elem_type>::value == false)),
  const mtGlue<uword, T1, T2, glue_rel_lt>
  >::result
operator<
(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_lt>( X, Y );
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_complex<typename T1::elem_type>::value == false) && (is_complex<typename T2::elem_type>::value == false)),
  const mtGlue<uword, T1, T2, glue_rel_gt>
  >::result
operator>
(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_gt>( X, Y );
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_complex<typename T1::elem_type>::value == false) && (is_complex<typename T2::elem_type>::value == false)),
  const mtGlue<uword, T1, T2, glue_rel_lteq>
  >::result
operator<=
(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_lteq>( X, Y );
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_complex<typename T1::elem_type>::value == false) && (is_complex<typename T2::elem_type>::value == false)),
  const mtGlue<uword, T1, T2, glue_rel_gteq>
  >::result
operator>=
(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_gteq>( X, Y );
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value),
  const mtGlue<uword, T1, T2, glue_rel_eq>
  >::result
operator==
(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_eq>( X, Y );
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value),
  const mtGlue<uword, T1, T2, glue_rel_noteq>
  >::result
operator!=
(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_noteq>( X, Y );
  }



//
//
//



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (is_complex<typename T1::elem_type>::value == false)),
  const mtOp<uword, T1, op_rel_lt_pre>
  >::result
operator<
(const typename T1::elem_type val, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_lt_pre>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (is_complex<typename T1::elem_type>::value == false)),
  const mtOp<uword, T1, op_rel_lt_post>
  >::result
operator<
(const T1& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_lt_post>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (is_complex<typename T1::elem_type>::value == false)),
  const mtOp<uword, T1, op_rel_gt_pre>
  >::result
operator>
(const typename T1::elem_type val, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_gt_pre>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (is_complex<typename T1::elem_type>::value == false)),
  const mtOp<uword, T1, op_rel_gt_post>
  >::result
operator>
(const T1& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_gt_post>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (is_complex<typename T1::elem_type>::value == false)),
  const mtOp<uword, T1, op_rel_lteq_pre>
  >::result
operator<=
(const typename T1::elem_type val, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_lteq_pre>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (is_complex<typename T1::elem_type>::value == false)),
  const mtOp<uword, T1, op_rel_lteq_post>
  >::result
operator<=
(const T1& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_lteq_post>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (is_complex<typename T1::elem_type>::value == false)),
  const mtOp<uword, T1, op_rel_gteq_pre>
  >::result
operator>=
(const typename T1::elem_type val, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_gteq_pre>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && (is_complex<typename T1::elem_type>::value == false)),
  const mtOp<uword, T1, op_rel_gteq_post>
  >::result
operator>=
(const T1& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_gteq_post>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const mtOp<uword, T1, op_rel_eq>
  >::result
operator==
(const typename T1::elem_type val, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_eq>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const mtOp<uword, T1, op_rel_eq>
  >::result
operator==
(const T1& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_eq>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const mtOp<uword, T1, op_rel_noteq>
  >::result
operator!=
(const typename T1::elem_type val, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_noteq>(X, val);
  }



template<typename T1>
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const mtOp<uword, T1, op_rel_noteq>
  >::result
operator!=
(const T1& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_noteq>(X, val);
  }



//! @}
