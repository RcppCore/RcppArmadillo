// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup fn_vecnorm
//! @{



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_type<T1>::value && resolves_to_vector<T1>::yes,
  typename T1::pod_type
  >::result
vecnorm
  (
  const T1&   X,
  const uword k = uword(2),
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> P(X);
  
  if(P.get_n_elem() == 0)  { return T(0); }
  
  if(k == uword(1))  { return op_norm::vec_norm_1(P); }
  if(k == uword(2))  { return op_norm::vec_norm_2(P); }
  
  arma_debug_check( (k == 0), "vecnorm(): unsupported vector norm type" );
  
  return op_norm::vec_norm_k(P, int(k));
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_type<T1>::value && resolves_to_vector<T1>::no,
  const mtOp<typename T1::pod_type, T1, op_vecnorm>
  >::result
vecnorm
  (
  const T1&   X,
  const uword k = uword(2),
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  const uword dim = 0;
  
  return mtOp<typename T1::pod_type, T1, op_vecnorm>(X, k, dim);
  }



template<typename T1>
arma_warn_unused
inline
const mtOp<typename T1::pod_type, T1, op_vecnorm>
vecnorm
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword k,
  const uword dim,
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return mtOp<typename T1::pod_type, T1, op_vecnorm>(X.get_ref(), k, dim);
  }



//



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_type<T1>::value && resolves_to_vector<T1>::yes,
  typename T1::pod_type
  >::result
vecnorm
  (
  const T1&   X,
  const char* method,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> P(X);
  
  if(P.get_n_elem() == 0)  { return T(0); }
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )  { return op_norm::vec_norm_max(P); }
  if( (sig == '-')                                 )  { return op_norm::vec_norm_min(P); }
  
  arma_stop_logic_error("vecnorm(): unsupported vector norm type");
  
  return T(0);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_type<T1>::value && resolves_to_vector<T1>::no, 
  const mtOp<typename T1::pod_type, T1, op_vecnorm_ext>
  >::result
vecnorm
  (
  const T1&   X,
  const char* method,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  uword method_id = 0;
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )  { method_id = 1; }
  if( (sig == '-')                                 )  { method_id = 2; }
  
  const uword dim = 0;
  
  return mtOp<typename T1::pod_type, T1, op_vecnorm_ext>(X, method_id, dim);
  }



template<typename T1>
arma_warn_unused
inline
const mtOp<typename T1::pod_type, T1, op_vecnorm_ext>
vecnorm
  (
  const Base<typename T1::elem_type,T1>& X,
  const char* method,
  const uword dim,
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  uword method_id = 0;
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )  { method_id = 1; }
  if( (sig == '-')                                 )  { method_id = 2; }
  
  return mtOp<typename T1::pod_type, T1, op_vecnorm_ext>(X.get_ref(), method_id, dim);
  }



//
// norms for sparse matrices



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::yes,
  typename T1::pod_type
  >::result
vecnorm
  (
  const T1&   X,
  const uword k = uword(2),
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return arma::norm(X, k);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::no,
  const mtSpOp<typename T1::pod_type, T1, spop_vecnorm>
  >::result
vecnorm
  (
  const T1&   X,
  const uword k = uword(2),
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  const uword dim = 0;
  
  return mtSpOp<typename T1::pod_type, T1, spop_vecnorm>(X, k, dim);
  }



template<typename T1>
arma_warn_unused
inline
const mtSpOp<typename T1::pod_type, T1, spop_vecnorm>
vecnorm
  (
  const SpBase<typename T1::elem_type,T1>& X,
  const uword k,
  const uword dim,
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return mtSpOp<typename T1::pod_type, T1, spop_vecnorm>(X.get_ref(), k, dim);
  }



//



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::yes,
  typename T1::pod_type
  >::result
vecnorm
  (
  const T1&   X,
  const char* method,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return arma::norm(X, method);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::no,
  const mtSpOp<typename T1::pod_type, T1, spop_vecnorm_ext>
  >::result
vecnorm
  (
  const T1&   X,
  const char* method,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk2 = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  uword method_id = 0;
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )  { method_id = 1; }
  if( (sig == '-')                                 )  { method_id = 2; }
  
  const uword dim = 0;
  
  return mtSpOp<typename T1::pod_type, T1, spop_vecnorm_ext>(X, method_id, dim);
  }



template<typename T1>
arma_warn_unused
inline
const mtSpOp<typename T1::pod_type, T1, spop_vecnorm_ext>
vecnorm
  (
  const SpBase<typename T1::elem_type,T1>& X,
  const char* method,
  const uword dim,
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  uword method_id = 0;
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )  { method_id = 1; }
  if( (sig == '-')                                 )  { method_id = 2; }
  
  return mtSpOp<typename T1::pod_type, T1, spop_vecnorm_ext>(X.get_ref(), method_id, dim);
  }



//! @}
