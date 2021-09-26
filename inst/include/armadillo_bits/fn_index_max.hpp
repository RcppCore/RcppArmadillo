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


//! \addtogroup fn_index_max
//! @{


template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_arma_type<T1>::value && resolves_to_vector<T1>::yes, uword>::result
index_max(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return X.index_max();
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value && resolves_to_vector<T1>::no, const mtOp<uword, T1, op_index_max> >::result
index_max(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_index_max>(X, 0, 0);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const mtOp<uword, T1, op_index_max> >::result
index_max(const T1& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_index_max>(X, dim, 0);
  }



template<typename T1>
arma_warn_unused
arma_inline
const mtOpCube<uword, T1, op_index_max>
index_max
  (
  const BaseCube<typename T1::elem_type, T1>& X,
  const uword dim = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_index_max>(X.get_ref(), dim, 0, 0);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::yes,
  typename T1::elem_type
  >::result
index_max(const T1& x)
  {
  arma_extra_debug_sigprint();
  
  return x.index_max();
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::no,
  Mat<uword>
  >::result
index_max(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<uword> out;
  
  op_index_max::apply(out, X, 0);
  
  return out;
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value,
  Mat<uword>
  >::result
index_max(const T1& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  Mat<uword> out;
  
  op_index_max::apply(out, X, dim);
  
  return out;
  }



arma_warn_unused
inline
uword
index_max(const SizeMat& s)
  {
  return (s.n_rows >= s.n_cols) ? uword(0) : uword(1);
  }



arma_warn_unused
inline
uword
index_max(const SizeCube& s)
  {
  const uword tmp_val   = (s.n_rows >= s.n_cols) ? s.n_rows : s.n_cols;
  const uword tmp_index = (s.n_rows >= s.n_cols) ? uword(0) : uword(1);
  
  return (tmp_val >= s.n_slices) ? tmp_index : uword(2);
  }



//! @}
