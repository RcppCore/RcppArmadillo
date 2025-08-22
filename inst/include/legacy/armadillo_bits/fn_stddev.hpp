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


//! \addtogroup fn_stddev
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
stddev(const T1& X, const uword norm_type = 0)
  {
  arma_debug_sigprint();
  
  return op_stddev::stddev_vec(X, norm_type);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_type<T1>::value && resolves_to_vector<T1>::no, 
  const mtOp<typename T1::pod_type, T1, op_stddev>
  >::result
stddev(const T1& X, const uword norm_type = 0)
  {
  arma_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_stddev>(X, norm_type, 0);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const mtOp<typename T1::pod_type, T1, op_stddev>
  >::result
stddev(const T1& X, const uword norm_type, const uword dim)
  {
  arma_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_stddev>(X, norm_type, dim);
  }



template<typename T>
arma_warn_unused
inline
typename arma_scalar_only<T>::result
stddev(const T&)
  {
  return T(0);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::yes,
  typename T1::pod_type
  >::result
stddev(const T1& X, const uword norm_type = 0)
  {
  arma_debug_sigprint();

  return op_sp_stddev::stddev_vec(X, norm_type);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::no, 
  const mtSpReduceOp<typename T1::pod_type, T1, op_sp_stddev>
  >::result
stddev(const T1& X, const uword norm_type = 0)
  {
  arma_debug_sigprint();
  
  return mtSpReduceOp<typename T1::pod_type, T1, op_sp_stddev>(X, norm_type, 0);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value,
  const mtSpReduceOp<typename T1::pod_type, T1, op_sp_stddev>
  >::result
stddev(const T1& X, const uword norm_type, const uword dim)
  {
  arma_debug_sigprint();
  
  return mtSpReduceOp<typename T1::pod_type, T1, op_sp_stddev>(X, norm_type, dim);
  }



//! @}
