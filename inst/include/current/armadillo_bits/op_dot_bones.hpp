// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup op_dot
//! @{

//! \brief
//! dot product operation 

struct op_dot
  : public traits_op_default
  {
  template<typename eT>
  arma_hot inline static
  typename arma_not_cx<eT>::result
  direct_dot_generic(const uword n_elem, const eT* const A, const eT* const B);
  
  // template<typename eT>
  // arma_hot inline static
  // typename arma_not_cx<eT>::result
  // direct_dot_generic_force_optimise(const uword n_elem, const eT* const A, const eT* const B);
  
  template<typename eT>
  arma_hot inline static
  typename arma_cx_only<eT>::result
  direct_dot_generic(const uword n_elem, const eT* const A, const eT* const B);
  
  template<typename eT>
  arma_hot inline static typename arma_blas_real_only<eT>::result
  direct_dot(const uword n_elem, const eT* const A, const eT* const B);
  
  template<typename eT>
  arma_hot inline static typename arma_blas_cx_only<eT>::result
  direct_dot(const uword n_elem, const eT* const A, const eT* const B);
  
  template<typename eT>
  arma_hot inline static typename arma_fp16_real_only<eT>::result
  direct_dot(const uword n_elem, const eT* const A, const eT* const B);
  
  template<typename eT>
  arma_hot inline static typename arma_fp16_cx_only<eT>::result
  direct_dot(const uword n_elem, const eT* const A, const eT* const B);
  
  template<typename eT>
  arma_hot inline static typename arma_integral_only<eT>::result
  direct_dot(const uword n_elem, const eT* const A, const eT* const B);
  
  
  template<typename eT>
  arma_hot inline static eT direct_dot(const uword n_elem, const eT* const A, const eT* const B, const eT* C);
  
  template<typename T1, typename T2>
  arma_hot inline static typename T1::elem_type apply(const T1& X, const T2& Y);
  
  template<typename T1, typename T2>
  arma_hot inline static typename  arma_not_cx<typename T1::elem_type>::result apply_proxy_linear(const Proxy<T1>& PA, const Proxy<T2>& PB);
  
  template<typename T1, typename T2>
  arma_hot inline static typename arma_cx_only<typename T1::elem_type>::result apply_proxy_linear(const Proxy<T1>& PA, const Proxy<T2>& PB);
  };



//! \brief
//! normalised dot product operation 

struct op_norm_dot
  : public traits_op_default
  {
  template<typename T1, typename T2>
  arma_hot inline static typename T1::elem_type apply(const T1& X, const T2& Y);
  };



//! \brief
//! complex conjugate dot product operation

struct op_cdot
  : public traits_op_default
  {
  template<typename eT>
  arma_hot inline static eT direct_cdot_generic(const uword n_elem, const eT* const A, const eT* const B);
  
  template<typename eT>
  arma_hot inline static eT direct_cdot(const uword n_elem, const eT* const A, const eT* const B, const typename arma_blas_real_or_cx_only<eT>::result* junk = 0);

  template<typename eT>
  arma_hot inline static eT direct_cdot(const uword n_elem, const eT* const A, const eT* const B, const typename arma_fp16_real_or_cx_only<eT>::result* junk = 0);

  template<typename T1, typename T2>
  arma_hot inline static typename T1::elem_type apply       (const T1& X, const T2& Y);
  
  template<typename T1, typename T2>
  arma_hot inline static typename T1::elem_type apply_unwrap(const T1& X, const T2& Y);
  
  template<typename T1, typename T2>
  arma_hot inline static typename T1::elem_type apply_proxy (const T1& X, const T2& Y);
  };



struct op_dot_mixed
  : public traits_op_default
  {
  template<typename T1, typename T2>
  arma_hot inline static
  typename promote_type<typename T1::elem_type, typename T2::elem_type>::result
  apply(const T1& A, const T2& B);
  };



//! @}
