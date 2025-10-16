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


//! \addtogroup op_mean
//! @{


//! Class for finding mean values of a matrix
struct op_mean
  : public traits_op_xvec
  {
  // dense matrices
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_mean>& in);
  
  template<typename T1>
  inline static void apply(Mat_noalias<typename T1::elem_type>& out, const Op<T1,op_mean>& in);
  
  template<typename eT>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim);
  
  template<typename eT>
  inline static void apply_noalias_promote(Mat<eT>& out, const Mat<eT>& X, const uword dim);
  
  // cubes
  
  template<typename T1>
  inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_mean>& in);
  
  template<typename eT>
  inline static void apply_noalias(Cube<eT>& out, const Cube<eT>& X, const uword dim);
  
  //
  
  template<typename eT>
  inline static eT direct_mean(const eT* X_mem, const uword N);
  
  template<typename eT>
  inline static eT direct_mean_robust(const eT old_mean, const eT* X_mem, const uword N);
  
  template<typename eT>
  inline static eT direct_mean_promote(const eT* X_mem, const uword N);
  
  template<typename eT>
  inline static eT direct_mean_robust_promote(const eT old_mean, const eT* X_mem, const uword N);
  
  //
  
  template<typename T1>
  inline static typename T1::elem_type mean_all(const T1& X);
  
  template<typename T1>
  inline static typename T1::elem_type mean_all(const Op<T1, op_omit>& X);
  
  template<typename eT, typename functor>
  inline static eT mean_all_omit(const eT* X_mem, const uword N, functor is_omitted);
  
  //
  
  template<typename eT>
  arma_inline static eT robust_mean(const eT A, const eT B);
  
  template<typename T>
  arma_inline static std::complex<T> robust_mean(const std::complex<T>& A, const std::complex<T>& B);
  };


//! @}
