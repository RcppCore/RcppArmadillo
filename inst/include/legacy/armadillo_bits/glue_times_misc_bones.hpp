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


//! \addtogroup glue_times_misc
//! @{



class dense_sparse_helper
  {
  public:
  
  template<typename eT>
  arma_inline static typename  arma_not_cx<eT>::result dot(const eT* A_mem, const SpMat<eT>& B, const uword col);
  
  template<typename eT>
  arma_inline static typename arma_cx_only<eT>::result dot(const eT* A_mem, const SpMat<eT>& B, const uword col);
  };



class glue_times_dense_sparse
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = T1::is_row;
    static constexpr bool is_col  = T2::is_col;
    static constexpr bool is_xvec = false;
    };
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const SpToDGlue<T1,T2,glue_times_dense_sparse>& expr);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(Mat<typename T1::elem_type>& out, const T1& x, const T2& y);
  
  template<typename T1, typename T2>
  inline static void apply_mixed(Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >& out, const T1& X, const T2& Y);
  };



class glue_times_sparse_dense
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = T1::is_row;
    static constexpr bool is_col  = T2::is_col;
    static constexpr bool is_xvec = false;
    };
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const SpToDGlue<T1,T2,glue_times_sparse_dense>& expr);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(Mat<typename T1::elem_type>& out, const T1& x, const T2& y);
  
  template<typename T1, typename T2>
  inline static void apply_noalias_trans(Mat<typename T1::elem_type>& out, const T1& x, const T2& y);
  
  template<typename T1, typename T2>
  inline static void apply_mixed(Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >& out, const T1& X, const T2& Y);
  };



//! @}
