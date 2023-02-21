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


//! \addtogroup SpBase
//! @{



template<typename elem_type, typename derived>
struct SpBase_eval_SpMat
  {
  arma_warn_unused inline const derived& eval() const;
  };


template<typename elem_type, typename derived>
struct SpBase_eval_expr
  {
  arma_warn_unused inline SpMat<elem_type> eval() const;   //!< force the immediate evaluation of a delayed expression
  };


template<typename elem_type, typename derived, bool condition>
struct SpBase_eval {};

template<typename elem_type, typename derived>
struct SpBase_eval<elem_type, derived, true>  { typedef SpBase_eval_SpMat<elem_type, derived> result; };

template<typename elem_type, typename derived>
struct SpBase_eval<elem_type, derived, false> { typedef SpBase_eval_expr<elem_type, derived>  result; };



template<typename elem_type, typename derived>
struct SpBase
  : public SpBase_eval<elem_type, derived, is_SpMat<derived>::value>::result
  {
  arma_inline const derived& get_ref() const;
  
  arma_inline bool is_alias(const SpMat<elem_type>& X) const;
  
  arma_warn_unused inline const SpOp<derived,spop_htrans>  t() const;  //!< Hermitian transpose
  arma_warn_unused inline const SpOp<derived,spop_htrans> ht() const;  //!< Hermitian transpose
  arma_warn_unused inline const SpOp<derived,spop_strans> st() const;  //!< simple transpose
  
  arma_cold inline void print(                           const std::string extra_text = "") const;
  arma_cold inline void print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_cold inline void raw_print(                           const std::string extra_text = "") const;
  arma_cold inline void raw_print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_cold inline void print_dense(                           const std::string extra_text = "") const;
  arma_cold inline void print_dense(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_cold inline void raw_print_dense(                           const std::string extra_text = "") const;
  arma_cold inline void raw_print_dense(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_cold inline void brief_print(                           const std::string extra_text = "") const;
  arma_cold inline void brief_print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_warn_unused inline elem_type min() const;
  arma_warn_unused inline elem_type max() const;
  
  inline elem_type min(uword& index_of_min_val) const;
  inline elem_type max(uword& index_of_max_val) const;
  
  inline elem_type min(uword& row_of_min_val, uword& col_of_min_val) const;
  inline elem_type max(uword& row_of_max_val, uword& col_of_max_val) const;
  
  arma_warn_unused inline uword index_min() const;
  arma_warn_unused inline uword index_max() const;
  
  arma_warn_unused inline bool is_symmetric() const;
  arma_warn_unused inline bool is_symmetric(const typename get_pod_type<elem_type>::result tol) const;
  
  arma_warn_unused inline bool is_hermitian() const;
  arma_warn_unused inline bool is_hermitian(const typename get_pod_type<elem_type>::result tol) const;
  
  arma_warn_unused inline bool is_zero(const typename get_pod_type<elem_type>::result tol = 0) const;
  
  arma_warn_unused inline bool is_trimatu() const;
  arma_warn_unused inline bool is_trimatl() const;
  arma_warn_unused inline bool is_diagmat() const;
  arma_warn_unused inline bool is_empty()   const;
  arma_warn_unused inline bool is_square()  const;
  arma_warn_unused inline bool is_vec()     const;
  arma_warn_unused inline bool is_colvec()  const;
  arma_warn_unused inline bool is_rowvec()  const;
  arma_warn_unused inline bool is_finite()  const;
  arma_warn_unused inline bool has_inf()    const;
  arma_warn_unused inline bool has_nan()    const;
  
  arma_warn_unused inline const SpOp<derived,spop_vectorise_col> as_col() const;
  arma_warn_unused inline const SpOp<derived,spop_vectorise_row> as_row() const;
  };



//! @}
