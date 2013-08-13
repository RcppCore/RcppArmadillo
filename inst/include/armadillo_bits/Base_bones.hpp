// Copyright (C) 2008-2012 Conrad Sanderson
// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup Base
//! @{



template<typename derived>
struct Base_blas_elem_type
  {
  arma_inline const Op<derived,op_inv> i(const bool slow = false) const;   //!< matrix inverse
  };


template<typename derived>
struct Base_other_elem_type
  {
  };


template<typename derived, bool condition>
struct Base_extra {};

template<typename derived>
struct Base_extra<derived, true>  { typedef Base_blas_elem_type<derived>  result; };

template<typename derived>
struct Base_extra<derived, false> { typedef Base_other_elem_type<derived> result; };



template<typename elem_type, typename derived>
struct Base_eval_Mat
  {
  const derived& eval() const;
  };


template<typename elem_type, typename derived>
struct Base_eval_expr
  {
  Mat<elem_type> eval() const;   //!< force the immediate evaluation of a delayed expression
  };


template<typename elem_type, typename derived, bool condition>
struct Base_eval {};

template<typename elem_type, typename derived>
struct Base_eval<elem_type, derived, true>  { typedef Base_eval_Mat<elem_type, derived>  result; };

template<typename elem_type, typename derived>
struct Base_eval<elem_type, derived, false> { typedef Base_eval_expr<elem_type, derived> result; };



//! Class for static polymorphism, modelled after the "Curiously Recurring Template Pattern" (CRTP).
//! Used for type-safe downcasting in functions that restrict their input(s) to be classes that are
//! derived from Base (e.g. Mat, Op, Glue, diagview, subview).
//! A Base object can be converted to a Mat object by the unwrap class.

template<typename elem_type, typename derived>
struct Base
  : public Base_extra<derived, is_supported_blas_type<elem_type>::value>::result
  , public Base_eval<elem_type, derived, is_Mat<derived>::value>::result
  {
  arma_inline const derived& get_ref() const;
  
  arma_inline const Op<derived,op_htrans>  t() const;  //!< Hermitian transpose
  arma_inline const Op<derived,op_htrans> ht() const;  //!< Hermitian transpose
  arma_inline const Op<derived,op_strans> st() const;  //!< simple transpose
  
  inline void print(const std::string extra_text = "") const;
  inline void print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void raw_print(const std::string extra_text = "") const;
  inline void raw_print(std::ostream& user_stream, const std::string extra_text = "") const;
  };



//! @}
