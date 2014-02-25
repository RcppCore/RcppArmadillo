// Copyright (C) 2012-2014 Conrad Sanderson
// Copyright (C) 2012-2014 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup SpBase
//! @{



template<typename elem_type, typename derived>
arma_inline
const derived&
SpBase<elem_type,derived>::get_ref() const
  {
  return static_cast<const derived&>(*this);
  }



template<typename elem_type, typename derived>
inline
const SpOp<derived, spop_htrans>
SpBase<elem_type,derived>::t() const
  {
  return SpOp<derived,spop_htrans>( (*this).get_ref() );
  }


template<typename elem_type, typename derived>
inline
const SpOp<derived, spop_htrans>
SpBase<elem_type,derived>::ht() const
  {
  return SpOp<derived, spop_htrans>( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
const SpOp<derived, spop_strans>
SpBase<elem_type,derived>::st() const
  {
  return SpOp<derived, spop_strans>( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
void
SpBase<elem_type,derived>::print(const std::string extra_text) const
  {
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  tmp.M.impl_print(extra_text);
  }



template<typename elem_type, typename derived>
inline
void
SpBase<elem_type,derived>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  tmp.M.impl_print(user_stream, extra_text);
  }
  


template<typename elem_type, typename derived>
inline
void
SpBase<elem_type,derived>::raw_print(const std::string extra_text) const
  {
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  tmp.M.impl_raw_print(extra_text);
  }



template<typename elem_type, typename derived>
inline
void
SpBase<elem_type,derived>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  tmp.M.impl_raw_print(user_stream, extra_text);
  }



template<typename elem_type, typename derived>
inline
void
SpBase<elem_type, derived>::print_dense(const std::string extra_text) const
  {
  const unwrap_spmat<derived> tmp( (*this).get_ref() );

  tmp.M.impl_print_dense(extra_text);
  }



template<typename elem_type, typename derived>
inline
void
SpBase<elem_type, derived>::print_dense(std::ostream& user_stream, const std::string extra_text) const
  {
  const unwrap_spmat<derived> tmp( (*this).get_ref() );

  tmp.M.impl_print_dense(user_stream, extra_text);
  }



template<typename elem_type, typename derived>
inline
void
SpBase<elem_type, derived>::raw_print_dense(const std::string extra_text) const
  {
  const unwrap_spmat<derived> tmp( (*this).get_ref() );

  tmp.M.impl_raw_print_dense(extra_text);
  }



template<typename elem_type, typename derived>
inline
void
SpBase<elem_type, derived>::raw_print_dense(std::ostream& user_stream, const std::string extra_text) const
  {
  const unwrap_spmat<derived> tmp( (*this).get_ref() );

  tmp.M.impl_raw_print_dense(user_stream, extra_text);
  }



//
// extra functions defined in SpBase_eval_SpMat

template<typename elem_type, typename derived>
inline
const derived&
SpBase_eval_SpMat<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return static_cast<const derived&>(*this);
  }



//
// extra functions defined in SpBase_eval_expr

template<typename elem_type, typename derived>
inline
SpMat<elem_type>
SpBase_eval_expr<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return SpMat<elem_type>( static_cast<const derived&>(*this) );
  }



//! @}
