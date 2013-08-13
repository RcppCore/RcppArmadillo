// Copyright (C) 2008-2012 Conrad Sanderson
// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup Base
//! @{



template<typename elem_type, typename derived>
arma_inline
const derived&
Base<elem_type,derived>::get_ref() const
  {
  return static_cast<const derived&>(*this);
  }



template<typename elem_type, typename derived>
arma_inline
const Op<derived,op_htrans>
Base<elem_type,derived>::t() const
  {
  return Op<derived,op_htrans>( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
arma_inline
const Op<derived,op_htrans>
Base<elem_type,derived>::ht() const
  {
  return Op<derived,op_htrans>( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
arma_inline
const Op<derived,op_strans>
Base<elem_type,derived>::st() const
  {
  return Op<derived,op_strans>( (*this).get_ref() );
  }




template<typename elem_type, typename derived>
inline
void
Base<elem_type,derived>::print(const std::string extra_text) const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
  const quasi_unwrap< typename Proxy<derived>::stored_type > tmp(P.Q);
  
  tmp.M.impl_print(extra_text);
  }



template<typename elem_type, typename derived>
inline
void
Base<elem_type,derived>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
  const quasi_unwrap< typename Proxy<derived>::stored_type > tmp(P.Q);
  
  tmp.M.impl_print(user_stream, extra_text);
  }
  


template<typename elem_type, typename derived>
inline
void
Base<elem_type,derived>::raw_print(const std::string extra_text) const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
  const quasi_unwrap< typename Proxy<derived>::stored_type > tmp(P.Q);
  
  tmp.M.impl_raw_print(extra_text);
  }



template<typename elem_type, typename derived>
inline
void
Base<elem_type,derived>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
  const quasi_unwrap< typename Proxy<derived>::stored_type > tmp(P.Q);
  
  tmp.M.impl_raw_print(user_stream, extra_text);
  }



//
// extra functions defined in Base_blas_elem_type

template<typename derived>
arma_inline
const Op<derived,op_inv>
Base_blas_elem_type<derived>::i(const bool slow) const
  {
  return Op<derived,op_inv>( static_cast<const derived&>(*this), ((slow == false) ? 0 : 1), 0 );
  }



//
// extra functions defined in Base_eval_Mat

template<typename elem_type, typename derived>
arma_inline
const derived&
Base_eval_Mat<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return static_cast<const derived&>(*this);
  }



//
// extra functions defined in Base_eval_expr

template<typename elem_type, typename derived>
arma_inline
Mat<elem_type>
Base_eval_expr<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return Mat<elem_type>( static_cast<const derived&>(*this) );
  }



//! @}
