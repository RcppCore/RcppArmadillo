// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup BaseCube
//! @{



template<typename elem_type, typename derived>
arma_inline
const derived&
BaseCube<elem_type,derived>::get_ref() const
  {
  return static_cast<const derived&>(*this);
  }



template<typename elem_type, typename derived>
inline
void
BaseCube<elem_type,derived>::print(const std::string extra_text) const
  {
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  tmp.M.impl_print(extra_text);
  }



template<typename elem_type, typename derived>
inline
void
BaseCube<elem_type,derived>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  tmp.M.impl_print(user_stream, extra_text);
  }
  


template<typename elem_type, typename derived>
inline
void
BaseCube<elem_type,derived>::raw_print(const std::string extra_text) const
  {
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  tmp.M.impl_raw_print(extra_text);
  }



template<typename elem_type, typename derived>
inline
void
BaseCube<elem_type,derived>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  tmp.M.impl_raw_print(user_stream, extra_text);
  }
  


//! @}
