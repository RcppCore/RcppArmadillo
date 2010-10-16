// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup operator_ostream
//! @{



template<typename eT, typename T1>
inline
std::ostream&
operator<< (std::ostream& o, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(X.get_ref());
  
  arma_ostream::print(o, tmp.M, true);
  
  return o;
  }



template<typename T1>
inline
std::ostream&
operator<< (std::ostream& o, const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(X.get_ref());
  
  arma_ostream::print(o, tmp.M, true);
  
  return o;
  }



//! Print the contents of a field to the specified stream.
template<typename T1>
inline
std::ostream&
operator<< (std::ostream& o, const field<T1>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_ostream::print(o, X);
  
  return o;
  }



//! Print the contents of a subfield to the specified stream
template<typename T1>
inline
std::ostream&
operator<< (std::ostream& o, const subview_field<T1>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_ostream::print(o, X);

  return o;
  }



//! @}
