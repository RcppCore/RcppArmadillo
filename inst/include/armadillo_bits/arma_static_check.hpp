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


//! \addtogroup arma_static_check
//! @{



template<bool ERROR___INCORRECT_OR_UNSUPPORTED_TYPE>
struct arma_type_check_cxx1998
  {
  arma_inline
  static
  void
  apply()
    {
    static const char
    junk[ ERROR___INCORRECT_OR_UNSUPPORTED_TYPE ? -1 : +1 ];
    }
  };



template<>
struct arma_type_check_cxx1998<false>
  {
  arma_inline
  static
  void
  apply()
    {
    }
  };



#if !defined(ARMA_USE_CXX11)

  #define arma_static_check(condition, message)  static const char message[ (condition) ? -1 : +1 ]
  
  #define arma_type_check(condition)  arma_type_check_cxx1998<condition>::apply()

#else

  #define arma_static_check(condition, message)  static_assert( !(condition), #message )
  
  #define arma_type_check(condition)  static_assert( !(condition), "error: incorrect or unsupported type" )

#endif



//! @}
