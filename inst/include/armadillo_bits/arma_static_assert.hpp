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


//! \addtogroup arma_static_assert
//! @{



//! Classes for primitive compile time assertions and checks (until the next version of C++)
template<bool x>
struct arma_static_assert
  {
  static const char
  static_error[  x ? +1 : -1 ];
  };



template<bool x>
struct arma_static_check
  {
  static const char
  static_error[  x ? -1 : +1 ];
  };



template<bool val>
struct arma_type_check
  {
  arma_inline
  static
  void
  apply()
    {
    arma_static_check<val> ERROR___INCORRECT_TYPE;
    ERROR___INCORRECT_TYPE = ERROR___INCORRECT_TYPE;
    }
  };



template<>
struct arma_type_check<false>
  {
  arma_inline
  static
  void
  apply()
    {
    }
  };


//! @}
