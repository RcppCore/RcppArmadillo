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


//! \addtogroup arma_static_assert
//! @{



//! Classes for primitive compile time assertions (until the next version of C++)
template<bool> struct arma_static_assert;
template<>     struct arma_static_assert<true> {};


template<bool val>
struct arma_type_check
  {
  arma_inline static void apply()
    {
    arma_static_assert<!val> ERROR___INCORRECT_TYPE;
    ERROR___INCORRECT_TYPE = ERROR___INCORRECT_TYPE;
    }
  };



//! @}
