// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_toeplitz
//! @{



template<typename T1>
inline
Glue<T1, T1, glue_toeplitz>
toeplitz(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T1, glue_toeplitz>( X.get_ref(), X.get_ref() );
  }



template<typename T1, typename T2>
inline
Glue<T1, T2, glue_toeplitz>
toeplitz(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_toeplitz>( X.get_ref(), Y.get_ref() );
  }



//! @}
