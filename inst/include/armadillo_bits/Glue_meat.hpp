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


//! \addtogroup Glue
//! @{



template<typename T1, typename T2, typename glue_type>
inline
Glue<T1,T2,glue_type>::Glue(const T1& in_A, const T2& in_B)
  : A(in_A)
  , B(in_B)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename glue_type>
inline
Glue<T1,T2,glue_type>::Glue(const T1& in_A, const T2& in_B, const uword in_aux_uword)
  : A(in_A)
  , B(in_B)
  , aux_uword(in_aux_uword)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename glue_type>
inline
Glue<T1,T2,glue_type>::~Glue()
  {
  arma_extra_debug_sigprint();
  }



//! @}
