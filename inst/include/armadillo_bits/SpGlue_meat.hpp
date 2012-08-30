// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup SpGlue
//! @{



template<typename T1, typename T2, typename spglue_type>
inline
SpGlue<T1,T2,spglue_type>::SpGlue(const T1& in_A, const T2& in_B)
  : A(in_A)
  , B(in_B)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename spglue_type>
inline
SpGlue<T1,T2,spglue_type>::SpGlue(const T1& in_A, const T2& in_B, const typename T1::elem_type in_aux)
  : A(in_A)
  , B(in_B)
  , aux(in_aux)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename spglue_type>
inline
SpGlue<T1,T2,spglue_type>::~SpGlue()
  {
  arma_extra_debug_sigprint();
  }



//! @}
