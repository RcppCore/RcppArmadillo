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


//! \addtogroup mtGlue
//! @{



template<typename out_eT, typename T1, typename T2, typename glue_type>
inline
mtGlue<out_eT,T1,T2,glue_type>::mtGlue(const T1& in_A, const T2& in_B)
  : A(in_A)
  , B(in_B)
  {
  arma_extra_debug_sigprint();
  }



template<typename out_eT, typename T1, typename T2, typename glue_type>
inline
mtGlue<out_eT,T1,T2,glue_type>::mtGlue(const T1& in_A, const T2& in_B, const uword in_aux_uword)
  : A(in_A)
  , B(in_B)
  , aux_uword(in_aux_uword)
  {
  arma_extra_debug_sigprint();
  }



template<typename out_eT, typename T1, typename T2, typename glue_type>
inline
mtGlue<out_eT,T1,T2,glue_type>::~mtGlue()
  {
  arma_extra_debug_sigprint();
  }



//! @}
