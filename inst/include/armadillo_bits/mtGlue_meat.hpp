// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
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
  , aux_u32(aux_u32)
  {
  arma_extra_debug_sigprint();
  }



template<typename out_eT, typename T1, typename T2, typename glue_type>
inline
mtGlue<out_eT,T1,T2,glue_type>::mtGlue(const T1& in_A, const T2& in_B, const u32 in_aux_u32)
  : A(in_A)
  , B(in_B)
  , aux_u32(in_aux_u32)
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
