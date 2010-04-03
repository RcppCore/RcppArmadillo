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


//! \addtogroup OpCube
//! @{



template<typename T1, typename op_type>
OpCube<T1, op_type>::OpCube(const BaseCube<typename T1::elem_type, T1>& in_m)
  : m(in_m.get_ref())
  , aux(aux)              // don't tell mum
  , aux_u32_a(aux_u32_a)
  , aux_u32_b(aux_u32_b)
  , aux_u32_c(aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename op_type>
OpCube<T1, op_type>::OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const typename T1::elem_type in_aux)
  : m(in_m.get_ref())
  , aux(in_aux)
  , aux_u32_a(aux_u32_a)
  , aux_u32_b(aux_u32_b)
  , aux_u32_c(aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename op_type>
OpCube<T1, op_type>::OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b)
  : m(in_m.get_ref())
  , aux(aux)
  , aux_u32_a(in_aux_u32_a)
  , aux_u32_b(in_aux_u32_b)
  , aux_u32_c(aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename op_type>
OpCube<T1, op_type>::OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c)
  : m(in_m.get_ref())
  , aux(aux)
  , aux_u32_a(in_aux_u32_a)
  , aux_u32_b(in_aux_u32_b)
  , aux_u32_c(in_aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename op_type>
OpCube<T1, op_type>::OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const typename T1::elem_type in_aux, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c)
  : m(in_m.get_ref())
  , aux(in_aux)
  , aux_u32_a(in_aux_u32_a)
  , aux_u32_b(in_aux_u32_b)
  , aux_u32_c(in_aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename op_type>
OpCube<T1, op_type>::OpCube(const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c)
  : m(m)
  , aux(aux)
  , aux_u32_a(in_aux_u32_a)
  , aux_u32_b(in_aux_u32_b)
  , aux_u32_c(in_aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename op_type>
OpCube<T1, op_type>::~OpCube()
  {
  arma_extra_debug_sigprint();
  }



//! @}
