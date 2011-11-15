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


//! \addtogroup Op
//! @{



template<typename T1, typename op_type>
inline
Op<T1, op_type>::Op(const T1& in_m)
  : m(in_m)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename op_type>
inline
Op<T1, op_type>::Op(const T1& in_m, const typename T1::elem_type in_aux)
  : m(in_m)
  , aux(in_aux)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename op_type>
inline
Op<T1, op_type>::Op(const T1& in_m, const typename T1::elem_type in_aux, const uword in_aux_uword_a, const uword in_aux_uword_b)
  : m(in_m)
  , aux(in_aux)
  , aux_uword_a(in_aux_uword_a)
  , aux_uword_b(in_aux_uword_b)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename op_type>
inline
Op<T1, op_type>::Op(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b)
  : m(in_m)
  , aux_uword_a(in_aux_uword_a)
  , aux_uword_b(in_aux_uword_b)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename op_type>
inline
Op<T1, op_type>::Op(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b, const uword in_aux_uword_c, const char)
  : m(in_m)
  , aux_uword_a(in_aux_uword_a)
  , aux_uword_b(in_aux_uword_b)
  , aux_uword_c(in_aux_uword_c)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename op_type>
inline
Op<T1, op_type>::~Op()
  {
  arma_extra_debug_sigprint();
  }



//! @}
