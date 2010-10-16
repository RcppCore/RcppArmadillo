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


//! \addtogroup eOpCube
//! @{



template<typename T1, typename eop_type>
eOpCube<T1, eop_type>::eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m)
  : P         (in_m.get_ref())
  , aux       (aux)
  , aux_u32_a (aux_u32_a)
  , aux_u32_b (aux_u32_b)
  , aux_u32_c (aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename eop_type>
eOpCube<T1, eop_type>::eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const typename T1::elem_type in_aux)
  : P         (in_m.get_ref())
  , aux       (in_aux)
  , aux_u32_a (aux_u32_a)
  , aux_u32_b (aux_u32_b)
  , aux_u32_c (aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename eop_type>
eOpCube<T1, eop_type>::eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b)
  : P         (in_m.get_ref())
  , aux       (aux)
  , aux_u32_a (in_aux_u32_a)
  , aux_u32_b (in_aux_u32_b)
  , aux_u32_c (aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename eop_type>
eOpCube<T1, eop_type>::eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c)
  : P         (in_m.get_ref())
  , aux       (aux)
  , aux_u32_a (in_aux_u32_a)
  , aux_u32_b (in_aux_u32_b)
  , aux_u32_c (in_aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename eop_type>
eOpCube<T1, eop_type>::eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const typename T1::elem_type in_aux, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c)
  : P         (in_m.get_ref())
  , aux       (in_aux)
  , aux_u32_a (in_aux_u32_a)
  , aux_u32_b (in_aux_u32_b)
  , aux_u32_c (in_aux_u32_c)
  {
  arma_extra_debug_sigprint();
  }



//! used by eop_randu, eop_randn, eop_zeros, eop_ones (i.e. element generators);
//! the proxy P is invalid for generators and must not be used.
template<typename T1, typename eop_type>
eOpCube<T1, eop_type>::eOpCube(const u32 in_n_rows, const u32 in_n_cols, const u32 in_n_slices)
  : P         (P)
  , aux       (aux)
  , aux_u32_a (in_n_rows)
  , aux_u32_b (in_n_cols)
  , aux_u32_c (in_n_slices)
  {
  arma_extra_debug_sigprint();
  }
  


template<typename T1, typename eop_type>
eOpCube<T1, eop_type>::~eOpCube()
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename eop_type>
arma_inline
u32
eOpCube<T1, eop_type>::get_n_rows() const
  {
  return (is_generator<eop_type>::value == false) ? P.get_n_rows() : aux_u32_a;
  }
  


template<typename T1, typename eop_type>
arma_inline
u32
eOpCube<T1, eop_type>::get_n_cols() const
  {
  return (is_generator<eop_type>::value == false) ? P.get_n_cols() : aux_u32_b;
  }



template<typename T1, typename eop_type>
arma_inline
u32
eOpCube<T1, eop_type>::get_n_elem_slice() const
  {
  return (is_generator<eop_type>::value == false) ? P.get_n_elem_slice() : (aux_u32_a * aux_u32_b);
  }



template<typename T1, typename eop_type>
arma_inline
u32
eOpCube<T1, eop_type>::get_n_slices() const
  {
  return (is_generator<eop_type>::value == false) ? P.get_n_slices() : aux_u32_c;
  }



template<typename T1, typename eop_type>
arma_inline
u32
eOpCube<T1, eop_type>::get_n_elem() const
  {
  return (is_generator<eop_type>::value == false) ? P.get_n_elem() : (aux_u32_a * aux_u32_b * aux_u32_c);
  }



template<typename T1, typename eop_type>
arma_inline
typename T1::elem_type
eOpCube<T1, eop_type>::operator[] (const u32 i) const
  {
  typedef typename T1::elem_type eT;
  
  if(is_generator<eop_type>::value == true)
    {
    return eop_aux::generate<eT,eop_type>();
    }
  else
    {
    return eop_core<eop_type>::process(P[i], aux);
    }
  }



template<typename T1, typename eop_type>
arma_inline
typename T1::elem_type
eOpCube<T1, eop_type>::at(const u32 row, const u32 col, const u32 slice) const
  {
  typedef typename T1::elem_type eT;
  
  if(is_generator<eop_type>::value == true)
    {
    return eop_aux::generate<eT,eop_type>();
    }
  else
    {
    return eop_core<eop_type>::process(P.at(row, col, slice), aux);
    }
  }



//! @}
