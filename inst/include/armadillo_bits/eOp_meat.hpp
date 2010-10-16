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


//! \addtogroup eOp
//! @{



template<typename T1, typename eop_type>
eOp<T1, eop_type>::eOp(const Base<typename T1::elem_type, T1>& in_m)
  : P(in_m.get_ref())
  , aux(aux)
  , aux_u32_a(aux_u32_a)
  , aux_u32_b(aux_u32_b)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename eop_type>
eOp<T1, eop_type>::eOp(const Base<typename T1::elem_type, T1>& in_m, const typename T1::elem_type in_aux)
  : P(in_m.get_ref())
  , aux(in_aux)
  , aux_u32_a(aux_u32_a)
  , aux_u32_b(aux_u32_b)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename eop_type>
eOp<T1, eop_type>::eOp(const Base<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b)
  : P(in_m.get_ref())
  , aux(aux)
  , aux_u32_a(in_aux_u32_a)
  , aux_u32_b(in_aux_u32_b)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename eop_type>
eOp<T1, eop_type>::eOp(const Base<typename T1::elem_type, T1>& in_m, const typename T1::elem_type in_aux, const u32 in_aux_u32_a, const u32 in_aux_u32_b)
  : P(in_m.get_ref())
  , aux(in_aux)
  , aux_u32_a(in_aux_u32_a)
  , aux_u32_b(in_aux_u32_b)
  {
  arma_extra_debug_sigprint();
  }



//! used by eop_randu, eop_randn, eop_zeros, eop_ones, eop_ones_diag (i.e. element generators);
//! the proxy P is invalid for generators and must not be used.
template<typename T1, typename eop_type>
eOp<T1, eop_type>::eOp(const u32 in_n_rows, const u32 in_n_cols)
  : P(P)
  , aux(aux)
  , aux_u32_a(in_n_rows)
  , aux_u32_b(in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename eop_type>
eOp<T1, eop_type>::~eOp()
  {
  arma_extra_debug_sigprint();
  }

  

template<typename T1, typename eop_type>
arma_inline
u32
eOp<T1, eop_type>::get_n_rows() const
  {
  return (is_generator<eop_type>::value == false) ? P.get_n_rows() : aux_u32_a;
  }
  


template<typename T1, typename eop_type>
arma_inline
u32
eOp<T1, eop_type>::get_n_cols() const
  {
  return (is_generator<eop_type>::value == false) ? P.get_n_cols() : aux_u32_b;
  }



template<typename T1, typename eop_type>
arma_inline
u32
eOp<T1, eop_type>::get_n_elem() const
  {
  return (is_generator<eop_type>::value == false) ? P.get_n_elem() : (aux_u32_a * aux_u32_b);
  }



template<typename T1, typename eop_type>
arma_inline
typename T1::elem_type
eOp<T1, eop_type>::operator[] (const u32 i) const
  {
  typedef typename T1::elem_type eT;
  
  if(is_generator<eop_type>::value == true)
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      return ((i % get_n_rows()) == (i / get_n_rows())) ? eT(1) : eT(0);
      }
    else
      {  
      return eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    return eop_core<eop_type>::process(P[i], aux);
    }
  }



template<typename T1, typename eop_type>
arma_inline
typename T1::elem_type
eOp<T1, eop_type>::at(const u32 row, const u32 col) const
  {
  typedef typename T1::elem_type eT;
  
  if(is_generator<eop_type>::value == true)
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      return (row == col) ? eT(1) : eT(0);
      }
    else
      {
      return eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    return eop_core<eop_type>::process(P.at(row, col), aux);
    }
  }



//! @}
