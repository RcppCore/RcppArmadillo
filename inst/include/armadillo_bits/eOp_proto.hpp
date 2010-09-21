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


//! \addtogroup eOp
//! @{



template<typename T1, typename eop_type>
class eOp : public Base<typename T1::elem_type, eOp<T1, eop_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef          Proxy<T1>                       proxy_type;
  
  arma_aligned const Proxy<T1> P;
  arma_aligned const elem_type aux;        //!< storage of auxiliary data, user defined format
  arma_aligned const u32       aux_u32_a;  //!< storage of auxiliary data, u32 format
  arma_aligned const u32       aux_u32_b;  //!< storage of auxiliary data, u32 format
  
  inline         ~eOp();
  inline explicit eOp(const Base<typename T1::elem_type, T1>& in_m);
  inline          eOp(const Base<typename T1::elem_type, T1>& in_m, const elem_type in_aux);
  inline          eOp(const Base<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b);
  inline          eOp(const Base<typename T1::elem_type, T1>& in_m, const elem_type in_aux, const u32 in_aux_u32_a, const u32 in_aux_u32_b);
  inline          eOp(const u32 in_n_rows, const u32 in_n_cols);
  
  arma_inline u32 get_n_rows() const;
  arma_inline u32 get_n_cols() const;
  arma_inline u32 get_n_elem() const;
  
  arma_inline elem_type operator[] (const u32 i)                  const;
  arma_inline elem_type at         (const u32 row, const u32 col) const;
  };



//! @}
