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


//! \addtogroup eOpCube
//! @{



template<typename T1, typename eop_type>
class eOpCube : public BaseCube<typename T1::elem_type, eOpCube<T1, eop_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  const ProxyCube<T1> P;
  
  const elem_type aux;        //!< storage of auxiliary data, user defined format
  const u32       aux_u32_a;  //!< storage of auxiliary data, u32 format
  const u32       aux_u32_b;  //!< storage of auxiliary data, u32 format
  const u32       aux_u32_c;  //!< storage of auxiliary data, u32 format
  
  inline         ~eOpCube();
  inline explicit eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  inline          eOpCube(const u32 in_n_rows, const u32 in_n_cols, const u32 in_n_slices);
  
  };



//! @}
