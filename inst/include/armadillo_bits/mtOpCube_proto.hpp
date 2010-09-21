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


//! \addtogroup mtOpCube
//! @{



template<typename out_eT, typename T1, typename op_type>
class mtOpCube : public BaseCube<out_eT, mtOpCube<out_eT, T1, op_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;

  typedef typename T1::elem_type                in_eT;

  inline explicit mtOpCube(const T1& in_m);
  inline          mtOpCube(const T1& in_m, const in_eT in_aux);
  inline          mtOpCube(const T1& in_m, const u32   in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  inline          mtOpCube(const T1& in_m, const in_eT in_aux,       const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  
  inline          mtOpCube(const char junk, const T1& in_m, const out_eT in_aux);
  
  inline         ~mtOpCube();
    
  
  arma_aligned const T1&    m;           //!< storage of reference to the operand (e.g. a matrix)
  arma_aligned const in_eT  aux;         //!< storage of auxiliary data, using the element type as used by T1
  arma_aligned const out_eT aux_out_eT;  //!< storage of auxiliary data, using the element type as specified by the out_eT template parameter
  arma_aligned const u32    aux_u32_a;   //!< storage of auxiliary data, u32 format
  arma_aligned const u32    aux_u32_b;   //!< storage of auxiliary data, u32 format
  arma_aligned const u32    aux_u32_c;   //!< storage of auxiliary data, u32 format
  
  };



//! @}
