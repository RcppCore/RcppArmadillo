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


//! \addtogroup OpCube
//! @{


//! Analog of the Op class, intended for cubes

template<typename T1, typename op_type>
class OpCube : public BaseCube<typename T1::elem_type, OpCube<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  
  inline explicit OpCube(const BaseCube<typename T1::elem_type, T1>& in_m);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux, const uword in_aux_uword_a, const uword in_aux_uword_b, const uword in_aux_uword_c);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b, const uword in_aux_uword_c);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b, const uword in_aux_uword_c, const uword in_aux_uword_d, const char junk);
  inline         ~OpCube();
  
  arma_aligned const T1&       m;            //!< storage of reference to the operand (e.g. a cube)
  arma_aligned       elem_type aux;          //!< storage of auxiliary data, user defined format
  arma_aligned       uword     aux_uword_a;  //!< storage of auxiliary data, uword format
  arma_aligned       uword     aux_uword_b;  //!< storage of auxiliary data, uword format
  arma_aligned       uword     aux_uword_c;  //!< storage of auxiliary data, uword format
  arma_aligned       uword     aux_uword_d;  //!< storage of auxiliary data, uword format
  
  };



//! @}
