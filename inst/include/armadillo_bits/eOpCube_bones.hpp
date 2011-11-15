// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
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
  
  static const bool prefer_at_accessor = ProxyCube<T1>::prefer_at_accessor;
  static const bool has_subview        = ProxyCube<T1>::has_subview;
  
  arma_aligned const ProxyCube<T1> P;
  arma_aligned       elem_type     aux;          //!< storage of auxiliary data, user defined format
  arma_aligned       uword         aux_uword_a;  //!< storage of auxiliary data, uword format
  arma_aligned       uword         aux_uword_b;  //!< storage of auxiliary data, uword format
  arma_aligned       uword         aux_uword_c;  //!< storage of auxiliary data, uword format
  
  inline         ~eOpCube();
  inline explicit eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b, const uword in_aux_uword_c);
  inline          eOpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux, const uword in_aux_uword_a, const uword in_aux_uword_b, const uword in_aux_uword_c);
  
  arma_inline uword get_n_rows()       const;
  arma_inline uword get_n_cols()       const;
  arma_inline uword get_n_elem_slice() const;
  arma_inline uword get_n_slices()     const;
  arma_inline uword get_n_elem()       const;
  
  arma_inline elem_type operator[] (const uword i)                                       const;
  arma_inline elem_type at         (const uword row, const uword col, const uword slice) const;
  };



//! @}
