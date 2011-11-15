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


//! \addtogroup eGlueCube
//! @{


template<typename T1, typename T2, typename eglue_type>
class eGlueCube : public BaseCube<typename T1::elem_type, eGlueCube<T1, T2, eglue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool prefer_at_accessor = (ProxyCube<T1>::prefer_at_accessor || ProxyCube<T2>::prefer_at_accessor);
  static const bool has_subview        = (ProxyCube<T1>::has_subview        || ProxyCube<T2>::has_subview       );
  
  arma_aligned const ProxyCube<T1> P1;
  arma_aligned const ProxyCube<T2> P2;
  
  arma_inline ~eGlueCube();
  arma_inline  eGlueCube(const T1& in_A, const T2& in_B);
  
  arma_inline uword get_n_rows()       const;
  arma_inline uword get_n_cols()       const;
  arma_inline uword get_n_elem_slice() const;
  arma_inline uword get_n_slices()     const;
  arma_inline uword get_n_elem()       const;
  
  arma_inline elem_type operator[] (const uword i)                                   const;
  arma_inline elem_type at         (const uword row, const uword col, const uword slice) const;
  };



//! @}
