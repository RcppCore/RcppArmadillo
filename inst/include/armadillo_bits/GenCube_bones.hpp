// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup GenCube
//! @{


//! support class for generator functions (eg. zeros, randu, randn, ...)
template<typename eT, typename gen_type>
class GenCube : public BaseCube<eT, GenCube<eT, gen_type> >
  {
  public:
  
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool prefer_at_accessor = false;
  
  arma_aligned const uword n_rows;
  arma_aligned const uword n_cols;
  arma_aligned const uword n_slices;
  
  arma_inline  GenCube(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices);
  arma_inline ~GenCube();
  
  arma_inline static eT generate();
  
  arma_inline eT operator[] (const uword i)                                       const;
  arma_inline eT at         (const uword row, const uword col, const uword slice) const;
  
  inline void apply              (Cube<eT>& out) const;
  inline void apply_inplace_plus (Cube<eT>& out) const;
  inline void apply_inplace_minus(Cube<eT>& out) const;
  inline void apply_inplace_schur(Cube<eT>& out) const;
  inline void apply_inplace_div  (Cube<eT>& out) const;
  };



//! @}
