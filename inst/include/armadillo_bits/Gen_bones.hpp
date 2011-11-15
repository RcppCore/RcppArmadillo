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


//! \addtogroup Gen
//! @{


//! support class for generator functions (eg. zeros, randu, randn, ...)
template<typename eT, typename gen_type>
class Gen : public Base<eT, Gen<eT, gen_type> >
  {
  public:
  
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool prefer_at_accessor = (is_same_type<gen_type, gen_ones_diag>::value) ? true : false;
  
  arma_aligned const uword n_rows;
  arma_aligned const uword n_cols;
  
  arma_inline  Gen(const uword in_n_rows, const uword in_n_cols);
  arma_inline ~Gen();
  
  arma_inline static eT generate();
  
  arma_inline eT operator[] (const uword i)                    const;
  arma_inline eT at         (const uword row, const uword col) const;
  
  inline void apply              (Mat<eT>& out) const;
  inline void apply_inplace_plus (Mat<eT>& out) const;
  inline void apply_inplace_minus(Mat<eT>& out) const;
  inline void apply_inplace_schur(Mat<eT>& out) const;
  inline void apply_inplace_div  (Mat<eT>& out) const;
  };



//! @}
