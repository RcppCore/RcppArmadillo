// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_repmat
//! @{



//! \brief
//! implementation of the 'repeat matrix' operation, used for constructing matrices
template<typename T1>
inline
void
op_repmat::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_repmat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& X     = tmp.M;
  
  arma_debug_check( (X.n_elem == 0), "op_repmat::apply(): given object has no elements" );

  const u32 copies_per_row = in.aux_u32_a;
  const u32 copies_per_col = in.aux_u32_b;
  
  out.set_size(X.n_rows * copies_per_row, X.n_cols * copies_per_col);
  
  for(u32 col = 0; col < out.n_cols; col += X.n_cols)
    {
    for(u32 row = 0; row < out.n_rows; row += X.n_rows)
      {
      out.submat(row, col, row+X.n_rows-1, col+X.n_cols-1) = X;
      }
    }
  }



//! @}
