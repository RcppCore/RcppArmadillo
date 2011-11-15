// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
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
  
  const uword copies_per_row = in.aux_uword_a;
  const uword copies_per_col = in.aux_uword_b;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  out.set_size(X_n_rows * copies_per_row, X_n_cols * copies_per_col);
  
  const uword out_n_rows = out.n_rows;
  const uword out_n_cols = out.n_cols;
  
  if( (out_n_rows > 0) && (out_n_cols > 0) )
    {
    for(uword col = 0; col < out_n_cols; col += X_n_cols)
      {
      for(uword row = 0; row < out_n_rows; row += X_n_rows)
        {
        out.submat(row, col, row+X_n_rows-1, col+X_n_cols-1) = X;
        }
      }
    }
  }



//! @}
