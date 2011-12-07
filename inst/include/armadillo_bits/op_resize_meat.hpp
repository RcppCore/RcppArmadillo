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



//! \addtogroup op_resize
//! @{



template<typename T1>
inline
void
op_resize::apply(Mat<typename T1::elem_type>& actual_out, const Op<T1,op_resize>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword out_n_rows = in.aux_uword_a;
  const uword out_n_cols = in.aux_uword_b;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  Mat<eT> B;
  
  const bool alias = (&actual_out == &A);
  
  Mat<eT>& out = alias ? B : actual_out;
  
  out.set_size(out_n_rows, out_n_cols);
  
  if( (out_n_rows > A_n_rows) || (out_n_cols > A_n_cols) )
    {
    out.zeros();
    }
  
  if(out.n_elem > 0)
    {
    const uword end_row = (std::min)(out_n_rows, A_n_rows) - 1;
    const uword end_col = (std::min)(out_n_cols, A_n_cols) - 1;
    
    out.submat(0, 0, end_row, end_col) = A.submat(0, 0, end_row, end_col);
    }
  
  if(alias)
    {
    actual_out.steal_mem(B);
    }
  
  }



template<typename T1>
inline
void
op_resize::apply(Cube<typename T1::elem_type>& actual_out, const OpCube<T1,op_resize>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword out_n_rows   = in.aux_uword_a;
  const uword out_n_cols   = in.aux_uword_b;
  const uword out_n_slices = in.aux_uword_c;
  
  const unwrap_cube<T1> tmp(in.m);
  const Cube<eT>& A   = tmp.M;
  
  const uword A_n_rows   = A.n_rows;
  const uword A_n_cols   = A.n_cols;
  const uword A_n_slices = A.n_slices;
  
  Cube<eT> B;
  
  const bool alias = (&actual_out == &A);
  
  Cube<eT>& out = alias ? B : actual_out;
  
  out.set_size(out_n_rows, out_n_cols, out_n_slices);
  
  if( (out_n_rows > A_n_rows) || (out_n_cols > A_n_cols) || (out_n_slices > A_n_slices) )
    {
    out.zeros();
    }
  
  if(out.n_elem > 0)
    {
    const uword end_row   = (std::min)(out_n_rows,   A_n_rows)   - 1;
    const uword end_col   = (std::min)(out_n_cols,   A_n_cols)   - 1;
    const uword end_slice = (std::min)(out_n_slices, A_n_slices) - 1;
    
    out.subcube(0, 0, 0, end_row, end_col, end_slice) = A.subcube(0, 0, 0, end_row, end_col, end_slice);
    }
  
  if(alias)
    {
    actual_out.steal_mem(B);
    }
  
  }



//! @}
