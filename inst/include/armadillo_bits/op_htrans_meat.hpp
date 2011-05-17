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


//! \addtogroup op_htrans
//! @{



template<typename eT>
arma_inline
void
op_htrans::apply_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_noalias(out, A);
  }



template<typename eT>
inline
void
op_htrans::apply_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const u32 A_n_rows = A.n_rows;
  const u32 A_n_cols = A.n_cols;
  
  out.set_size(A_n_cols, A_n_rows);
  
  for(u32 in_row = 0; in_row < A_n_rows; ++in_row)
    {
    const u32 out_col = in_row;
  
    for(u32 in_col = 0; in_col < A_n_cols; ++in_col)
      {
      const u32 out_row = in_col;
      out.at(out_row, out_col) = std::conj( A.at(in_row, in_col) );
      }
    }
  
  }



template<typename eT>
arma_inline
void
op_htrans::apply(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply(out, A);
  }



template<typename eT>
inline
void
op_htrans::apply(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(&out != &A)
    {
    op_htrans::apply_noalias(out, A);
    }
  else
    {
    if(out.n_rows == out.n_cols)
      {
      arma_extra_debug_print("doing in-place hermitian transpose of a square matrix");
      
      const u32 n_rows = out.n_rows;
      const u32 n_cols = out.n_cols;
      
      for(u32 col=0; col<n_cols; ++col)
        {
        eT* coldata = out.colptr(col);
        
        out.at(col,col) = std::conj( out.at(col,col) );
        
        for(u32 row=(col+1); row<n_rows; ++row)
          {
          const eT val1 = std::conj(coldata[row]);
          const eT val2 = std::conj(out.at(col,row));
          
          out.at(col,row) = val1;
          coldata[row]    = val2;
          }
        }
      }
    else
      {
      const Mat<eT> A_copy = A;
      op_htrans::apply_noalias(out, A_copy);
      }
    }
  
  }



template<typename T1>
inline
void
op_htrans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  op_htrans::apply(out, A);
  }



template<typename T1>
inline
void
op_htrans::apply(Mat<typename T1::elem_type>& out, const Op< Op<T1, op_trimat>, op_htrans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m.m);
  const Mat<eT>& A = tmp.M;
  
  const bool upper = in.m.aux_u32_a;
  
  op_trimat::apply_htrans(out, A, upper);
  }



template<typename T1>
inline
void
op_htrans2::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  op_htrans::apply(out, tmp.M);
  
  arrayops::inplace_mul( out.memptr(), in.aux, out.n_elem );
  }



//! @}
