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


//! \addtogroup op_htrans
//! @{



//! Immediate transpose of a complex matrix
template<typename T>
inline
void
op_htrans::apply_noalias(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(A.n_cols, A.n_rows);
  
  for(u32 in_row = 0; in_row<A.n_rows; ++in_row)
    {
    const u32 out_col = in_row;
  
    for(u32 in_col = 0; in_col<A.n_cols; ++in_col)
      {
      const u32 out_row = in_col;
      out.at(out_row, out_col) = std::conj( A.at(in_row, in_col) );
      }
    }
  
  }



//! Immediate transpose of a complex matrix
template<typename T>
inline
void
op_htrans::apply(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;

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



template<typename T, typename T1>
inline
void
op_htrans::apply(Mat< std::complex<T> >& out, const Op<T1,op_htrans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  isnt_same_type<eT,typename T1::elem_type>::check();
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  op_htrans::apply(out, A);
  }



//! @}
