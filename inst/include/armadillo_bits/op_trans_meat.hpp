// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_trans
//! @{



//! Immediate transpose of a dense matrix
template<typename eT>
inline
void
op_trans::apply_noalias(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  const u32 A_n_cols = A.n_cols;
  const u32 A_n_rows = A.n_rows;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (A_n_cols == 1) || (A_n_rows == 1) )
    {
    syslib::copy_elem( out.memptr(), A.mem, A.n_elem );
    }
  else
    {
    for(u32 in_row = 0; in_row<A_n_rows; ++in_row)
      {
      const u32 out_col = in_row;
    
      for(u32 in_col = 0; in_col<A_n_cols; ++in_col)
        {
        const u32 out_row = in_col;
        out.at(out_row, out_col) = A.at(in_row, in_col);
        }
      }
    }
  
  }



template<typename eT>
inline
void
op_trans::apply(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  if(&out != &A)
    {
    op_trans::apply_noalias(out, A);
    }
  else
    {
    if(out.n_rows == out.n_cols)
      {
      arma_extra_debug_print("op_trans::apply(): doing in-place transpose of a square matrix");
      
      const u32 n_rows = out.n_rows;
      const u32 n_cols = out.n_cols;
      
      for(u32 col=0; col<n_cols; ++col)
        {
        eT* coldata = out.colptr(col);
        
        for(u32 row=(col+1); row<n_rows; ++row)
          {
          std::swap( out.at(col,row), coldata[row] );
          }
        }
      }
    else
      {
      const Mat<eT> A_copy = A;
      op_trans::apply_noalias(out, A_copy);
      }
    }
  }



template<typename T1>
inline
void
op_trans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  op_trans::apply(out, A);
  }



// inline void op_trans::apply_inplace(mat &X)
//   {
//   arma_extra_debug_sigprint();
//   
//   if((X.n_rows == 1) || (X.n_cols == 1))
//     {
//     const u32 old_n_rows = X.n_rows;
//     access::rw(X.n_rows) = X.n_cols;
//     access::rw(X.n_cols) = old_n_rows;
//     }
//   else
//   if(X.n_rows == X.n_cols)
//     {
//     for(u32 col=0; col < X.n_cols; ++col)
//       {
//       double* X_coldata = X.colptr(col);
//       
//       for(u32 row=(col+1); row < X.n_rows; ++row)
//         {
//         std::swap( A.at(col,row), A_coldata[row] );
//         }
//       }
//     }
//   else
//     {
//     mat tmp = trans(X);
//     
//     if(X.mem != X.mem_local)
//       {
//       double* old_mem = X.memptr();
//       access::rw(X.mem) = tmp.memptr();
//       access::rw(tmp.mem) = old_mem;
//       }
//     else
//       {
//       X = tmp;
//       }
//     }
//   
//   }




template<typename T1>
inline
void
op_trans2::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trans2>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  op_trans::apply(out, tmp.M);
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = out.n_elem;
  const eT  k       = in.aux;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] *= k;
    }
  
  }



//! @}
