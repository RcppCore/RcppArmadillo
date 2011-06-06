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


//! \addtogroup op_strans
//! @{



//! for tiny square matrices (size <= 4x4)
template<typename eT>
inline
void
op_strans::apply_noalias_tinysq(Mat<eT>& out, const Mat<eT>& A)
  {
  const eT* Am   = A.memptr();
        eT* outm = out.memptr();
  
  switch(A.n_rows)
    {
    case 1:
      {
      outm[0] = Am[0];
      }
      break;
      
    case 2:
      {
      outm[pos<false,0,0>::n2] = Am[pos<true,0,0>::n2];
      outm[pos<false,1,0>::n2] = Am[pos<true,1,0>::n2];
      
      outm[pos<false,0,1>::n2] = Am[pos<true,0,1>::n2];
      outm[pos<false,1,1>::n2] = Am[pos<true,1,1>::n2];
      }
      break;
    
    case 3:
      {
      outm[pos<false,0,0>::n3] = Am[pos<true,0,0>::n3];
      outm[pos<false,1,0>::n3] = Am[pos<true,1,0>::n3];
      outm[pos<false,2,0>::n3] = Am[pos<true,2,0>::n3];
      
      outm[pos<false,0,1>::n3] = Am[pos<true,0,1>::n3];
      outm[pos<false,1,1>::n3] = Am[pos<true,1,1>::n3];
      outm[pos<false,2,1>::n3] = Am[pos<true,2,1>::n3];
      
      outm[pos<false,0,2>::n3] = Am[pos<true,0,2>::n3];
      outm[pos<false,1,2>::n3] = Am[pos<true,1,2>::n3];
      outm[pos<false,2,2>::n3] = Am[pos<true,2,2>::n3];
      }
      break;
    
    case 4:
      {
      outm[pos<false,0,0>::n4] = Am[pos<true,0,0>::n4];
      outm[pos<false,1,0>::n4] = Am[pos<true,1,0>::n4];
      outm[pos<false,2,0>::n4] = Am[pos<true,2,0>::n4];
      outm[pos<false,3,0>::n4] = Am[pos<true,3,0>::n4];
      
      outm[pos<false,0,1>::n4] = Am[pos<true,0,1>::n4];
      outm[pos<false,1,1>::n4] = Am[pos<true,1,1>::n4];
      outm[pos<false,2,1>::n4] = Am[pos<true,2,1>::n4];
      outm[pos<false,3,1>::n4] = Am[pos<true,3,1>::n4];
      
      outm[pos<false,0,2>::n4] = Am[pos<true,0,2>::n4];
      outm[pos<false,1,2>::n4] = Am[pos<true,1,2>::n4];
      outm[pos<false,2,2>::n4] = Am[pos<true,2,2>::n4];
      outm[pos<false,3,2>::n4] = Am[pos<true,3,2>::n4];
      
      outm[pos<false,0,3>::n4] = Am[pos<true,0,3>::n4];
      outm[pos<false,1,3>::n4] = Am[pos<true,1,3>::n4];
      outm[pos<false,2,3>::n4] = Am[pos<true,2,3>::n4];
      outm[pos<false,3,3>::n4] = Am[pos<true,3,3>::n4];
      }
      break;
    
    default:
      ;
    }
  
  }



//! Immediate transpose of a dense matrix
template<typename eT>
inline
void
op_strans::apply_noalias(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  const u32 A_n_cols = A.n_cols;
  const u32 A_n_rows = A.n_rows;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (A_n_cols == 1) || (A_n_rows == 1) )
    {
    arrayops::copy( out.memptr(), A.mem, A.n_elem );
    }
  else
    {
    if( (A_n_rows <= 4) && (A_n_rows == A_n_cols) )
      {
      op_strans::apply_noalias_tinysq(out, A);
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
  }



template<typename eT>
inline
void
op_strans::apply(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  if(&out != &A)
    {
    op_strans::apply_noalias(out, A);
    }
  else
    {
    if(out.n_rows == out.n_cols)
      {
      arma_extra_debug_print("op_strans::apply(): doing in-place transpose of a square matrix");
      
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
      const Mat<eT> AA = A;
      op_strans::apply_noalias(out, AA);
      }
    }
  }



template<typename T1>
inline
void
op_strans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_strans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  op_strans::apply(out, A);
  }



// inline void op_strans::apply_inplace(mat &X)
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




//! @}
