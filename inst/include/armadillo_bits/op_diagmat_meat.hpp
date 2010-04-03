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


//! \addtogroup op_diagmat
//! @{



template<typename T1>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.m);
  const Mat<eT>& A = tmp.M;
  
  if(A.is_vec() == true)
    {
    // generate a diagonal matrix out of a vector
    
    const u32 N     = A.n_elem;
    const eT* A_mem = A.memptr();
    
    if(&out != &A)
      {
      // no aliasing
      out.zeros(N,N);
      
      for(u32 i=0; i<N; ++i)
        {
        out.at(i,i) = A_mem[i];
        }
      }
    else
      {
      // aliasing
      
      const podarray<eT> tmp(A_mem, N);
      
      const eT* tmp_mem = tmp.memptr();
      
      out.zeros(N,N);
      
      for(u32 i=0; i<N; ++i)
        {
        out.at(i,i) = tmp_mem[i];
        }
      }
    }
  else
    {
    // generate a diagonal matrix out of a matrix
    
    arma_debug_check( (A.is_square() == false), "diagmat(): given matrix is not square" );
    
    const u32 N = A.n_rows;
    
    out.set_size(N,N);
    
    for(u32 col=0; col<N; ++col)
      {
      for(u32 row=0;     row<col; ++row) { out.at(row,col) = eT(0); }
      
      out.at(col,col) = A.at(col,col);
      
      for(u32 row=col+1; row<N;   ++row) { out.at(row,col) = eT(0); }
      }
    }
  }



//! @}
