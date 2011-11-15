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
    
    const uword N     = A.n_elem;
    const eT* A_mem = A.memptr();
    
    if(&out != &A)
      {
      // no aliasing
      out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)
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
      
      for(uword i=0; i<N; ++i)
        {
        out.at(i,i) = tmp_mem[i];
        }
      }
    }
  else
    {
    // generate a diagonal matrix out of a matrix
    
    arma_debug_check( (A.is_square() == false), "diagmat(): given matrix is not square" );
    
    const uword N = A.n_rows;
    
    if(&out != &A)
      {
      // no aliasing
      
      out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)
        {
        out.at(i,i) = A.at(i,i);
        }
      }
    else
      {
      // aliasing
      
      for(uword i=0; i<N; ++i)
        {
        eT* colptr = out.colptr(i);
        
        // clear above the diagonal
        arrayops::inplace_set(colptr, eT(0), i);
        
        // clear below the diagonal
        arrayops::inplace_set(colptr+(i+1), eT(0), N-1-i);
        }
      }
    }
  }



//! @}
