// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_trimat
//! @{



template<typename T1>
inline
void
op_trimat::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trimat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( ((A.n_elem > 0) && (A.is_square() == false)), "trimatu()/trimatl(): given matrix must be square" );
  
  const u32 N = A.n_rows;
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    if(in.aux_u32_a == 0) 
      {
      // upper triangular: copy the diagonal and the elements above the diagonal
      for(u32 i=0; i<N; ++i)
        {
        const eT* A_data   = A.colptr(i);
              eT* out_data = out.colptr(i);
        
        arrayops::copy( out_data, A_data, i+1 );
        }
      }
    else
      {
      // lower triangular: copy the diagonal and the elements below the diagonal
      for(u32 i=0; i<N; ++i)
        {
        const eT* A_data   = A.colptr(i);
              eT* out_data = out.colptr(i);
        
        arrayops::copy( &out_data[i], &A_data[i], N-i );
        }
      }
    }
  
  
  if(in.aux_u32_a == 0) 
    {
    // upper triangular: set all elements below the diagonal to zero
    
    for(u32 i=0; i<N; ++i)
      {
      eT* data = out.colptr(i);
      
      arrayops::inplace_set( &data[i+1], eT(0), (N-(i+1)) );
      }
    }
  else
    {
    // lower triangular: set all elements above the diagonal to zero
    
    for(u32 i=1; i<N; ++i)
      {
      eT* data = out.colptr(i);
      
      arrayops::inplace_set( data, eT(0), i );
      }
    }
  }



//! @}
