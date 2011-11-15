// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// Copyright (C) 2011 Alcatel Lucent
// Copyright (C) 2011 Gerhard Schreiber
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup glue_toeplitz
//! @{



template<typename T1, typename T2>
inline
void
glue_toeplitz::apply(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_toeplitz>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if( ((void*)(&in.A)) == ((void*)(&in.B)) )
    {
    arma_extra_debug_print("glue_toeplitz::apply(): one argument version");
    
    const unwrap_check<T1>  tmp(in.A, out);
    const Mat<eT>& A      = tmp.M;
    
    arma_debug_check( (A.is_vec() == false), "toeplitz(): input argument must be a vector" );
    
    const uword N     = A.n_elem;
    const eT* A_mem = A.memptr();
    
    out.set_size(N,N);
    
    for(uword col=0; col<N; ++col)
      {
      eT* col_mem = out.colptr(col);
      
      uword i;
      
      i = col;
      for(uword row=0; row<col; ++row, --i)
        {
        col_mem[row] = A_mem[i];
        }
      
      i = 0;
      for(uword row=col; row<N; ++row, ++i)
        {
        col_mem[row] = A_mem[i];
        }      
      }
    }
  else
    {
    arma_extra_debug_print("glue_toeplitz::apply(): two argument version");
    
    const unwrap_check<T1> tmp1(in.A, out);
    const unwrap_check<T2> tmp2(in.B, out);
    
    const Mat<eT>& A = tmp1.M;
    const Mat<eT>& B = tmp2.M;
    
    arma_debug_check( ( (A.is_vec() == false) || (B.is_vec() == false) ), "toeplitz(): input arguments must be vectors" );
    
    const uword A_N = A.n_elem;
    const uword B_N = B.n_elem;
    
    const eT* A_mem = A.memptr();
    const eT* B_mem = B.memptr();
    
    out.set_size(A_N, B_N);
    
    if( out.is_empty() )
      {
      return;
      }
    
    for(uword col=0; col<B_N; ++col)
      {
      eT* col_mem = out.colptr(col);
      
      uword i = 0;
      for(uword row=col; row<A_N; ++row, ++i)
        {
        col_mem[row] = A_mem[i];
        }
      }
    
    for(uword row=0; row<A_N; ++row)
      {
      uword i = 1;
      for(uword col=(row+1); col<B_N; ++col, ++i)
        {
        out.at(row,col) = B_mem[i];
        }
      }
    
    }
  
  
  }



template<typename T1, typename T2>
inline
void
glue_toeplitz_circ::apply(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_toeplitz_circ>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if( ((void*)(&in.A)) == ((void*)(&in.B)) )
    {
    arma_extra_debug_print("glue_toeplitz_circ::apply(): one argument version");
    
    const unwrap_check<T1>  tmp(in.A, out);
    const Mat<eT>& A      = tmp.M;
    
    arma_debug_check( (A.is_vec() == false), "toeplitz(): input argument must be a vector" );
    
    const uword N     = A.n_elem;
    const eT* A_mem = A.memptr();
    
    out.set_size(N,N);
    
    
    if(A.is_colvec())
      {
      // A is interpreted as colvec
      
      for(uword col=0; col<N; ++col)
        {
        eT* col_mem = out.colptr(col);
        
        
        uword i = col;
        
        for(uword row=0; row<col; ++row, --i)
          {
          col_mem[row] = A_mem[N-i];
          }
        
        
        i = 0;
        
        for(uword row=col; row<N; ++row, ++i)
          {
          col_mem[row] = A_mem[i];
          }      
        }
      
      }
    else
      {
      // A is interpreted as rowvec - probably not the computationally most efficient way to do this ;-)
      
      for(uword row=0; row<N; ++row)
        {
        uword i = row;
         
        for(uword col=0; col<row; ++col, --i)
          {
          out.at(row,col) = A_mem[N-i];
          }
         
        i = 0;
         
        for(uword col=row; col<N; ++col, ++i)
          {
          out.at(row,col) = A_mem[i];
          }
        }
      }
    }
  }



//! @}
