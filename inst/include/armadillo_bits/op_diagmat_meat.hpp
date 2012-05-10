// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
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
  
  const Proxy<T1> P(X.m);
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  const bool P_is_vec = (n_rows == 1) || (n_cols == 1);
  
  
  if(P.is_alias(out) == false)
    {
    if(P_is_vec)    // generate a diagonal matrix out of a vector
      {
      const uword N = (n_rows == 1) ? n_cols : n_rows;
      
      out.zeros(N, N);
      
      if(Proxy<T1>::prefer_at_accessor == false)
        {
        typename Proxy<T1>::ea_type P_ea = P.get_ea();
        
        for(uword i=0; i < N; ++i) { out.at(i,i) = P_ea[i]; }
        }
      else
        {
        if(n_rows == 1)
          {
          for(uword i=0; i < N; ++i) { out.at(i,i) = P.at(0,i); }
          }
        else
          {
          for(uword i=0; i < N; ++i) { out.at(i,i) = P.at(i,0); }
          }
        }
      }
    else   // generate a diagonal matrix out of a matrix
      {
      arma_debug_check( (n_rows != n_cols), "diagmat(): given matrix is not square" );
      
      out.zeros(n_rows, n_rows);
      
      for(uword i=0; i < n_rows; ++i) { out.at(i,i) = P.at(i,i); }
      }
    }
  else   // we have aliasing
    {
    if(P_is_vec)   // generate a diagonal matrix out of a vector
      {
      const uword N = (n_rows == 1) ? n_cols : n_rows;
      
      podarray<eT> tmp(N);
      eT* tmp_mem = tmp.memptr();
      
      if(Proxy<T1>::prefer_at_accessor == false)
        {
        typename Proxy<T1>::ea_type P_ea = P.get_ea();
        
        for(uword i=0; i < N; ++i) { tmp_mem[i] = P_ea[i]; }
        }
      else
        {
        if(n_rows == 1)
          {
          for(uword i=0; i < N; ++i) { tmp_mem[i] = P.at(0,i); }
          }
        else
          {
          for(uword i=0; i < N; ++i) { tmp_mem[i] = P.at(i,0); }
          }
        }
      
      out.zeros(N, N);
      
      for(uword i=0; i < N; ++i) { out.at(i,i) = tmp_mem[i]; }
      }
    else   // generate a diagonal matrix out of a matrix
      {
      // NOTE: we're assuming that the output matrix is the same as the matrix provided by the Proxy,
      // NOTE: and the alias is not due to a matrix using auxiliary memory;
      // NOTE: this assumption is currently valid for matrices, but not for vectors;
      // NOTE: as we've checked that at this point in code we're dealing with a matrix,
      // NOTE: the assumption is thus currently valid
      
      arma_debug_check( (n_rows != n_cols), "diagmat(): given matrix is not square" );
      
      for(uword i=0; i < n_rows; ++i)
        {
        eT* colptr = out.colptr(i);
        
        // clear above the diagonal
        arrayops::inplace_set(colptr, eT(0), i);
        
        // clear below the diagonal
        arrayops::inplace_set(colptr+(i+1), eT(0), n_rows-1-i);
        }
      }
    }
  }



//! @}
