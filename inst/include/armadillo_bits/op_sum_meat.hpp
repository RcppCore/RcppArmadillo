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


//! \addtogroup op_sum
//! @{

//! \brief
//! Immediate sum of elements of a matrix along a specified dimension (either rows or columns).
//! The result is stored in a dense matrix that has either one column or one row.
//! See the sum() function for more details.
template<typename T1>
inline
void
op_sum::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sum>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  arma_debug_check( (dim > 1), "sum(): incorrect usage. dim must be 0 or 1");
  
  const Proxy<T1> P(in.m);
  
  const bool is_alias = P.is_alias(out);
  
  if( (is_Mat< typename Proxy<T1>::stored_type>::value == true) || is_alias )
    {
    const unwrap_check< typename Proxy<T1>::stored_type > tmp(P.Q, is_alias);
    
    const Mat<eT>& X = tmp.M;
    
    const uword X_n_rows = X.n_rows;
    const uword X_n_cols = X.n_cols;
    
    if(dim == 0)  // traverse across rows (i.e. find the sum in each column)
      {
      out.set_size(1, X_n_cols);
      
      eT* out_mem = out.memptr();
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        out_mem[col] = arrayops::accumulate( X.colptr(col), X_n_rows );
        }
      }
    else  // traverse across columns (i.e. find the sum in each row)
      {
      out.set_size(X_n_rows, 1);
      
      eT* out_mem = out.memptr();
        
      for(uword row=0; row < X_n_rows; ++row)
        {
        eT val = eT(0);
        
        uword i,j;
        for(i=0, j=1; j < X_n_cols; i+=2, j+=2)
          {
          val += X.at(row,i);
          val += X.at(row,j);
          }
        
        if(i < X_n_cols)
          {
          val += X.at(row,i);
          }
        
        out_mem[row] = val;
        }
      }
    }
  else
    {
    const uword P_n_rows = P.get_n_rows();
    const uword P_n_cols = P.get_n_cols();
    
    if(dim == 0)  // traverse across rows (i.e. find the sum in each column)
      {
      out.set_size(1, P_n_cols);
      
      eT* out_mem = out.memptr();
      
      for(uword col=0; col < P_n_cols; ++col)
        {
        eT val = eT(0);
        
        uword i,j;
        for(i=0, j=1; j < P_n_rows; i+=2, j+=2)
          {
          val += P.at(i,col);
          val += P.at(j,col);
          }
        
        if(i < P_n_rows)
          {
          val += P.at(i,col);
          }
        
        out_mem[col] = val;
        }
      }
    else  // traverse across columns (i.e. find the sum in each row)
      {
      out.set_size(P_n_rows, 1);
      
      eT* out_mem = out.memptr();
      
      for(uword row=0; row < P_n_rows; ++row)
        {
        eT val = eT(0);
        
        uword i,j;
        for(i=0, j=1; j < P_n_cols; i+=2, j+=2)
          {
          val += P.at(row,i);
          val += P.at(row,j);
          }
        
        if(i < P_n_cols)
          {
          val += P.at(row,i);
          }
        
        out_mem[row] = val;
        }
      }
    }
  }



//! @}
