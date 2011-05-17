// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_cumsum
//! @{


template<typename T1>
inline
void
op_cumsum_mat::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cumsum_mat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "cumsum(): incorrect usage. dim must be 0 or 1");
  
  out.copy_size(X);
  
  const u32 X_n_rows = X.n_rows;
  const u32 X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_cumsum::apply(), dim = 0");
    
    for(u32 col=0; col<X_n_cols; ++col)
      {
            eT* out_colmem = out.colptr(col);
      const eT* X_colmem   = X.colptr(col);
      
      eT acc = eT(0);
      
      for(u32 row=0; row<X_n_rows; ++row)
        {
        acc += X_colmem[row];
        
        out_colmem[row] = acc;
        }
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_cumsum::apply(), dim = 1");
    
    for(u32 row=0; row<X_n_rows; ++row)
      {
      eT acc = eT(0);
      
      for(u32 col=0; col<X_n_cols; ++col)
        {
        acc += X.at(row,col);
        
        out.at(row,col) = acc;
        }
      }
    }
  }



template<typename T1>
inline
void
op_cumsum_vec::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cumsum_vec>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  const u32 n_elem = X.n_elem;
  
  out.copy_size(X);
  
        eT* out_mem = out.memptr();
  const eT* X_mem   = X.memptr();
  
  eT acc = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    acc += X_mem[i];
    
    out_mem[i] = acc;
    }
  }



//! @}

