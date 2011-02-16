// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_prod
//! @{

//! \brief
//! Immediate product of elements of a matrix along a specified dimension (either rows or columns).
//! The result is stored in a dense matrix that has either one column or one row.
//! See the prod() function for more details.
template<typename T1>
inline
void
op_prod::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_prod>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "prod(): incorrect usage. dim must be 0 or 1");
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& X     = tmp.M;
  
  arma_debug_check( (X.n_elem < 1), "prod(): give object has no elements");
  
  const u32 X_n_rows = X.n_rows;
  const u32 X_n_cols = X.n_cols;
    
  if(dim == 0)  // traverse across rows (i.e. find the product in each column)
    {
    out.set_size(1, X_n_cols);
    
    for(u32 col=0; col<X_n_cols; ++col)
      {
      out.at(0,col) = arrayops::product(X.colptr(col), X_n_rows);
      }
    }
  else  // traverse across columns (i.e. find the product in each row)
    {
    out.set_size(X_n_rows, 1);
    
    for(u32 row=0; row < X_n_rows; ++row)
      {
      eT val = X.at(row,0);
      
      for(u32 col=1; col < X_n_cols; ++col)
        {
        val *= X.at(row,col);
        }
    
      out.at(row,0) = val;
      }
    
    }
  
  }



//! @}
