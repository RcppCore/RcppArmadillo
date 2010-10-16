// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_stddev
//! @{


//! \brief
//! For each row or for each column, find the standard deviation.
//! The result is stored in a dense matrix that has either one column or one row.
//! The dimension for which the standard deviations are found is set via the stddev() function.
template<typename eT>
inline
void
op_stddev::apply(Mat< typename get_pod_type<eT>::result >& out, const Mat<eT>& X, const u32 norm_type, const u32 dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (X.n_elem == 0), "stddev(): given matrix has no elements" );
  
  arma_debug_check( (norm_type > 1), "stddev(): incorrect usage. norm_type must be 0 or 1");
  arma_debug_check( (dim > 1),       "stddev(): incorrect usage. dim must be 0 or 1"      );
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_stddev::apply(), dim = 0");
    
    out.set_size(1, X.n_cols);
    
    for(u32 col=0; col<X.n_cols; ++col)
      {
      out[col] = std::sqrt( op_var::direct_var( X.colptr(col), X.n_rows, norm_type ) );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_stddev::apply(), dim = 1");
    
    const u32 n_rows = X.n_rows;
    const u32 n_cols = X.n_cols;
    
    out.set_size(n_rows, 1);
    
    podarray<eT> tmp(n_cols);
    
    eT* tmp_mem = tmp.memptr();
    
    for(u32 row=0; row<n_rows; ++row)
      {
      for(u32 col=0; col<n_cols; ++col)
        {
        tmp_mem[col] = X.at(row,col);
        }
      
      out[row] = std::sqrt( op_var::direct_var(tmp_mem, n_cols, norm_type) );
      }
    
    }
  
  }



//! @}
