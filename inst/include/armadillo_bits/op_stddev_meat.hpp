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


//! \addtogroup op_stddev
//! @{


//! \brief
//! For each row or for each column, find the standard deviation.
//! The result is stored in a dense matrix that has either one column or one row.
//! The dimension for which the standard deviations are found is set via the stddev() function.
template<typename T1>
inline
void
op_stddev::apply(Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_stddev>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type  in_eT;
  typedef typename T1::pod_type  out_eT;
  
  const unwrap_check_mixed<T1> tmp(in.m, out);
  const Mat<in_eT>&        X = tmp.M;
  
  const u32 norm_type = in.aux_u32_a;
  const u32 dim       = in.aux_u32_b;
  
  arma_debug_check( (X.n_elem == 0), "stddev(): given matrix has no elements"             );
  arma_debug_check( (norm_type > 1), "stddev(): incorrect usage. norm_type must be 0 or 1");
  arma_debug_check( (dim > 1),       "stddev(): incorrect usage. dim must be 0 or 1"      );
  
  const u32 X_n_rows = X.n_rows;
  const u32 X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_stddev::apply(), dim = 0");
    
    out.set_size(1, X_n_cols);
    
    for(u32 col=0; col<X_n_cols; ++col)
      {
      out[col] = std::sqrt( op_var::direct_var( X.colptr(col), X_n_rows, norm_type ) );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_stddev::apply(), dim = 1");
    
    out.set_size(X_n_rows, 1);
    
    podarray<in_eT> tmp(X_n_cols);
    
    in_eT* tmp_mem = tmp.memptr();
    
    for(u32 row=0; row<X_n_rows; ++row)
      {
      for(u32 col=0; col<X_n_cols; ++col)
        {
        tmp_mem[col] = X.at(row,col);
        }
      
      out[row] = std::sqrt( op_var::direct_var(tmp_mem, X_n_cols, norm_type) );
      }
    }
  }



//! @}
