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


//! \addtogroup op_mean
//! @{


//! find the mean value of an array
template<typename eT>
arma_pure
inline
eT
op_mean::direct_mean(const eT* const X, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  return arrayops::accumulate(X, n_elem) / T(n_elem);
  }



//! find the mean value of a subview
template<typename eT>
inline 
eT
op_mean::direct_mean(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const u32 X_n_elem = X.n_elem;
        eT  val      = eT(0);
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    val += X[i];
    }
  
  return val / T(X_n_elem);
  }



//! find the mean value of a diagview
template<typename eT>
inline 
eT
op_mean::direct_mean(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const u32 X_n_elem = X.n_elem;
        eT  val      = eT(0);
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    val += X[i];
    }
  
  return val / T(X_n_elem);
  }



//! \brief
//! For each row or for each column, find the mean value.
//! The result is stored in a dense matrix that has either one column or one row.
//! The dimension, for which the means are found, is set via the mean() function.
template<typename T1>
inline
void
op_mean::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_mean>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type            eT;
  typedef typename get_pod_type<eT>::result  T;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem == 0), "mean(): given matrix has no elements" );
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "mean(): incorrect usage. dim must be 0 or 1");
  
  const u32 X_n_rows = X.n_rows;
  const u32 X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_mean::apply(), dim = 0");
    
    out.set_size(1, X_n_cols);
    
    for(u32 col=0; col<X_n_cols; ++col)
      {
      out[col] = op_mean::direct_mean( X.colptr(col), X_n_rows );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_mean::apply(), dim = 1");
    
    out.set_size(X_n_rows, 1);
    
    for(u32 row=0; row<X_n_rows; ++row)
      {
      eT val = eT(0);
      
      for(u32 col=0; col<X_n_cols; ++col)
        {
        val += X.at(row,col);
        }
      
      out[row] = val / T(X_n_cols);
      }
    }
  }


//! @}

