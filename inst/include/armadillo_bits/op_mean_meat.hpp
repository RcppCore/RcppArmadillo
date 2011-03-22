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
  
  const eT result = arrayops::accumulate(X, n_elem) / T(n_elem);
  
  return arma_isfinite(result) ? result : op_mean::direct_mean_robust(X, n_elem);
  }



template<typename eT>
inline
eT
op_mean::direct_mean(const Mat<eT>& X, const u32 row)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const u32 X_n_cols = X.n_cols;
  
  eT val = eT(0);
  
  for(u32 col=0; col<X_n_cols; ++col)
    {
    val += X.at(row,col);
    }
  
  const eT result = val / T(X_n_cols);
  
  return arma_isfinite(result) ? result : direct_mean_robust(X, row);
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
  
  const eT result = val / T(X_n_elem);
  
  return arma_isfinite(result) ? result : direct_mean_robust(X);
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
  
  const eT result = val / T(X_n_elem);
  
  return arma_isfinite(result) ? result : direct_mean_robust(X);
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
      out[row] = op_mean::direct_mean( X, row );
      }
    }
  }



template<typename eT>
arma_pure
inline
eT
op_mean::direct_mean_robust(const eT* const X, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  // use an adapted form of the mean finding algorithm from the running_stat class
  
  typedef typename get_pod_type<eT>::result T;
  
  u32 i,j;
  
  eT r_mean = eT(0);
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    const eT Xi = X[i];
    const eT Xj = X[j];
    
    r_mean = r_mean + (Xi - r_mean)/T(j);    // we need i+1, and j is equivalent to i+1 here
    r_mean = r_mean + (Xj - r_mean)/T(j+1);
    }
  
  
  if(i < n_elem)
    {
    const eT Xi = X[i];
    
    r_mean = r_mean + (Xi - r_mean)/T(i+1);
    }
  
  return r_mean;
  }



template<typename eT>
inline
eT
op_mean::direct_mean_robust(const Mat<eT>& X, const u32 row)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const u32 X_n_cols = X.n_cols;
  
  eT r_mean = eT(0);
  
  for(u32 col=0; col<X_n_cols; ++col)
    {
    r_mean = r_mean + (X.at(row,col) - r_mean)/T(col+1);
    }
  
  return r_mean;
  }



template<typename eT>
inline 
eT
op_mean::direct_mean_robust(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const u32 X_n_elem = X.n_elem;
  
  eT r_mean = eT(0);
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    r_mean = r_mean + (X[i] - r_mean)/T(i+1);
    }
  
  return r_mean;
  }



template<typename eT>
inline 
eT
op_mean::direct_mean_robust(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const u32 X_n_elem = X.n_elem;
  
  eT r_mean = eT(0);
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    r_mean = r_mean + (X[i] - r_mean)/T(i+1);
    }
  
  return r_mean;
  }



//! @}

