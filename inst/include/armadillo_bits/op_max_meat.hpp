// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_max
//! @{


//! find the maximum value in an array
template<typename eT>
arma_pure
inline 
eT
op_max::direct_max(const eT* const X, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  eT max_val = (n_elem != 1) ? priv::most_neg<eT>() : X[0];
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    const eT X_i = X[i];
    const eT X_j = X[j];
    
    if(X_i > max_val)
      {
      max_val = X_i;
      }
    
    if(X_j > max_val)
      {
      max_val = X_j;
      }
    }
  
  
  if(i < n_elem)
    {
    const eT X_i = X[i];
    
    if(X_i > max_val)
      {
      max_val = X_i;
      }
    }
  
  return max_val;
  }



//! find the maximum value in a subview
template<typename eT>
inline 
eT
op_max::direct_max(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_elem = X.n_elem;
        eT  max_val  = (X_n_elem != 1) ? priv::most_neg<eT>() : X[0];
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    eT tmp_val = X[i];
    
    if(tmp_val > max_val)
      {
      max_val = tmp_val;
      }
    }
  
  return max_val;
  }



//! find the maximum value in a diagview
template<typename eT>
inline 
eT
op_max::direct_max(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_elem = X.n_elem;
        eT  max_val  = (X_n_elem != 1) ? priv::most_neg<eT>() : X[0];
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    eT tmp_val = X[i];
    
    if(tmp_val > max_val)
      {
      max_val = tmp_val;
      }
    }
  
  return max_val;
  }



//! \brief
//! For each row or for each column, find the maximum value.
//! The result is stored in a dense matrix that has either one column or one row.
//! The dimension, for which the maxima are found, is set via the max() function.
template<typename T1>
inline
void
op_max::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_max>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& X     = tmp.M;
  
  arma_debug_check( (X.n_elem == 0), "max(): given matrix has no elements" );
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "max(): incorrect usage. dim must be 0 or 1");
  
  const u32 X_n_rows = X.n_rows;
  const u32 X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_max::apply(), dim = 0");
    
    out.set_size(1, X_n_cols);
    
    for(u32 col=0; col<X_n_cols; ++col)
      {
      out[col] = op_max::direct_max( X.colptr(col), X_n_rows );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_max::apply(), dim = 1");
    
    out.set_size(X_n_rows, 1);
    
    for(u32 row=0; row<X_n_rows; ++row)
      {
      eT max_val = (X_n_cols != 1) ? priv::most_neg<eT>() : X.at(row,0);
      
      for(u32 col=0; col<X_n_cols; ++col)
        {
        const eT tmp_val = X.at(row,col);
        
        if(tmp_val > max_val)
          {
          max_val = tmp_val;
          }
        }
      
      out[row] = max_val;
      }
    }
  }



//! Find the maximum value in an array (version for complex numbers)
template<typename T>
inline 
std::complex<T>
op_max::direct_max(const std::complex<T>* const X, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  u32 index   = 0;
  T   max_val = (n_elem != 1) ? priv::most_neg<T>() : std::abs(X[0]);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const T tmp_val = std::abs(X[i]);
    
    if(tmp_val > max_val)
      {
      max_val = tmp_val;
      index   = i;
      }
    }
  
  return X[index];
  }



//! Find the maximum value in a subview (version for complex numbers)
template<typename T>
inline 
std::complex<T>
op_max::direct_max(const subview< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_elem = X.n_elem;
        u32 index    = 0;
        T   max_val  = (X_n_elem != 1) ? priv::most_neg<T>() : std::abs(X[0]);
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    const T tmp_val = std::abs(X[i]);
    
    if(tmp_val > max_val)
      {
      max_val = tmp_val;
      index   = i;
      }
    }
  
  return X[index];
  }



//! Find the maximum value in a diagview (version for complex numbers)
template<typename T>
inline 
std::complex<T>
op_max::direct_max(const diagview< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_elem = X.n_elem;
        u32 index    = 0;
        T   max_val  = (X_n_elem != 1) ? priv::most_neg<T>() : std::abs(X[0]);
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    const T tmp_val = std::abs(X[i]);
    
    if(tmp_val > max_val)
      {
      max_val = tmp_val;
      index   = i;
      }
    }
  
  return X[index];
  }



//! Implementation for complex numbers
template<typename T, typename T1>
inline void op_max::apply(Mat< std::complex<T> >& out, const Op<T1,op_max>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem == 0), "max(): given matrix has no elements" );
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "max(): incorrect usage. dim must be 0 or 1");
  
  const u32 X_n_rows = X.n_rows;
  const u32 X_n_cols = X.n_cols;
  
  if(dim == 0)  // column-wise max
    {
    arma_extra_debug_print("op_max::apply(), dim = 0");
    
    out.set_size(1, X_n_cols);
    
    for(u32 col=0; col<X_n_cols; ++col)
      {
      out[col] = op_max::direct_max( X.colptr(col), X_n_rows );
      }
    }
  else
  if(dim == 1)  // row-wise max
    {
    arma_extra_debug_print("op_max::apply(), dim = 1");
    
    out.set_size(X_n_rows, 1);
    
    for(u32 row=0; row<X_n_rows; ++row)
      {
      u32 index   = 0;
      T   max_val = (X_n_cols != 1) ? priv::most_neg<T>() : std::abs(X.at(row,0));
      
      for(u32 col=0; col<X_n_cols; ++col)
        {
        const T tmp_val = std::abs(X.at(row,col));
        
        if(tmp_val > max_val)
          {
          max_val = tmp_val;
          index   = col;
          }
        }
      
      out[row] = X.at(row,index);
      }
    
    }
  
  }



//! @}
