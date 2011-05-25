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


//! \addtogroup op_min
//! @{



template<typename eT>
arma_pure
inline 
eT
op_min::direct_min(const eT* const X, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  eT min_val = (n_elem != 1) ? priv::most_pos<eT>() : X[0];
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    const eT X_i = X[i];
    const eT X_j = X[j];
    
    if(X_i < min_val)
      {
      min_val = X_i;
      }
    
    if(X_j < min_val)
      {
      min_val = X_j;
      }
    }
  
  
  if(i < n_elem)
    {
    const eT X_i = X[i];
    
    if(X_i < min_val)
      {
      min_val = X_i;
      }
    }
  
  return min_val;
  }



template<typename eT>
inline 
eT
op_min::direct_min(const eT* const X, const u32 n_elem, u32& index_of_min_val)
  {
  arma_extra_debug_sigprint();
  
  eT min_val = (n_elem != 1) ? priv::most_pos<eT>() : X[0];
  
  u32 best_index = 0;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    const eT X_i = X[i];
    const eT X_j = X[j];
    
    if(X_i < min_val)
      {
      min_val    = X_i;
      best_index = i;
      }
    
    if(X_j < min_val)
      {
      min_val    = X_j;
      best_index = j;
      }
    }
  
  
  if(i < n_elem)
    {
    const eT X_i = X[i];
    
    if(X_i < min_val)
      {
      min_val    = X_i;
      best_index = i;
      }
    }
  
  index_of_min_val = best_index;
  
  return min_val;
  }  



template<typename eT>
inline 
eT
op_min::direct_min(const Mat<eT>& X, const u32 row)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_cols = X.n_cols;
  
  eT min_val = (X_n_cols != 1) ? priv::most_pos<eT>() : X.at(row,0);
  
  for(u32 col=0; col<X_n_cols; ++col)
    {
    const eT tmp_val = X.at(row,col);
    
    if(tmp_val < min_val)
      {
      min_val = tmp_val;
      }
    }
  
  return min_val;
  }



template<typename eT>
inline 
eT
op_min::direct_min(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_elem = X.n_elem;
  
  eT min_val = (X_n_elem != 1) ? priv::most_pos<eT>() : X[0];
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    eT tmp_val = X[i];
    
    if(tmp_val < min_val)
      {
      min_val = tmp_val;
      }
    }
  
  return min_val;
  }



template<typename eT>
inline 
eT
op_min::direct_min(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_elem = X.n_elem;
  
  eT min_val = (X_n_elem != 1) ? priv::most_pos<eT>() : X[0];
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    eT tmp_val = X[i];
    
    if(tmp_val < min_val)
      {
      min_val = tmp_val;
      }
    }
  
  return min_val;
  }



//! \brief
//! For each row or for each column, find the minimum value.
//! The result is stored in a dense matrix that has either one column or one row.
//! The dimension, for which the minima are found, is set via the min() function.
template<typename T1>
inline void op_min::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_min>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>&     X = tmp.M;
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "min(): incorrect usage. dim must be 0 or 1");
  
  const u32 X_n_rows = X.n_rows;
  const u32 X_n_cols = X.n_cols;
  
  if(dim == 0)  // min in each column
    {
    arma_extra_debug_print("op_min::apply(), dim = 0");
    
    out.set_size( (X_n_rows > 0) ? 1 : 0, X_n_cols );
    
    if(X_n_rows > 0)
      {
      eT* out_mem = out.memptr();
      
      for(u32 col=0; col<X_n_cols; ++col)
        {
        out_mem[col] = op_min::direct_min( X.colptr(col), X_n_rows );
        }
      }
    }
  else
  if(dim == 1)  // min in each row
    {
    arma_extra_debug_print("op_min::apply(), dim = 1");
    
    out.set_size( X_n_rows, (X_n_cols > 0) ? 1 : 0 );
    
    if(X_n_cols > 0)
      {
      eT* out_mem = out.memptr();
      
      for(u32 row=0; row<X_n_rows; ++row)
        {
        out_mem[row] = op_min::direct_min( X, row );
        }
      }
    }
  }



template<typename T>
inline 
std::complex<T>
op_min::direct_min(const std::complex<T>* const X, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  u32 index   = 0;
  T   min_val = (n_elem != 1) ? priv::most_pos<T>() : std::abs(X[0]);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const T tmp_val = std::abs(X[i]);
    
    if(tmp_val < min_val)
      {
      min_val = tmp_val;
      index   = i;
      }
    }
  
  return X[index];
  }



template<typename T>
inline 
std::complex<T>
op_min::direct_min(const std::complex<T>* const X, const u32 n_elem, u32& index_of_min_val)
  {
  arma_extra_debug_sigprint();
  
  u32 index   = 0;
  T   min_val = (n_elem != 1) ? priv::most_pos<T>() : std::abs(X[0]);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const T tmp_val = std::abs(X[i]);
    
    if(tmp_val < min_val)
      {
      min_val = tmp_val;
      index   = i;
      }
    }
  
  index_of_min_val = index;
  
  return X[index];
  }



template<typename T>
inline 
std::complex<T>
op_min::direct_min(const Mat< std::complex<T> >& X, const u32 row)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_cols = X.n_cols;
  
  u32 index   = 0;
  T   min_val = (X_n_cols != 1) ? priv::most_pos<T>() : std::abs(X.at(row,0));
  
  for(u32 col=0; col<X_n_cols; ++col)
    {
    const T tmp_val = std::abs(X.at(row,col));
    
    if(tmp_val < min_val)
      {
      min_val = tmp_val;
      index   = col;
      }
    }
  
  return X.at(row,index);
  }



template<typename T>
inline 
std::complex<T>
op_min::direct_min(const subview< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_elem = X.n_elem;
        u32 index    = 0;
        T   min_val  = (X_n_elem != 1) ? priv::most_pos<T>() : std::abs(X[0]);
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    const T tmp_val = std::abs(X[i]);
    
    if(tmp_val < min_val)
      {
      min_val = tmp_val;
      index   = i;
      }
    }
  
  return X[index];
  }



template<typename T>
inline 
std::complex<T>
op_min::direct_min(const diagview< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  const u32 X_n_elem = X.n_elem;
        u32 index    = 0;
        T   min_val  = (X_n_elem != 1) ? priv::most_pos<T>() : std::abs(X[0]);
  
  for(u32 i=0; i<X_n_elem; ++i)
    {
    const T tmp_val = std::abs(X[i]);
    
    if(tmp_val < min_val)
      {
      min_val = tmp_val;
      index   = i;
      }
    }
  
  return X[index];
  }



//! @}
