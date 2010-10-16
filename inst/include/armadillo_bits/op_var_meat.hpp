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


//! \addtogroup op_var
//! @{


//! find the variance of an array
template<typename eT>
inline
eT
op_var::direct_var(const eT* const X, const u32 n_elem, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  eT acc1 = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    acc1 += X[i];
    }
  
  const eT div_val = (n_elem > 0) ? eT(n_elem) : eT(1);
  acc1 /= div_val;
  
  eT acc2 = eT(0);
  eT acc3 = eT(0);

  for(u32 i=0; i<n_elem; ++i)
    {
    const eT tmp = acc1 - X[i];
  
    acc2 += tmp*tmp;
    acc3 += tmp;
    }
  
  
  const eT norm_val = (norm_type == 0) ? ( (n_elem > 1) ? eT(n_elem-1) : eT(1) ) : eT(n_elem);
  const eT var_val  = (acc2 - acc3*acc3/div_val) / norm_val;
  
  return var_val;
  }



//! find the variance of an array (version for complex numbers)
template<typename T>
inline
T
op_var::direct_var(const std::complex<T>* const X, const u32 n_elem, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  eT acc1 = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    acc1 += X[i];
    }
  
  const T div_val = (n_elem > 0) ? T(n_elem) : T(1);
  acc1 /= div_val;
  
  T  acc2 =  T(0);
  eT acc3 = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const eT tmp = acc1 - X[i];
    
    acc2 += std::norm(tmp);
    acc3 += tmp;
    }
  
  const T norm_val = (norm_type == 0) ? ( (n_elem > 1) ? T(n_elem-1) : T(1) ) : T(n_elem);
  const T var_val  = (acc2 - std::norm(acc3)/div_val) / norm_val;
  
  return var_val;
  }



//! find the variance of a subview_row
template<typename eT>
inline 
typename get_pod_type<eT>::result
op_var::direct_var(const subview_row<eT>& X, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  const u32 n_elem = X.n_elem;
  
  podarray<eT> tmp(n_elem);
  
  eT* tmp_mem = tmp.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    tmp_mem[i] = X[i];
    }
  
  return op_var::direct_var(tmp_mem, n_elem, norm_type);
  }



//! find the variance of a subview_col
template<typename eT>
inline 
typename get_pod_type<eT>::result
op_var::direct_var(const subview_col<eT>& X, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  return op_var::direct_var(X.colptr(0), X.n_elem, norm_type);
  }



//! find the variance of a diagview
template<typename eT>
inline 
typename get_pod_type<eT>::result
op_var::direct_var(const diagview<eT>& X, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  const u32 n_elem = X.n_elem;
  
  podarray<eT> tmp(n_elem);
  
  eT* tmp_mem = tmp.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    tmp_mem[i] = X[i];
    }
  
  return op_var::direct_var(tmp_mem, n_elem, norm_type);
  }



//! \brief
//! For each row or for each column, find the variance.
//! The result is stored in a dense matrix that has either one column or one row.
//! The dimension, for which the variances are found, is set via the var() function.
template<typename eT>
inline
void
op_var::apply(Mat< typename get_pod_type<eT>::result >& out, const Mat<eT>& X, const u32 norm_type, const u32 dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (X.n_elem == 0), "var(): given matrix has no elements" );
  
  arma_debug_check( (norm_type > 1), "var(): incorrect usage. norm_type must be 0 or 1");
  arma_debug_check( (dim > 1),       "var(): incorrect usage. dim must be 0 or 1"      );
  
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_var::apply(), dim = 0");
    
    out.set_size(1, X.n_cols);
    
    for(u32 col=0; col<X.n_cols; ++col)
      {
      out[col] = op_var::direct_var( X.colptr(col), X.n_rows, norm_type );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_var::apply(), dim = 1");
    
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
      
      out[row] = op_var::direct_var(tmp_mem, n_cols, norm_type);
      }
    
    }
  
  }



//! @}
