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


//! \addtogroup op_var
//! @{


//! find the variance of an array
template<typename eT>
inline
eT
op_var::direct_var(const eT* const X, const u32 n_elem, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  if(n_elem > 1)
    {
    const eT acc1 = op_mean::direct_mean(X, n_elem);
    
    eT acc2 = eT(0);
    eT acc3 = eT(0);
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT Xi = X[i];
      const eT Xj = X[j];
      
      const eT tmpi = acc1 - Xi;
      const eT tmpj = acc1 - Xj;
      
      acc2 += tmpi*tmpi + tmpj*tmpj;
      acc3 += tmpi      + tmpj;
      }
    
    if(i < n_elem)
      {
      const eT Xi = X[i];
      
      const eT tmpi = acc1 - Xi;
      
      acc2 += tmpi*tmpi;
      acc3 += tmpi;
      }
    
    const eT norm_val = (norm_type == 0) ? eT(n_elem-1) : eT(n_elem);
    const eT var_val  = (acc2 - acc3*acc3/eT(n_elem)) / norm_val;
    
    return arma_isfinite(var_val) ? var_val : op_var::direct_var_robust(X, n_elem, norm_type);
    }
  else
    {
    return eT(0);
    }
  }



//! find the variance of an array (version for complex numbers)
template<typename T>
inline
T
op_var::direct_var(const std::complex<T>* const X, const u32 n_elem, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  if(n_elem > 1)
    {
    const eT acc1 = op_mean::direct_mean(X, n_elem);
    
    T  acc2 =  T(0);
    eT acc3 = eT(0);
    
    for(u32 i=0; i<n_elem; ++i)
      {
      const eT tmp = acc1 - X[i];
      
      acc2 += std::norm(tmp);
      acc3 += tmp;
      }
    
    const T norm_val = (norm_type == 0) ? T(n_elem-1) : T(n_elem);
    const T var_val  = (acc2 - std::norm(acc3)/T(n_elem)) / norm_val;
    
    return arma_isfinite(var_val) ? var_val : op_var::direct_var_robust(X, n_elem, norm_type);
    }
  else
    {
    return T(0);
    }
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
template<typename T1>
inline
void
op_var::apply(Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_var>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type  in_eT;
  typedef typename T1::pod_type  out_eT;
  
  const unwrap_check_mixed<T1> tmp(in.m, out);
  const Mat<in_eT>&        X = tmp.M;
  
  const u32 norm_type = in.aux_u32_a;
  const u32 dim       = in.aux_u32_b;
  
  arma_debug_check( (X.n_elem == 0), "var(): given matrix has no elements"             );
  arma_debug_check( (norm_type > 1), "var(): incorrect usage. norm_type must be 0 or 1");
  arma_debug_check( (dim > 1),       "var(): incorrect usage. dim must be 0 or 1"      );
  
  const u32 X_n_rows = X.n_rows;
  const u32 X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_var::apply(), dim = 0");
    
    out.set_size(1, X_n_cols);
    
    for(u32 col=0; col<X_n_cols; ++col)
      {
      out[col] = op_var::direct_var( X.colptr(col), X_n_rows, norm_type );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_var::apply(), dim = 1");
    
    out.set_size(X_n_rows, 1);
    
    podarray<in_eT> tmp(X_n_cols);
    
    in_eT* tmp_mem = tmp.memptr();
    
    for(u32 row=0; row<X_n_rows; ++row)
      {
      for(u32 col=0; col<X_n_cols; ++col)
        {
        tmp_mem[col] = X.at(row,col);
        }
      
      out[row] = op_var::direct_var(tmp_mem, X_n_cols, norm_type);
      }
    }
  }



//! find the variance of an array (robust but slow)
template<typename eT>
inline
eT
op_var::direct_var_robust(const eT* const X, const u32 n_elem, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  if(n_elem > 1)
    {
    eT r_mean = X[0];
    eT r_var  = eT(0);
    
    for(u32 i=1; i<n_elem; ++i)
      {
      const eT tmp      = X[i] - r_mean;
      const eT i_plus_1 = eT(i+1);
      
      r_var  = eT(i-1)/eT(i) * r_var + (tmp*tmp)/i_plus_1;
      
      r_mean = r_mean + tmp/i_plus_1;
      }
    
    return (norm_type == 0) ? r_var : (eT(n_elem-1)/eT(n_elem)) * r_var;
    }
  else
    {
    return eT(0);
    }
  }



//! find the variance of an array (version for complex numbers) (robust but slow)
template<typename T>
inline
T
op_var::direct_var_robust(const std::complex<T>* const X, const u32 n_elem, const u32 norm_type)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  if(n_elem > 1)
    {
    eT r_mean = X[0];
     T r_var  = T(0);
    
    for(u32 i=1; i<n_elem; ++i)
      {
      const eT tmp      = X[i] - r_mean;
      const  T i_plus_1 = T(i+1);
      
      r_var  = T(i-1)/T(i) * r_var + std::norm(tmp)/i_plus_1;
      
      r_mean = r_mean + tmp/i_plus_1;
      }
    
    return (norm_type == 0) ? r_var : (T(n_elem-1)/T(n_elem)) * r_var;
    }
  else
    {
    return T(0);
    }
  }



//! @}

