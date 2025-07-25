// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup op_var
//! @{



template<typename T1>
inline
void
op_var::apply(Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_var>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type out_eT;
  
  const uword norm_type = in.aux_uword_a;
  const uword dim       = in.aux_uword_b;
  
  arma_conform_check( (norm_type > 1), "var(): parameter 'norm_type' must be 0 or 1" );
  arma_conform_check( (dim > 1),       "var(): parameter 'dim' must be 0 or 1"       );
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<out_eT> tmp;
    
    op_var::apply_noalias(tmp, U.M, norm_type, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_var::apply_noalias(out, U.M, norm_type, dim);
    }
  }



template<typename in_eT>
inline
void
op_var::apply_noalias(Mat<typename get_pod_type<in_eT>::result>& out, const Mat<in_eT>& X, const uword norm_type, const uword dim)
  {
  arma_debug_sigprint();
  
  typedef typename get_pod_type<in_eT>::result out_eT;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_debug_print("op_var::apply_noalias(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols);
    
    if(X_n_rows > 0)
      {
      out_eT* out_mem = out.memptr();
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        out_mem[col] = op_var::direct_var( X.colptr(col), X_n_rows, norm_type );
        }
      }
    }
  else
  if(dim == 1)
    {
    arma_debug_print("op_var::apply_noalias(): dim = 1");
    
    out.set_size(X_n_rows, (X_n_cols > 0) ? 1 : 0);
    
    if(X_n_cols > 0)
      {
      out_eT* out_mem = out.memptr();
      
      podarray<in_eT> tmp;
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        tmp.copy_row(X, row);
        
        out_mem[row] = op_var::direct_var( tmp.memptr(), tmp.n_elem, norm_type );
        }
      }
    }
  }



template<typename T1>
inline
typename T1::pod_type
op_var::var_vec(const T1& X, const uword norm_type)
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  arma_conform_check( (norm_type > 1), "var(): parameter 'norm_type' must be 0 or 1" );
  
  const quasi_unwrap<T1> U(X);
  
  if(U.M.n_elem == 0)
    {
    arma_conform_check(true, "var(): object has no elements");
    
    return Datum<T>::nan;
    }
  
  return op_var::direct_var(U.M.memptr(), U.M.n_elem, norm_type);
  }



//! find the variance of an array
template<typename eT>
inline
eT
op_var::direct_var(const eT* X, const uword n_elem, const uword norm_type)
  {
  arma_debug_sigprint();
  
  if(n_elem >= 2)
    {
    const eT acc1 = op_mean::direct_mean(X, n_elem);
    
    if(arma_isnonfinite(acc1))  { return Datum<eT>::nan; }
    
    eT acc2 = eT(0);
    eT acc3 = eT(0);
    
    uword i,j;
    
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



//! find the variance of an array (robust but slow)
template<typename eT>
inline
eT
op_var::direct_var_robust(const eT* X, const uword n_elem, const uword norm_type)
  {
  arma_debug_sigprint();
  
  if(n_elem > 1)
    {
    eT r_mean = X[0];
    eT r_var  = eT(0);
    
    for(uword i=1; i<n_elem; ++i)
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



//! find the variance of an array (version for complex numbers)
template<typename T>
inline
T
op_var::direct_var(const std::complex<T>* const X, const uword n_elem, const uword norm_type)
  {
  arma_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  if(n_elem >= 2)
    {
    const eT acc1 = op_mean::direct_mean(X, n_elem);
    
    if(arma_isnonfinite(acc1))  { return Datum<T>::nan; }
    
    T  acc2 =  T(0);
    eT acc3 = eT(0);
    
    for(uword i=0; i<n_elem; ++i)
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



//! find the variance of an array (version for complex numbers) (robust but slow)
template<typename T>
inline
T
op_var::direct_var_robust(const std::complex<T>* const X, const uword n_elem, const uword norm_type)
  {
  arma_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  if(n_elem > 1)
    {
    eT r_mean = X[0];
     T r_var  = T(0);
    
    for(uword i=1; i<n_elem; ++i)
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
