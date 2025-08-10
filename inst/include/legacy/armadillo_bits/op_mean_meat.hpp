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


//! \addtogroup op_mean
//! @{



template<typename T1>
inline
void
op_mean::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_mean>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  arma_conform_check( (dim > 1), "mean(): parameter 'dim' must be 0 or 1" );
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_mean::apply_noalias(tmp, U.M, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_mean::apply_noalias(out, U.M, dim);
    }
  }



template<typename eT>
inline
void
op_mean::apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim)
  {
  arma_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols);
    
    if(X_n_rows == 0)  { return; }
    
    eT* out_mem = out.memptr();
    
    for(uword col=0; col < X_n_cols; ++col)
      {
      out_mem[col] = op_mean::direct_mean( X.colptr(col), X_n_rows );
      }
    }
  else
  if(dim == 1)
    {
    out.zeros(X_n_rows, (X_n_cols > 0) ? 1 : 0);
    
    if(X_n_cols == 0)  { return; }
    
    eT* out_mem = out.memptr();
    
    for(uword col=0; col < X_n_cols; ++col)
      {
      arrayops::inplace_plus(out_mem, X.colptr(col), X_n_rows);
      }
    
    out /= T(X_n_cols);
    
    if(out.internal_has_nonfinite())
      {
      podarray<eT> tmp;
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        const eT old_mean = out_mem[row];
        
        if(arma_isnonfinite(old_mean))
          {
          tmp.copy_row(X, row);
          
          out_mem[row] = op_mean::direct_mean_robust(old_mean, tmp.memptr(), tmp.n_elem);
          }
        }
      }
    }
  }



//



template<typename T1>
inline
void
op_mean::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_mean>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  arma_conform_check( (dim > 2), "mean(): parameter 'dim' must be 0 or 1 or 2" );
  
  const unwrap_cube<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Cube<eT> tmp;
    
    op_mean::apply_noalias(tmp, U.M, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_mean::apply_noalias(out, U.M, dim);
    }
  }



template<typename eT>
inline
void
op_mean::apply_noalias(Cube<eT>& out, const Cube<eT>& X, const uword dim)
  {
  arma_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const uword X_n_rows   = X.n_rows;
  const uword X_n_cols   = X.n_cols;
  const uword X_n_slices = X.n_slices;
  
  if(dim == 0)
    {
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols, X_n_slices);
    
    if(X_n_rows == 0)  { return; }
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      eT* out_mem = out.slice_memptr(slice);
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        out_mem[col] = op_mean::direct_mean( X.slice_colptr(slice,col), X_n_rows );
        }
      }
    }
  else
  if(dim == 1)
    {
    out.zeros(X_n_rows, (X_n_cols > 0) ? 1 : 0, X_n_slices);
    
    if(X_n_cols == 0)  { return; }
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      eT* out_mem = out.slice_memptr(slice);
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        arrayops::inplace_plus(out_mem, X.slice_colptr(slice,col), X_n_rows);
        }
      
      for(uword row=0; row < X_n_rows; ++row)  { out_mem[row] /= T(X_n_cols); }
      
      if(arrayops::is_finite(out_mem, X_n_rows) == false)
        {
        const Mat<eT> tmp_mat('j', X.slice_memptr(slice), X_n_rows, X_n_cols);
        
        podarray<eT> tmp_vec;
        
        for(uword row=0; row < X_n_rows; ++row)
          {
          const eT old_mean = out_mem[row];
          
          if(arma_isnonfinite(old_mean))
            {
            tmp_vec.copy_row(tmp_mat, row);
            
            out_mem[row] = op_mean::direct_mean_robust(old_mean, tmp_vec.memptr(), tmp_vec.n_elem);
            }
          }
        }
      }
    }
  else
  if(dim == 2)
    {
    out.zeros(X_n_rows, X_n_cols, (X_n_slices > 0) ? 1 : 0);
    
    if(X_n_slices == 0)  { return; }
    
    eT* out_mem = out.memptr();
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      arrayops::inplace_plus(out_mem, X.slice_memptr(slice), X.n_elem_slice );
      }
    
    out /= T(X_n_slices);
    
    if(out.internal_has_nonfinite())
      {
      podarray<eT> tmp(X_n_slices);
      
      for(uword col=0; col < X_n_cols; ++col)
      for(uword row=0; row < X_n_rows; ++row)
        {
        const eT old_mean = out.at(row,col,0);
        
        if(arma_isnonfinite(old_mean))
          {
          for(uword slice=0; slice < X_n_slices; ++slice)  { tmp[slice] = X.at(row,col,slice); }
          
          out.at(row,col,0) = op_mean::direct_mean_robust(old_mean, tmp.memptr(), tmp.n_elem);
          }
        }
      }
    }
  }



//



template<typename eT>
inline
eT
op_mean::direct_mean(const eT* X_mem, const uword N)
  {
  arma_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const eT mean = arrayops::accumulate(X_mem, N) / T(N);
  
  return arma_isfinite(mean) ? mean : op_mean::direct_mean_robust(mean, X_mem, N);
  }



template<typename eT>
inline
eT
op_mean::direct_mean_robust(const eT old_mean, const eT* X_mem, const uword N)
  {
  arma_debug_sigprint();
  
  // use an adapted form of the mean finding algorithm from the running_stat class
  
  typedef typename get_pod_type<eT>::result T;
  
  if(arrayops::is_finite(X_mem, N) == false)  { return old_mean; }
  
  eT r_mean = eT(0);
  
  for(uword i=0; i < N; ++i)
    {
    r_mean = r_mean + (X_mem[i] - r_mean) / T(i+1);
    }
  
  return r_mean;
  }



//



template<typename T1>
inline
typename T1::elem_type
op_mean::mean_all(const T1& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(X);
  
  if(U.M.n_elem == 0)
    {
    arma_conform_check(true, "mean(): object has no elements");
    
    return Datum<eT>::nan;
    }
  
  return op_mean::direct_mean(U.M.memptr(), U.M.n_elem);
  }



template<typename T1>
inline
typename T1::elem_type
op_mean::mean_all(const Op<T1, op_omit>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword omit_mode = in.aux_uword_a;
  
  if(arma_config::fast_math_warn)
    {
    if(omit_mode == 1)  { arma_warn(1, "omit_nan(): detection of NaN is not reliable in fast math mode"); }
    if(omit_mode == 2)  { arma_warn(1, "omit_nonfinite(): detection of non-finite values is not reliable in fast math mode"); }
    }
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.M.n_elem == 0)
    {
    arma_conform_check(true, "mean(): object has no elements");
    
    return Datum<eT>::nan;
    }
  
  auto is_omitted_1 = [](const eT& x) -> bool { return arma_isnan(x);       };
  auto is_omitted_2 = [](const eT& x) -> bool { return arma_isnonfinite(x); };
  
  eT result = eT(0);
  
  if(omit_mode == 1)  { result = op_mean::mean_all_omit(U.M.memptr(), U.M.n_elem, is_omitted_1); }
  if(omit_mode == 2)  { result = op_mean::mean_all_omit(U.M.memptr(), U.M.n_elem, is_omitted_2); }
  
  return result;
  }



template<typename eT, typename functor>
inline
eT
op_mean::mean_all_omit(const eT* X_mem, const uword N, functor is_omitted)
  {
  arma_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  uword count = 0;
  eT    acc   = eT(0);
  
  for(uword i=0; i < N; ++i)
    {
    const eT val = X_mem[i];
    
    if(is_omitted(val) == false)  { acc += val; ++count; }
    }
  
  acc /= T(count);
  
  if(arma_isfinite(acc))  { return acc; }
  
  // handle possible overflow
  
  eT r_mean = eT(0);
  
  count = 0;
  
  for(uword i=0; i < N; ++i)
    {
    const eT val = X_mem[i];
    
    if(is_omitted(val) == false)
      {
      r_mean = r_mean + (val - r_mean) / T(count+1);  // kept as count+1 to use same algorithm as op_mean::direct_mean_robust()
      
      ++count;
      }
    }
  
  return r_mean;
  }



//



template<typename eT>
arma_inline
eT
op_mean::robust_mean(const eT A, const eT B)
  {
  return (arma_isfinite(A) && arma_isfinite(B)) ? eT( A + (B - A)/eT(2) ) : eT( (A+B)/eT(2) );
  }



template<typename T>
arma_inline
std::complex<T>
op_mean::robust_mean(const std::complex<T>& A, const std::complex<T>& B)
  {
  typedef typename std::complex<T> eT;
  
  return (arma_isfinite(A) && arma_isfinite(B)) ? eT( A + (B - A)/T(2) ) : eT( (A+B)/T(2) );
  }



//! @}
