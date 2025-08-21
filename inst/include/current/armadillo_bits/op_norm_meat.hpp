// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup op_norm
//! @{



template<typename T1>
inline
typename T1::elem_type
op_norm::vec_norm_1(const Proxy<T1>& P, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  constexpr bool use_direct_mem = (is_Mat<typename Proxy<T1>::stored_type>::value) || (is_subview_col<typename Proxy<T1>::stored_type>::value) || (arma_config::openmp && Proxy<T1>::use_mp);
  
  if(use_direct_mem)
    {
    const quasi_unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    
    return op_norm::vec_norm_1_direct_std(tmp.M);
    }
  
  typedef typename T1::elem_type eT;
  
  eT acc = eT(0);
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    
    uword i,j;
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      acc1 += std::abs(A[i]);
      acc2 += std::abs(A[j]);
      }
    
    if(i < N)
      {
      acc1 += std::abs(A[i]);
      }
    
    acc = (acc1 + acc2);
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows == 1)
      {
      for(uword col=0; col<n_cols; ++col)
        {
        acc += std::abs(P.at(0,col));
        }
      }
    else
      {
      eT acc1 = eT(0);
      eT acc2 = eT(0);
      
      for(uword col=0; col<n_cols; ++col)
        {
        uword i,j;
        
        for(i=0, j=1; j<n_rows; i+=2, j+=2)
          {
          acc1 += std::abs(P.at(i,col));
          acc2 += std::abs(P.at(j,col));
          }
        
        if(i < n_rows)
          {
          acc1 += std::abs(P.at(i,col));
          }
        }
      
      acc = (acc1 + acc2);
      }
    }
  
  return acc;
  }



template<typename T1>
inline
typename T1::pod_type
op_norm::vec_norm_1(const Proxy<T1>& P, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  T acc = T(0);
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    for(uword i=0; i<N; ++i)
      {
      const std::complex<T>& X = A[i];
      
      const T a = X.real();
      const T b = X.imag();
      
      acc += std::sqrt( (a*a) + (b*b) );
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows == 1)
      {
      for(uword col=0; col<n_cols; ++col)
        {
        const std::complex<T>& X = P.at(0,col);
        
        const T a = X.real();
        const T b = X.imag();
        
        acc += std::sqrt( (a*a) + (b*b) );
        }
      }
    else
      {
      for(uword col=0; col<n_cols; ++col)
      for(uword row=0; row<n_rows; ++row)
        {
        const std::complex<T>& X = P.at(row,col);
        
        const T a = X.real();
        const T b = X.imag();
        
        acc += std::sqrt( (a*a) + (b*b) );
        }
      }
    }
  
  if( (acc != T(0)) && arma_isfinite(acc) )
    {
    return acc;
    }
  else
    {
    arma_debug_print("detected possible underflow or overflow");
    
    typedef typename promote_type<eT, float>::result acc_eT;
    
    typedef typename get_pod_type<acc_eT>::result acc_T;
    
    arma_debug_type_print<eT>("eT");
    arma_debug_type_print<acc_eT>("acc_eT");
    
    const quasi_unwrap<typename Proxy<T1>::stored_type> R(P.Q);
    
    const uword N     = R.M.n_elem;
    const eT*   R_mem = R.M.memptr();
    
    acc_T max_val = priv::most_neg<acc_T>();
    
    for(uword i=0; i<N; ++i)
      {
      const std::complex<T>& X = R_mem[i];
      
      const acc_T a = std::abs(acc_T(X.real()));
      const acc_T b = std::abs(acc_T(X.imag()));
      
      if(a > max_val)  { max_val = a; }
      if(b > max_val)  { max_val = b; }
      }
    
    if(max_val == acc_T(0))  { return T(0); }
    
    acc_T alt_acc = acc_T(0);
    
    for(uword i=0; i<N; ++i)
      {
      const std::complex<T>& X = R_mem[i];
      
      const acc_T a = acc_T(X.real()) / max_val;
      const acc_T b = acc_T(X.imag()) / max_val;
      
      alt_acc += std::sqrt( (a*a) + (b*b) );
      }
    
    return T( alt_acc * max_val );
    }
  }



template<typename eT>
inline
eT
op_norm::vec_norm_1_direct_std(const Mat<eT>& X, const typename arma_blas_real_only<eT>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = X.n_elem;
  const eT*   A = X.memptr();
  
  eT out_val = eT(0);
  
  #if defined(ARMA_USE_ATLAS)
    {
    arma_debug_print("atlas::cblas_asum()");
    out_val = atlas::cblas_asum(N,A);
    }
  #elif defined(ARMA_USE_BLAS)
    {
    if(has_blas_float_bug<eT>::value)
      {
      out_val = op_norm::vec_norm_1_direct_mem(N,A);
      }
    else
      {
      arma_debug_print("blas::asum()");
      out_val = blas::asum(N,A);
      }
    }
  #else
    {
    out_val = op_norm::vec_norm_1_direct_mem(N,A);
    }
  #endif
  
  return (out_val <= eT(0)) ? eT(0) : out_val;
  }



template<typename eT>
inline
eT
op_norm::vec_norm_1_direct_std(const Mat<eT>& X, const typename arma_fp16_real_only<eT>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);

  const uword N = X.n_elem;
  const eT*   A = X.memptr();

  // fp16 support must be direct non-BLAS
  eT out_val = op_norm::vec_norm_1_direct_mem(N,A);
  
  return (out_val <= eT(0)) ? eT(0) : out_val;
  }



template<typename eT>
inline
eT
op_norm::vec_norm_1_direct_mem(const uword N, const eT* A)
  {
  arma_debug_sigprint();
  
  #if (defined(ARMA_SIMPLE_LOOPS) || defined(__FAST_MATH__))
    {
    eT acc1 = eT(0);
    
    if(memory::is_aligned(A))
      {
      memory::mark_as_aligned(A);
      
      for(uword i=0; i<N; ++i)  { acc1 += std::abs(A[i]); }
      }
    else
      {
      for(uword i=0; i<N; ++i)  { acc1 += std::abs(A[i]); }
      }
    
    return acc1;
    }
  #else
    {
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    
    uword j;
    
    for(j=1; j<N; j+=2)
      {
      const eT tmp_i = (*A);  A++;
      const eT tmp_j = (*A);  A++;
      
      acc1 += std::abs(tmp_i);
      acc2 += std::abs(tmp_j);
      }
    
    if((j-1) < N)
      {
      acc1 += std::abs(*A);
      }
    
    return (acc1 + acc2);
    }
  #endif
  }



template<typename T1>
inline
typename T1::elem_type
op_norm::vec_norm_2(const Proxy<T1>& P, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  constexpr bool use_direct_mem = (is_Mat<typename Proxy<T1>::stored_type>::value) || (is_subview_col<typename Proxy<T1>::stored_type>::value) || (arma_config::openmp && Proxy<T1>::use_mp);
  
  if(use_direct_mem)
    {
    const quasi_unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    
    return op_norm::vec_norm_2_direct_std(tmp.M);
    }
  
  typedef typename T1::elem_type eT;
  
  eT acc = eT(0);
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    
    uword i,j;
    
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      const eT tmp_i = A[i];
      const eT tmp_j = A[j];
      
      acc1 += (tmp_i * tmp_i);
      acc2 += (tmp_j * tmp_j);
      }
    
    if(i < N)
      {
      const eT tmp_i = A[i];
      
      acc1 += (tmp_i * tmp_i);
      }
    
    acc = (acc1 + acc2);
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows == 1)
      {
      for(uword col=0; col<n_cols; ++col)
        {
        const eT tmp = P.at(0,col);
        
        acc += (tmp * tmp);
        }
      }
    else
      {
      eT acc1 = eT(0);
      eT acc2 = eT(0);
      
      for(uword col=0; col<n_cols; ++col)
        {
        uword i,j;
        for(i=0, j=1; j<n_rows; i+=2, j+=2)
          {
          const eT tmp_i = P.at(i,col);
          const eT tmp_j = P.at(j,col);
          
          acc1 += (tmp_i * tmp_i);
          acc2 += (tmp_j * tmp_j);
          }
        
        if(i < n_rows)
          {
          const eT tmp_i = P.at(i,col);
          
          acc1 += (tmp_i * tmp_i);
          }
        }
      
      acc = (acc1 + acc2);
      }
    }
  
  const eT sqrt_acc = std::sqrt(acc);
  
  if( (sqrt_acc != eT(0)) && arma_isfinite(sqrt_acc) )
    {
    return sqrt_acc;
    }
  else
    {
    arma_debug_print("detected possible underflow or overflow");
    
    const quasi_unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    
    return op_norm::vec_norm_2_direct_robust(tmp.M);
    }
  }



template<typename T1>
inline
typename T1::pod_type
op_norm::vec_norm_2(const Proxy<T1>& P, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  T acc = T(0);
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    for(uword i=0; i<N; ++i)
      {
      const std::complex<T>& X = A[i];
      
      const T a = X.real();
      const T b = X.imag();
      
      acc += (a*a) + (b*b);
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows == 1)
      {
      for(uword col=0; col<n_cols; ++col)
        {
        const std::complex<T>& X = P.at(0,col);
        
        const T a = X.real();
        const T b = X.imag();
        
        acc += (a*a) + (b*b);
        }
      }
    else
      {
      for(uword col=0; col<n_cols; ++col)
      for(uword row=0; row<n_rows; ++row)
        {
        const std::complex<T>& X = P.at(row,col);
        
        const T a = X.real();
        const T b = X.imag();
        
        acc += (a*a) + (b*b);
        }
      }
    }
  
  const T sqrt_acc = std::sqrt(acc);
  
  if( (sqrt_acc != T(0)) && arma_isfinite(sqrt_acc) )
    {
    return sqrt_acc;
    }
  else
    {
    arma_debug_print("detected possible underflow or overflow");
    
    typedef typename promote_type<eT, float>::result acc_eT;
    
    typedef typename get_pod_type<acc_eT>::result acc_T;
    
    arma_debug_type_print<eT>("eT");
    arma_debug_type_print<acc_eT>("acc_eT");
    
    const quasi_unwrap<typename Proxy<T1>::stored_type> R(P.Q);
    
    const uword N     = R.M.n_elem;
    const eT*   R_mem = R.M.memptr();
    
    acc_T max_val = priv::most_neg<acc_T>();
    
    for(uword i=0; i<N; ++i)
      {
      const acc_T val_i = std::abs(acc_eT(R_mem[i]));
      
      if(val_i > max_val)  { max_val = val_i; }
      }
    
    if(max_val == acc_T(0))  { return T(0); }
    
    acc_T alt_acc = acc_T(0);
    
    for(uword i=0; i<N; ++i)
      {
      const acc_T val_i = std::abs(acc_eT(R_mem[i])) / max_val;
      
      alt_acc += val_i * val_i;
      }
    
    return T( std::sqrt(alt_acc) * max_val );
    }
  }



template<typename eT>
inline
eT
op_norm::vec_norm_2_direct_std(const Mat<eT>& X, const typename arma_blas_real_only<eT>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = X.n_elem;
  const eT*   A = X.memptr();
  
  eT out_val = eT(0);
  
  #if defined(ARMA_USE_ATLAS)
    {
    arma_debug_print("atlas::cblas_nrm2()");
    out_val = atlas::cblas_nrm2(N,A);
    }
  #elif defined(ARMA_USE_BLAS)
    {
    if(has_blas_float_bug<eT>::value)
      {
      out_val = op_norm::vec_norm_2_direct_mem(N,A);
      }
    else
      {
      arma_debug_print("blas::nrm2()");
      out_val = blas::nrm2(N,A);
      }
    }
  #else
    {
    out_val = op_norm::vec_norm_2_direct_mem(N,A);
    }
  #endif
  
  if( (out_val != eT(0)) && arma_isfinite(out_val) )
    {
    return (out_val < eT(0)) ? eT(0) : out_val;
    }
  else
    {
    arma_debug_print("detected possible underflow or overflow");
    
    return op_norm::vec_norm_2_direct_robust(X);
    }
  }



template<typename eT>
inline
eT
op_norm::vec_norm_2_direct_std(const Mat<eT>& X, const typename arma_fp16_real_only<eT>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = X.n_elem;
  const eT*   A = X.memptr();
  
  // fp16 support must be non-BLAS
  eT out_val = op_norm::vec_norm_2_direct_mem(N,A);
  
  if( (out_val != eT(0)) && arma_isfinite(out_val) )
    {
    return (out_val < eT(0)) ? eT(0) : out_val;
    }
  else
    {
    arma_debug_print("detected possible underflow or overflow");
    
    return op_norm::vec_norm_2_direct_robust(X);
    }
  }



template<typename eT>
inline
eT
op_norm::vec_norm_2_direct_mem(const uword N, const eT* A)
  {
  arma_debug_sigprint();
  
  eT acc = eT(0);
  
  #if (defined(ARMA_SIMPLE_LOOPS) || defined(__FAST_MATH__))
    {
    eT acc1 = eT(0);
    
    if(memory::is_aligned(A))
      {
      memory::mark_as_aligned(A);
      
      for(uword i=0; i<N; ++i)  { const eT tmp_i = A[i];  acc1 += (tmp_i * tmp_i); }
      }
    else
      {
      for(uword i=0; i<N; ++i)  { const eT tmp_i = A[i];  acc1 += (tmp_i * tmp_i); }
      }
    
    acc = acc1;
    }
  #else
    {
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    
    uword j;
    
    for(j=1; j<N; j+=2)
      {
      const eT tmp_i = (*A);  A++;
      const eT tmp_j = (*A);  A++;
      
      acc1 += (tmp_i * tmp_i);
      acc2 += (tmp_j * tmp_j);
      }
    
    if((j-1) < N)
      {
      const eT tmp_i = (*A);
      
      acc1 += (tmp_i * tmp_i);
      }
    
    acc = acc1 + acc2;
    }
  #endif
  
  return std::sqrt(acc);
  }



template<typename eT>
inline
eT
op_norm::vec_norm_2_direct_robust(const Mat<eT>& X)
  {
  arma_debug_sigprint();
  
  typedef typename promote_type<eT, float>::result acc_eT;
  
  arma_debug_type_print<eT>("eT");
  arma_debug_type_print<acc_eT>("acc_eT");
  
  const uword N = X.n_elem;
  
  if(is_fp16<eT>::yes)
    {
    // try straightforward type promotion before the default slow algorithm
    
    const eT* X_mem = X.memptr();
    
    acc_eT acc1 = acc_eT(0);
    acc_eT acc2 = acc_eT(0);
    
    uword j;
    
    for(j=1; j<N; j+=2)
      {
      const acc_eT tmp_i = acc_eT(*X_mem);  X_mem++;
      const acc_eT tmp_j = acc_eT(*X_mem);  X_mem++;
      
      acc1 += (tmp_i * tmp_i);
      acc2 += (tmp_j * tmp_j);
      }
    
    if((j-1) < N)
      {
      const acc_eT tmp_i = acc_eT(*X_mem);
      
      acc1 += (tmp_i * tmp_i);
      }
    
    const acc_eT sqrt_acc = std::sqrt(acc1 + acc2);
    
    if( (sqrt_acc != acc_eT(0)) && arma_isfinite(sqrt_acc) )  { return eT(sqrt_acc); }
    }
  
  const eT* A = X.memptr();
  
  acc_eT max_val = priv::most_neg<acc_eT>();
  
  uword j;
  
  for(j=1; j<N; j+=2)
    {
    acc_eT val_i = acc_eT(*A);  A++;
    acc_eT val_j = acc_eT(*A);  A++;
    
    val_i = std::abs(val_i);
    val_j = std::abs(val_j);
    
    if(val_i > max_val)  { max_val = val_i; }
    if(val_j > max_val)  { max_val = val_j; }
    }
  
  if((j-1) < N)
    {
    const acc_eT val_i = std::abs(acc_eT(*A));
    
    if(val_i > max_val)  { max_val = val_i; }
    }
  
  if(max_val == acc_eT(0))  { return eT(0); }
  
  const eT* B = X.memptr();
  
  acc_eT acc1 = acc_eT(0);
  acc_eT acc2 = acc_eT(0);
  
  for(j=1; j<N; j+=2)
    {
    acc_eT val_i = acc_eT(*B);  B++;
    acc_eT val_j = acc_eT(*B);  B++;
    
    val_i /= max_val;
    val_j /= max_val;
    
    acc1 += (val_i * val_i);
    acc2 += (val_j * val_j);
    }
  
  if((j-1) < N)
    {
    const acc_eT val_i = acc_eT(*B) / max_val;
    
    acc1 += (val_i * val_i);
    }
  
  const acc_eT out_val = std::sqrt(acc1 + acc2) * max_val;
  
  return (out_val <= acc_eT(0)) ? eT(0) : eT(out_val);
  }



template<typename T1>
inline
typename T1::pod_type
op_norm::vec_norm_k(const Proxy<T1>& P, const int k)
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    for(uword i=0; i<N; ++i)
      {
      acc += arma_pow(std::abs(A[i]), k);
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows != 1)
      {
      for(uword col=0; col < n_cols; ++col)
      for(uword row=0; row < n_rows; ++row)
        {
        acc += arma_pow(std::abs(P.at(row,col)), k);
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
        {
        acc += arma_pow(std::abs(P.at(0,col)), k);
        }
      }
    }
  
  return std::pow(acc, T(1)/T(k));
  }



template<typename T1>
inline
typename T1::pod_type
op_norm::vec_norm_max(const Proxy<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const uword N = P.get_n_elem();
  
  T max_val = (N != 1) ? priv::most_neg<T>() : std::abs(P[0]);
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    uword i,j;
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      const T tmp_i = std::abs(A[i]);
      const T tmp_j = std::abs(A[j]);
      
      if(max_val < tmp_i) { max_val = tmp_i; }
      if(max_val < tmp_j) { max_val = tmp_j; }
      }
    
    if(i < N)
      {
      const T tmp_i = std::abs(A[i]);
      
      if(max_val < tmp_i) { max_val = tmp_i; }
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows != 1)
      {
      for(uword col=0; col < n_cols; ++col)
      for(uword row=0; row < n_rows; ++row)
        {
        const T tmp = std::abs(P.at(row,col));
        
        if(max_val < tmp) { max_val = tmp; }
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
        {
        const T tmp = std::abs(P.at(0,col));
        
        if(max_val < tmp) { max_val = tmp; }
        }
      }
    }
  
  return max_val;
  }



template<typename T1>
inline
typename T1::pod_type
op_norm::vec_norm_min(const Proxy<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const uword N = P.get_n_elem();
  
  T min_val = (N != 1) ? priv::most_pos<T>() : std::abs(P[0]);
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    uword i,j;
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      const T tmp_i = std::abs(A[i]);
      const T tmp_j = std::abs(A[j]);
      
      if(min_val > tmp_i) { min_val = tmp_i; }
      if(min_val > tmp_j) { min_val = tmp_j; }
      }
    
    if(i < N)
      {
      const T tmp_i = std::abs(A[i]);
      
      if(min_val > tmp_i) { min_val = tmp_i; }
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows != 1)
      {
      for(uword col=0; col < n_cols; ++col)
      for(uword row=0; row < n_rows; ++row)
        {
        const T tmp = std::abs(P.at(row,col));
        
        if(min_val > tmp) { min_val = tmp; }
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
        {
        const T tmp = std::abs(P.at(0,col));
        
        if(min_val > tmp) { min_val = tmp; }
        }
      }
    }
  
  return min_val;
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
op_norm::mat_norm_1(const Mat<eT>& X)
  {
  arma_debug_sigprint();
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(X), 0), 1) );
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
op_norm::mat_norm_2(const Mat<eT>& X, const typename arma_blas_real_or_cx_only<eT>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename get_pod_type<eT>::result T;
  
  if(X.internal_has_nonfinite())  { arma_warn(1, "norm(): given matrix has non-finite elements"); }
  
  Col<T> S;
  
  arma::svd(S, X);
  
  const T out_val = (S.n_elem > 0) ? S[0] : T(0);
  
  return (out_val <= T(0)) ? T(0) : out_val;
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
op_norm::mat_norm_2(const Mat<eT>& X, const typename arma_fp16_real_or_cx_only<eT>::result* junk)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename get_pod_type<eT>::result T;
  
  typedef typename promote_type<eT, float>::result promoted_eT;
  
  arma_debug_type_print<eT>("eT");
  arma_debug_type_print<promoted_eT>("promoted_eT");
  
  const Mat<promoted_eT> XX = conv_to< Mat<promoted_eT> >::from(X);
  
  return T(op_norm::mat_norm_2(XX));
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
op_norm::mat_norm_inf(const Mat<eT>& X)
  {
  arma_debug_sigprint();
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(X), 1), 0) );
  }



//! @}
