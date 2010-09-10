// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_norm
//! @{



template<typename T1>
arma_hot
inline
typename T1::pod_type
norm_unwrap(const Base<typename T1::elem_type, T1>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  
  const eT* A_mem = A.memptr();
  const u32 N     = A.n_elem;
  
  switch(k)
    {
    case 1:
      return arrayops::norm_1( A_mem, N );
      break;
    
    case 2:
      return arrayops::norm_2( A_mem, N );
      break;
    
    default:
      arma_debug_check( (k == 0), "norm(): k must be greater than zero" );
      return arrayops::norm_k( A_mem, N, int(k) );
    }
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
norm_unwrap(const Base<typename T1::elem_type, T1>& X, const char* method)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;

  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  
  const eT* A_mem = A.memptr();
  const u32 N     = A.n_elem;
  
  const char sig = method[0];
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // max norm
    {
    return arrayops::norm_max(A_mem, N);
    }
  else
  if(sig == '-')   // min norm
    {
    return arrayops::norm_min(A_mem, N);
    }
  else
    {
    arma_stop("norm(): unknown norm type");
    
    return T(0);
    }
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_1_proxy(const Proxy<T1>& A)
  {
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  const u32 N = A.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    acc += std::abs(A[i]);
    acc += std::abs(A[j]);
    }
  
  if(i < N)
    {
    acc += std::abs(A[i]);
    }
  
  return acc;
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_2_proxy(const Proxy<T1>& A, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  const u32 N = A.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    const T tmp_i = A[i];
    const T tmp_j = A[j];
    
    acc += tmp_i * tmp_i;
    acc += tmp_j * tmp_j;
    }
  
  if(i < N)
    {
    const T tmp_i = A[i];
    
    acc += tmp_i * tmp_i;
    }
  
  return std::sqrt(acc);
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_2_proxy(const Proxy<T1>& A, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  const u32 N = A.n_elem;
  
  for(u32 i=0; i<N; ++i)
    {
    const T tmp = std::abs(A[i]);
    acc += tmp*tmp;
    }
  
  return std::sqrt(acc);
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_k_proxy(const Proxy<T1>& A, const int k)
  {
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  const u32 N = A.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    acc += std::pow(std::abs(A[i]), k);
    acc += std::pow(std::abs(A[j]), k);
    }
  
  if(i < N)
    {
    acc += std::pow(std::abs(A[i]), k);
    }
  
  return std::pow(acc, T(1)/T(k));
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_max_proxy(const Proxy<T1>& A)
  {
  typedef typename T1::pod_type T;
  
  const u32 N = A.n_elem;
  
  T max_val = std::abs(A[0]);
  
  u32 i,j;
  
  for(i=1, j=2; j<N; i+=2, j+=2)
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
  
  return max_val;
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_min_proxy(const Proxy<T1>& A)
  {
  typedef typename T1::pod_type T;
  
  const u32 N = A.n_elem;
  
  T min_val = std::abs(A[0]);
  
  u32 i,j;
  
  for(i=1, j=2; j<N; i+=2, j+=2)
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
  
  return min_val;
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
norm_proxy(const Base<typename T1::elem_type, T1>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  
  switch(k)
    {
    case 1:
      return norm_1_proxy(A);
      break;
    
    case 2:
      return norm_2_proxy(A);
      break;
    
    default:
      {
      arma_debug_check( (k == 0), "norm(): k must be greater than zero"   );
      
      return norm_k_proxy(A, k);
      }
    }
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
norm_proxy(const Base<typename T1::elem_type, T1>& X, const char* method)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  
  const char sig = method[0];
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // max norm
    {
    return norm_max_proxy(A);
    }
  else
  if(sig == '-')   // min norm
    {
    return norm_min_proxy(A);
    }
  else
    {
    arma_stop("norm(): unknown norm type");
    
    return T(0);
    }
  
  }



template<typename T1>
arma_inline
arma_warn_unused
typename T1::pod_type
norm
  (
  const Base<typename T1::elem_type,T1>& X,
  const u32 k,
  const typename arma_float_or_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return (is_Mat<T1>::value == true) ? norm_unwrap(X, k) : norm_proxy(X, k);
  }



template<typename T1>
arma_inline
arma_warn_unused
typename T1::pod_type
norm
  (
  const Base<typename T1::elem_type,T1>& X,
  const char* method,
  const typename arma_float_or_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return (is_Mat<T1>::value == true) ? norm_unwrap(X, method) : norm_proxy(X, method);
  }



//! @}
