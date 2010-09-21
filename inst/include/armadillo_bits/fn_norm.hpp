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
arma_inline
typename T1::pod_type
norm_1(const Proxy<T1>& A)
  {
  typedef typename T1::pod_type T;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  T acc = T(0);
  
        ea_type P = A.get_ea();
  const u32     N = A.get_n_elem();
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    acc += std::abs(P[i]);
    acc += std::abs(P[j]);
    }
  
  if(i < N)
    {
    acc += std::abs(P[i]);
    }
  
  return acc;
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_2(const Proxy<T1>& A, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  typedef typename T1::pod_type       T;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  T acc = T(0);
  
        ea_type P = A.get_ea();
  const u32     N = A.get_n_elem();
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    const T tmp_i = P[i];
    const T tmp_j = P[j];
    
    acc += tmp_i * tmp_i;
    acc += tmp_j * tmp_j;
    }
  
  if(i < N)
    {
    const T tmp_i = P[i];
    
    acc += tmp_i * tmp_i;
    }
  
  return std::sqrt(acc);
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_2(const Proxy<T1>& A, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  typedef typename T1::pod_type       T;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  T acc = T(0);
  
        ea_type P = A.get_ea();
  const u32     N = A.get_n_elem();
  
  for(u32 i=0; i<N; ++i)
    {
    const T tmp = std::abs(P[i]);
    acc += tmp*tmp;
    }
  
  return std::sqrt(acc);
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_k(const Proxy<T1>& A, const int k)
  {
  typedef typename T1::pod_type       T;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  T acc = T(0);
  
        ea_type P = A.get_ea();
  const u32     N = A.get_n_elem();
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    acc += std::pow(std::abs(P[i]), k);
    acc += std::pow(std::abs(P[j]), k);
    }
  
  if(i < N)
    {
    acc += std::pow(std::abs(P[i]), k);
    }
  
  return std::pow(acc, T(1)/T(k));
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_max(const Proxy<T1>& A)
  {
  typedef typename T1::pod_type       T;
  typedef typename Proxy<T1>::ea_type ea_type;
  
        ea_type P = A.get_ea();
  const u32     N = A.get_n_elem();
  
  T max_val = std::abs(P[0]);
  
  u32 i,j;
  
  for(i=1, j=2; j<N; i+=2, j+=2)
    {
    const T tmp_i = std::abs(P[i]);
    const T tmp_j = std::abs(P[j]);
    
    if(max_val < tmp_i) { max_val = tmp_i; }
    if(max_val < tmp_j) { max_val = tmp_j; }
    }
  
  if(i < N)
    {
    const T tmp_i = std::abs(P[i]);
    
    if(max_val < tmp_i) { max_val = tmp_i; }
    }
  
  return max_val;
  }



template<typename T1>
arma_hot
arma_inline
typename T1::pod_type
norm_min(const Proxy<T1>& A)
  {
  typedef typename T1::pod_type       T;
  typedef typename Proxy<T1>::ea_type ea_type;
  
        ea_type P = A.get_ea();
  const u32     N = A.get_n_elem();
  
  T min_val = std::abs(P[0]);
  
  u32 i,j;
  
  for(i=1, j=2; j<N; i+=2, j+=2)
    {
    const T tmp_i = std::abs(P[i]);
    const T tmp_j = std::abs(P[j]);
    
    if(min_val > tmp_i) { min_val = tmp_i; }
    if(min_val > tmp_j) { min_val = tmp_j; }
    }
  
  if(i < N)
    {
    const T tmp_i = std::abs(P[i]);
    
    if(min_val > tmp_i) { min_val = tmp_i; }
    }
  
  return min_val;
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
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_check(    (A.get_n_elem() == 0),                            "norm(): given object has no elements"  );
  arma_debug_check( !( (A.get_n_rows() == 1) || (A.get_n_cols() == 1) ), "norm(): given object must be a vector" );
  
  switch(k)
    {
    case 1:
      return norm_1(A);
      break;
    
    case 2:
      return norm_2(A);
      break;
    
    default:
      {
      arma_debug_check( (k == 0), "norm(): k must be greater than zero"   );
      
      return norm_k(A, k);
      }
    }
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
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_check(    (A.get_n_elem() == 0),                            "norm(): given object has no elements"  );
  arma_debug_check( !( (A.get_n_rows() == 1) || (A.get_n_cols() == 1) ), "norm(): given object must be a vector" );
  
  const char sig = method[0];
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // max norm
    {
    return norm_max(A);
    }
  else
  if(sig == '-')   // min norm
    {
    return norm_min(A);
    }
  else
    {
    arma_stop("norm(): unknown norm type");
    
    return T(0);
    }
  
  }



//! @}
