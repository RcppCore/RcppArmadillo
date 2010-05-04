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
norm_unwrap(const Base<typename T1::elem_type,T1>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;

  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  arma_debug_check(    (k == 0),                             "norm(): k must be greater than zero"   );

  const eT* A_mem = A.memptr();
  const u32 N     = A.n_elem;

  if(k==1)
    {
    T acc = T(0);
    
    for(u32 i=0; i<N; ++i)
      {
      acc += std::abs(A_mem[i]);
      }
    
    return acc;
    }
  else
  if(k==2)
    {
    if(is_complex<eT>::value == false)
      {
      eT acc = eT(0);
      
      for(u32 i=0; i<N; ++i)
        {
        const eT tmp = A_mem[i];
        acc += tmp*tmp;
        }
      
      return std::sqrt(access::tmp_real(acc));
      }
    else
      {
      T acc = T(0);
      
      for(u32 i=0; i<N; ++i)
        {
        const T tmp = std::abs(A_mem[i]);
        acc += tmp*tmp;
        }
      
      return std::sqrt(acc);
      }
    }
  else
    {
    T acc = T(0);
    
    for(u32 i=0; i<N; ++i)
      {
      acc += std::pow(std::abs(A_mem[i]), int(k));
      }
    
    return std::pow(acc, T(1)/T(k));
    }
  
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
norm_unwrap(const Base<typename T1::elem_type,T1>& X, const char* method)
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
    T max_val = std::abs(A_mem[0]);
    
    for(u32 i=1; i<N; ++i)
      {
      const T tmp_val = std::abs(A_mem[i]);
      
      if(tmp_val > max_val)
        {
        max_val = tmp_val; 
        }
      }
    
    return max_val;
    }
  else
  if(sig == '-')   // min norm
    {
    T min_val = std::abs(A_mem[0]);
    
    for(u32 i=1; i<N; ++i)
      {
      const T tmp_val = std::abs(A_mem[i]);
      
      if(tmp_val < min_val)
        {
        min_val = tmp_val; 
        }
      }
    
    return min_val;
    }
  else
    {
    arma_stop("norm(): unknown norm type");
    
    return T(0);
    }
  
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
norm_proxy(const Base<typename T1::elem_type,T1>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  arma_debug_check(    (k == 0),                             "norm(): k must be greater than zero"   );
  
  const u32 N = A.n_elem;
  
  if(k==1)
    {
    T acc = T(0);
    
    for(u32 i=0; i<N; ++i)
      {
      acc += std::abs(A[i]);
      }
    
    return acc;
    }
  else
  if(k==2)
    {
    if(is_complex<eT>::value == false)
      {
      eT acc = eT(0);
      
      for(u32 i=0; i<N; ++i)
        {
        const eT tmp = A[i];
        acc += tmp*tmp;
        }
      
      return std::sqrt(access::tmp_real(acc));
      }
    else
      {
      T acc = T(0);
      
      for(u32 i=0; i<N; ++i)
        {
        const T tmp = std::abs(A[i]);
        acc += tmp*tmp;
        }
      
      return std::sqrt(acc);
      }
    }
  else
    {
    T acc = T(0);
    
    for(u32 i=0; i<N; ++i)
      {
      acc += std::pow(std::abs(A[i]), int(k));
      }
    
    return std::pow(acc, T(1)/T(k));
    }
  
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
norm_proxy(const Base<typename T1::elem_type,T1>& X, const char* method)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  
  const u32 N = A.n_elem;
  
  const char sig = method[0];
  
  if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // max norm
    {
    T max_val = std::abs(A[0]);
    
    for(u32 i=1; i<N; ++i)
      {
      const T tmp_val = std::abs(A[i]);
      
      if(tmp_val > max_val)
        {
        max_val = tmp_val; 
        }
      }
    
    return max_val;
    }
  else
  if(sig == '-')   // min norm
    {
    T min_val = std::abs(A[0]);
    
    for(u32 i=1; i<N; ++i)
      {
      const T tmp_val = std::abs(A[i]);
      
      if(tmp_val < min_val)
        {
        min_val = tmp_val; 
        }
      }
    
    return min_val;
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
norm(const Base<typename T1::elem_type,T1>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  if(is_Mat<T1>::value == true)
    {
    return norm_unwrap(X, k);
    }
  else
    {
    return norm_proxy(X, k);
    }
  }



template<typename T1>
arma_inline
arma_warn_unused
typename T1::pod_type
norm(const Base<typename T1::elem_type,T1>& X, const char* method)
  {
  arma_extra_debug_sigprint();
  
  if(is_Mat<T1>::value == true)
    {
    return norm_unwrap(X, method);
    }
  else
    {
    return norm_proxy(X, method);
    }
  }



//! @}
