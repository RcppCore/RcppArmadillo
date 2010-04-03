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
typename T1::elem_type
norm_unwrap(const Base<typename T1::elem_type,T1>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;

  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  arma_debug_check(    (k == 0),                             "norm(): k must be greater than zero"   );

  const eT* A_mem = A.memptr();
  const u32 N     = A.n_elem;

  if(k==1)
    {
    eT acc = eT(0);
    
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
      
      return std::sqrt(acc);
      }
    else
      {
      eT acc = eT(0);
      
      for(u32 i=0; i<N; ++i)
        {
        acc += std::abs(A_mem[i]);
        }
      
      return std::sqrt(acc);
      }
    }
  else
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<N; ++i)
      {
      acc += std::pow(std::abs(A_mem[i]), int(k));
      }
    
    return std::pow(acc, eT(1)/eT(k));
    }
  
  }



template<typename T1>
arma_hot
inline
typename T1::elem_type
norm_proxy(const Base<typename T1::elem_type,T1>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements"  );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): given object must be a vector" );
  arma_debug_check(    (k == 0),                             "norm(): k must be greater than zero"   );
  
  const u32 N = A.n_elem;
  
  if(k==1)
    {
    eT acc = eT(0);
    
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
      
      return std::sqrt(acc);
      }
    else
      {
      eT acc = eT(0);
      
      for(u32 i=0; i<N; ++i)
        {
        acc += std::abs(A[i]);
        }
      
      return std::sqrt(acc);
      
      }
    }
  else
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<N; ++i)
      {
      acc += std::pow(std::abs(A[i]), int(k));
      }
    
    return std::pow(acc, eT(1)/eT(k));
    }
  
  }



template<typename T1>
arma_inline
arma_warn_unused
typename T1::elem_type
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



//! @}
