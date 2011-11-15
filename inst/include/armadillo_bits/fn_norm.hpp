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


//! \addtogroup fn_norm
//! @{



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_1(const Proxy<T1>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;

  T acc = T(0);
    
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type P = A.get_ea();
    
    const uword N = A.get_n_elem();
    
    uword i,j;
    
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      acc += std::abs(P[i]);
      acc += std::abs(P[j]);
      }
    
    if(i < N)
      {
      acc += std::abs(P[i]);
      }
    }
  else
    {
    const uword n_rows = A.get_n_rows();
    const uword n_cols = A.get_n_cols();
    
    for(uword col=0; col<n_cols; ++col)
      {
      uword i,j;
      
      for(i=0, j=1; j<n_rows; i+=2, j+=2)
        {
        acc += std::abs(A.at(i,col));
        acc += std::abs(A.at(j,col));
        }
      
      if(i < n_rows)
        {
        acc += std::abs(A.at(i,col));
        }
      }
    }
    
  return acc;
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_2(const Proxy<T1>& A, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type P = A.get_ea();
    
    const uword N = A.get_n_elem();
    
    uword i,j;
    
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
    }
  else
    {
    const uword n_rows = A.get_n_rows();
    const uword n_cols = A.get_n_cols();
    
    for(uword col=0; col<n_cols; ++col)
      {
      uword i,j;
      
      for(i=0, j=1; j<n_rows; i+=2, j+=2)
        {
        const T tmp_i = A.at(i,col);
        const T tmp_j = A.at(j,col);
        
        acc += tmp_i * tmp_i;
        acc += tmp_j * tmp_j;
        }
      
      if(i < n_rows)
        {
        const T tmp_i = A.at(i,col);
        
        acc += tmp_i * tmp_i;
        }
      }
    }
  
  return std::sqrt(acc);
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_2(const Proxy<T1>& A, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type P = A.get_ea();
    
    const uword N = A.get_n_elem();
    
    for(uword i=0; i<N; ++i)
      {
      const T tmp = std::abs(P[i]);
      acc += tmp*tmp;
      }
    }
  else
    {
    const uword n_rows = A.get_n_rows();
    const uword n_cols = A.get_n_cols();
    
    for(uword col=0; col<n_cols; ++col)
    for(uword row=0; row<n_rows; ++row)
      {
      const T tmp = std::abs(A.at(row,col));
      acc += tmp*tmp;
      }
    }
  
  return std::sqrt(acc);
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_k(const Proxy<T1>& A, const int k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type P = A.get_ea();
    
    const uword N = A.get_n_elem();
    
    uword i,j;
    
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      acc += std::pow(std::abs(P[i]), k);
      acc += std::pow(std::abs(P[j]), k);
      }
    
    if(i < N)
      {
      acc += std::pow(std::abs(P[i]), k);
      }
    }
  else
    {
    const uword n_rows = A.get_n_rows();
    const uword n_cols = A.get_n_cols();
    
    for(uword col=0; col<n_cols; ++col)
    for(uword row=0; row<n_rows; ++row)
      {
      acc += std::pow(std::abs(A.at(row,col)), k);
      }
    }
  
  return std::pow(acc, T(1)/T(k));
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_max(const Proxy<T1>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type       T;
  typedef typename Proxy<T1>::ea_type ea_type;
  
        ea_type P = A.get_ea();
  const uword     N = A.get_n_elem();
  
  T max_val = (N != 1) ? priv::most_neg<T>() : std::abs(P[0]);
  
  uword i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
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
inline
typename T1::pod_type
arma_vec_norm_min(const Proxy<T1>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type       T;
  typedef typename Proxy<T1>::ea_type ea_type;
  
        ea_type P = A.get_ea();
  const uword     N = A.get_n_elem();
  
  T min_val = (N != 1) ? priv::most_pos<T>() : std::abs(P[0]);
  
  uword i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
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
inline
typename T1::pod_type
arma_mat_norm_1(const Proxy<T1>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const unwrap<typename Proxy<T1>::stored_type> tmp(A.Q);
  const Mat<eT>& X = tmp.M;
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(X)), 1) );
  }



template<typename T1>
inline
typename T1::pod_type
arma_mat_norm_2(const Proxy<T1>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const unwrap<typename Proxy<T1>::stored_type> tmp(A.Q);
  const Mat<eT>& X = tmp.M;
  
  Col<T> S;
  svd(S, X);
  
  return (S.n_elem > 0) ? max(S) : T(0);
  }



template<typename T1>
inline
typename T1::pod_type
arma_mat_norm_inf(const Proxy<T1>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const unwrap<typename Proxy<T1>::stored_type> tmp(A.Q);
  const Mat<eT>& X = tmp.M;
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(X),1) ) );
  }



template<typename T1>
inline
arma_warn_unused
typename T1::pod_type
norm
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword k,
  const typename arma_float_or_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> A(X.get_ref());
  
  if(A.get_n_elem() == 0)
    {
    return T(0);
    }
  
  const bool is_vec = (A.get_n_rows() == 1) || (A.get_n_cols() == 1);
  
  if(is_vec == true)
    {
    switch(k)
      {
      case 1:
        return arma_vec_norm_1(A);
        break;
      
      case 2:
        return arma_vec_norm_2(A);
        break;
      
      default:
        {
        arma_debug_check( (k == 0), "norm(): k must be greater than zero"   );
        return arma_vec_norm_k(A, int(k));
        }
      }
    }
  else
    {
    switch(k)
      {
      case 1:
        return arma_mat_norm_1(A);
        break;
      
      case 2:
        return arma_mat_norm_2(A);
        break;
      
      default:
        arma_stop("norm(): unsupported matrix norm type");
        return T(0);
      }
    }
  }



template<typename T1>
inline
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
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> A(X.get_ref());
  
  if(A.get_n_elem() == 0)
    {
    return T(0);
    }
  
  const char sig    = method[0];
  const bool is_vec = (A.get_n_rows() == 1) || (A.get_n_cols() == 1);
  
  if(is_vec == true)
    {
    if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // max norm
      {
      return arma_vec_norm_max(A);
      }
    else
    if(sig == '-')   // min norm
      {
      return arma_vec_norm_min(A);
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return arma_vec_norm_2(A);
      }
    else
      {
      arma_stop("norm(): unsupported vector norm type");
      return T(0);
      }
    }
  else
    {
    if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // inf norm
      {
      return arma_mat_norm_inf(A);
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return arma_vec_norm_2(A);
      }
    else
      {
      arma_stop("norm(): unsupported matrix norm type");
      return T(0);
      }
    }
  }



//! @}
