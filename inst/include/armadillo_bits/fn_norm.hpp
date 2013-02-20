// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
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
arma_vec_norm_1(const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    T acc1 = T(0);
    T acc2 = T(0);
    
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
    
    acc = acc1 + acc2;
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
      for(uword col=0; col<n_cols; ++col)
        {
        uword i,j;
        
        for(i=0, j=1; j<n_rows; i+=2, j+=2)
          {
          acc += std::abs(P.at(i,col));
          acc += std::abs(P.at(j,col));
          }
        
        if(i < n_rows)
          {
          acc += std::abs(P.at(i,col));
          }
        }
      }
    }
    
  return acc;
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_2
  (
  const Proxy<T1>& P,
  const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    T acc1 = T(0);
    T acc2 = T(0);
    
    uword i,j;
    
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      const T tmp_i = A[i];
      const T tmp_j = A[j];
      
      acc1 += tmp_i * tmp_i;
      acc2 += tmp_j * tmp_j;
      }
    
    if(i < N)
      {
      const T tmp_i = A[i];
      
      acc1 += tmp_i * tmp_i;
      }
    
    acc = acc1 + acc2;
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows == 1)
      {
      for(uword col=0; col<n_cols; ++col)
        {
        const T tmp = P.at(0,col);
        
        acc += tmp * tmp;
        }
      }
    else
      {
      for(uword col=0; col<n_cols; ++col)
        {
        uword i,j;
        for(i=0, j=1; j<n_rows; i+=2, j+=2)
          {
          const T tmp_i = P.at(i,col);
          const T tmp_j = P.at(j,col);
          
          acc += tmp_i * tmp_i;
          acc += tmp_j * tmp_j;
          }
        
        if(i < n_rows)
          {
          const T tmp_i = P.at(i,col);
          
          acc += tmp_i * tmp_i;
          }
        }
      }
    }
  
  return std::sqrt(acc);
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_2
  (
  const Proxy<T1>& P,
  const typename arma_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    for(uword i=0; i<N; ++i)
      {
      const T tmp = std::abs(A[i]);
      acc += tmp*tmp;
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
        const T tmp = std::abs(P.at(0,col));
        acc += tmp*tmp;
        }
      }
    else
      {
      for(uword col=0; col<n_cols; ++col)
      for(uword row=0; row<n_rows; ++row)
        {
        const T tmp = std::abs(P.at(row,col));
        acc += tmp*tmp;
        }
      }
    }
  
  return std::sqrt(acc);
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_k
  (
  const Proxy<T1>& P,
  const int k
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  T acc = T(0);
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    const uword N = P.get_n_elem();
    
    uword i,j;
    
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      acc += std::pow(std::abs(A[i]), k);
      acc += std::pow(std::abs(A[j]), k);
      }
    
    if(i < N)
      {
      acc += std::pow(std::abs(A[i]), k);
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
        acc += std::pow(std::abs(P.at(row,col)), k);
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
        {
        acc += std::pow(std::abs(P.at(0,col)), k);
        }
      }
    }
  
  return std::pow(acc, T(1)/T(k));
  }



template<typename T1>
arma_hot
inline
typename T1::pod_type
arma_vec_norm_max(const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const uword N = P.get_n_elem();
  
  T max_val = (N != 1) ? priv::most_neg<T>() : std::abs(P[0]);
  
  if(Proxy<T1>::prefer_at_accessor == false)
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
arma_hot
inline
typename T1::pod_type
arma_vec_norm_min(const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const uword N = P.get_n_elem();
  
  T min_val = (N != 1) ? priv::most_pos<T>() : std::abs(P[0]);
  
  if(Proxy<T1>::prefer_at_accessor == false)
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



template<typename T1>
inline
typename T1::pod_type
arma_mat_norm_1(const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(P.Q), 0), 1) );
  }



template<typename T1>
inline
typename T1::pod_type
arma_mat_norm_2(const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type   T;
  
  // TODO: is the SVD based approach only valid for square matrices?
  
  Col<T> S;
  svd(S, P.Q);
  
  return (S.n_elem > 0) ? max(S) : T(0);
  }



template<typename T1>
inline
typename T1::pod_type
arma_mat_norm_inf(const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(P.Q), 1), 0) );
  }



template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_type<T1>::value, typename T1::pod_type >::result
norm
  (
  const T1& X,
  const uword k,
  const typename arma_float_or_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> P(X);
  
  if(P.get_n_elem() == 0)
    {
    return T(0);
    }
  
  const bool is_vec = (P.get_n_rows() == 1) || (P.get_n_cols() == 1);
  
  if(is_vec == true)
    {
    switch(k)
      {
      case 1:
        return arma_vec_norm_1(P);
        break;
      
      case 2:
        return arma_vec_norm_2(P);
        break;
      
      default:
        {
        arma_debug_check( (k == 0), "norm(): k must be greater than zero"   );
        return arma_vec_norm_k(P, int(k));
        }
      }
    }
  else
    {
    switch(k)
      {
      case 1:
        return arma_mat_norm_1(P);
        break;
      
      case 2:
        return arma_mat_norm_2(P);
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
typename enable_if2< is_arma_type<T1>::value, typename T1::pod_type >::result
norm
  (
  const T1& X,
  const char* method,
  const typename arma_float_or_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> P(X);
  
  if(P.get_n_elem() == 0)
    {
    return T(0);
    }
  
  const char sig    = method[0];
  const bool is_vec = (P.get_n_rows() == 1) || (P.get_n_cols() == 1);
  
  if(is_vec == true)
    {
    if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // max norm
      {
      return arma_vec_norm_max(P);
      }
    else
    if(sig == '-')   // min norm
      {
      return arma_vec_norm_min(P);
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return arma_vec_norm_2(P);
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
      return arma_mat_norm_inf(P);
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return arma_vec_norm_2(P);
      }
    else
      {
      arma_stop("norm(): unsupported matrix norm type");
      return T(0);
      }
    }
  }



//
// norms for sparse matrices



template<typename T1>
inline
typename T1::pod_type
arma_mat_norm_1(const SpProxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(P.Q), 0), 1) );
  }



// template<typename T1>
// inline
// typename T1::pod_type
// arma_mat_norm_2(const SpProxy<T1>& P)
//   {
//   arma_extra_debug_sigprint();
//   
//   // TODO: norm = sqrt( largest eigenvalue of (A^H)*A ), where ^H is the conjugate transpose
//   // TODO: can use ARPACK or directly implement the Arnoldi iteration
//   // http://math.stackexchange.com/questions/4368/computing-the-largest-eigenvalue-of-a-very-large-sparse-matrix
//   }



template<typename T1>
inline
typename T1::pod_type
arma_mat_norm_inf(const SpProxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(P.Q), 1), 0) );
  }



template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_sparse_type<T1>::value, typename T1::pod_type >::result
norm
  (
  const T1& X,
  const uword k,
  const typename arma_float_or_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const SpProxy<T1> P(X);
  
  if(P.get_n_nonzero() == 0)
    {
    return T(0);
    }
  
  const bool is_vec = (P.get_n_rows() == 1) || (P.get_n_cols() == 1);
  
  if(is_vec == true)
    {
    const unwrap_spmat<typename SpProxy<T1>::stored_type> tmp(P.Q);
    const SpMat<eT>& A = tmp.M;
    
    // create a fake dense vector to allow reuse of code for dense vectors
    Col<eT> fake_vector( access::rwp(A.values), A.n_nonzero, false );
    
    const Proxy< Col<eT> > P_fake_vector(fake_vector);
    
    switch(k)
      {
      case 1:
        return arma_vec_norm_1(P_fake_vector);
        break;
      
      case 2:
        return arma_vec_norm_2(P_fake_vector);
        break;
      
      default:
        {
        arma_debug_check( (k == 0), "norm(): k must be greater than zero"   );
        return arma_vec_norm_k(P_fake_vector, int(k));
        }
      }
    }
  else
    {
    switch(k)
      {
      case 1:
        return arma_mat_norm_1(P);
        break;
      
      // case 2:
      //   return arma_mat_norm_2(P);
      //   break;
      
      default:
        arma_stop("norm(): unsupported or unimplemented norm type for sparse matrices");
        return T(0);
      }
    }
  }



template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_sparse_type<T1>::value, typename T1::pod_type >::result
norm
  (
  const T1& X,
  const char* method,
  const typename arma_float_or_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const SpProxy<T1> P(X);
  
  if(P.get_n_nonzero() == 0)
    {
    return T(0);
    }
  
  
  const unwrap_spmat<typename SpProxy<T1>::stored_type> tmp(P.Q);
  const SpMat<eT>& A = tmp.M;
  
  // create a fake dense vector to allow reuse of code for dense vectors
  Col<eT> fake_vector( access::rwp(A.values), A.n_nonzero, false );
  
  const Proxy< Col<eT> > P_fake_vector(fake_vector);
  
  
  const char sig    = method[0];
  const bool is_vec = (P.get_n_rows() == 1) || (P.get_n_cols() == 1);
  
  if(is_vec == true)
    {
    if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // max norm
      {
      return arma_vec_norm_max(P_fake_vector);
      }
    else
    if(sig == '-')   // min norm
      {
      const T val = arma_vec_norm_min(P_fake_vector);
      
      if( P.get_n_nonzero() < P.get_n_elem() )
        {
        return (std::min)(T(0), val);
        }
      else
        {
        return val;
        }
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return arma_vec_norm_2(P_fake_vector);
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
      return arma_mat_norm_inf(P);
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return arma_vec_norm_2(P_fake_vector);
      }
    else
      {
      arma_stop("norm(): unsupported matrix norm type");
      return T(0);
      }
    }
  }



//! @}
