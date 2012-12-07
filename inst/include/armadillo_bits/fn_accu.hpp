// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// Copyright (C) 2012 Ryan Curtin
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_accu
//! @{



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_proxy_linear(const Proxy<T1>& P)
  {
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
        ea_type A      = P.get_ea();
  const uword   n_elem = P.get_n_elem();
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  uword i,j;
  for(i=0, j=1; j < n_elem; i+=2, j+=2)
    {
    val1 += A[i];
    val2 += A[j];
    }
  
  if(i < n_elem)
    {
    val1 += A[i];   // equivalent to: val1 += A[n_elem-1];
    }
  
  return (val1 + val2);
  }



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_proxy_at(const Proxy<T1>& P)
  {
  typedef typename T1::elem_type eT;
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  eT val = eT(0);
    
  if(n_rows != 1)
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      val += P.at(row,col);
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
      {
      val += P.at(0,col);
      }
    }
  
  return val;
  }



//! accumulate the elements of a matrix
template<typename T1>
arma_hot
inline
typename enable_if2< is_arma_type<T1>::value, typename T1::elem_type >::result
accu(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> P(X);
  
  return (Proxy<T1>::prefer_at_accessor == false) ? accu_proxy_linear(P) : accu_proxy_at(P);
  }



//! explicit handling of Hamming norm (also known as zero norm)
template<typename T1>
inline
arma_warn_unused
uword
accu(const mtOp<uword,T1,op_rel_noteq>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const eT val = X.aux;
  
  const Proxy<T1> P(X.m);
  
  uword n_nonzero = 0;
    
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
          ea_type A      = P.get_ea();
    const uword   n_elem = P.get_n_elem();
    
    for(uword i=0; i<n_elem; ++i)
      {
      if(A[i] != val) { ++n_nonzero; }
      }
    }
  else
    {
    const uword P_n_cols = P.get_n_cols();
    const uword P_n_rows = P.get_n_rows();
    
    if(P_n_rows == 1)
      {
      for(uword col=0; col < P_n_cols; ++col)
        {
        if(P.at(0,col) != val) { ++n_nonzero; }
        }
      }
    else
      {
      for(uword col=0; col < P_n_cols; ++col)
      for(uword row=0; row < P_n_rows; ++row)
        {
        if(P.at(row,col) != val) { ++n_nonzero; }
        }
      }
    }
  
  return n_nonzero;
  }



//! accumulate the elements of a subview (submatrix)
template<typename eT>
arma_hot
arma_pure
arma_warn_unused
inline
eT
accu(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();  
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  eT val = eT(0);
  
  if(X_n_rows == 1)
    {
    const Mat<eT>& A = X.m;
    
    const uword start_row = X.aux_row1;
    const uword start_col = X.aux_col1;
    
    const uword end_col_p1 = start_col + X_n_cols;
    
    uword i,j;
    for(i=start_col, j=start_col+1; j < end_col_p1; i+=2, j+=2)
      {
      val += A.at(start_row, i);
      val += A.at(start_row, j);
      }
    
    if(i < end_col_p1)
      {
      val += A.at(start_row, i);
      }
    }
  else
  if(X_n_cols == 1)
    {
    val = arrayops::accumulate( X.colptr(0), X_n_rows );
    }
  else
    {
    for(uword col=0; col < X_n_cols; ++col)
      {
      val += arrayops::accumulate( X.colptr(col), X_n_rows );
      }
    }
  
  return val;
  }



template<typename eT>
arma_hot
arma_pure
arma_warn_unused
inline
eT
accu(const subview_col<eT>& X)
  {
  arma_extra_debug_sigprint();  
  
  return arrayops::accumulate( X.colptr(0), X.n_rows );
  }



//! accumulate the elements of a cube
template<typename T1>
arma_hot
arma_warn_unused
inline
typename T1::elem_type
accu(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> A(X.get_ref());
  
  if(ProxyCube<T1>::prefer_at_accessor == false)
    {
          ea_type P      = A.get_ea();
    const uword   n_elem = A.get_n_elem();
    
    eT val1 = eT(0);
    eT val2 = eT(0);
    
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      val1 += P[i];
      val2 += P[j];
      }
    
    if(i < n_elem)
      {
      val1 += P[i];
      }
    
    return val1 + val2;
    }
  else
    {
    const uword n_rows   = A.get_n_rows();
    const uword n_cols   = A.get_n_cols();
    const uword n_slices = A.get_n_slices();
    
    eT val = eT(0);
    
    for(uword slice=0; slice<n_slices; ++slice)
    for(uword col=0; col<n_cols; ++col)
    for(uword row=0; row<n_rows; ++row)
      {
      val += A.at(row,col,slice);
      }
    
    return val;
    }
  }



template<typename T>
arma_inline
arma_warn_unused
const typename arma_scalar_only<T>::result &
accu(const T& x)
  {
  return x;
  }



//! accumulate values in a sparse object
template<typename T1>
arma_hot
inline
arma_warn_unused
typename enable_if2<is_arma_sparse_type<T1>::value, typename T1::elem_type>::result
accu(const T1& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> p(x);
  
  if(SpProxy<T1>::must_use_iterator == false)
    {
    // direct counting
    return arrayops::accumulate(p.get_values(), p.get_n_nonzero());
    }
  else
    {
    typename SpProxy<T1>::const_iterator_type it     = p.begin();
    typename SpProxy<T1>::const_iterator_type it_end = p.end();
    
    eT result = eT(0);
    
    while(it != it_end)
      {
      result += (*it);
      ++it;
      }
    
    return result;
    }
  }



//! @}
