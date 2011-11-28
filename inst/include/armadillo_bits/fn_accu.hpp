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


//! \addtogroup fn_accu
//! @{



//! accumulate the elements of a matrix
template<typename T1>
arma_hot
inline
typename T1::elem_type
accu(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> A(X.get_ref());
  
  if(Proxy<T1>::prefer_at_accessor == false)
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
    const uword n_rows = A.get_n_rows();
    const uword n_cols = A.get_n_cols();
    
    eT val = eT(0);
    
    for(uword col=0; col<n_cols; ++col)
      {
      for(uword row=0; row<n_rows; ++row)
        {
        val += A.at(row,col);
        }
      }
    
    return val;
    }
  }



//! explicit handling of Hamming norm (also known as zero norm)
template<typename T1>
arma_inline
arma_warn_unused
uword
accu(const mtOp<uword,T1,op_rel_noteq>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  const uword n_elem = A.get_n_elem();
  const eT    val    = X.aux;
  
  uword n_nonzero = 0;
  for(uword i=0; i<n_elem; ++i)
    {
    if(A[i] != val)
      {
      ++n_nonzero;
      }
    }
  
  return n_nonzero;
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
      {
      for(uword col=0; col<n_cols; ++col)
        {
        for(uword row=0; row<n_rows; ++row)
          {
          val += A.at(row,col,slice);
          }
        }
      }
    
    return val;
    }
  }



//! accumulate the elements of a diagview
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();  
  
  const uword n_elem = X.n_elem;
  
  eT val = eT(0);
  
  for(uword i=0; i<n_elem; ++i)
    {
    val += X[i];
    }
  
  return val;
  }



//! accumulate the elements of a subview (submatrix)
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const subview<eT>& S)
  {
  arma_extra_debug_sigprint();  
  
  const uword S_n_rows = S.n_rows;
  const uword S_n_cols = S.n_cols;
  const uword S_n_elem = S.n_elem;
  
  eT val = eT(0);
  
  if(S_n_elem > 0)
    {
    for(uword col=0; col<S_n_cols; ++col)
      {
      val += arrayops::accumulate( S.colptr(col), S_n_rows );
      }
    }
  
  return val;
  }



//! accumulate the elements of a subview_row
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const subview_row<eT>& S)
  {
  arma_extra_debug_sigprint();  
  
  const Mat<eT>& X = S.m;
  
  const uword n_elem     = S.n_elem;
  const uword row        = S.aux_row1;
  const uword start_col  = S.aux_col1;
  const uword end_col_p1 = start_col + S.n_cols;
  
  eT val = eT(0);
    
  if(n_elem > 0)
    {
    uword i,j;
    
    for(i=start_col, j=start_col+1; j < end_col_p1; i+=2, j+=2)
      {
      val += X.at(row,i);
      val += X.at(row,j);
      }
    
    if(i < end_col_p1)
      {
      val += X.at(row,i);
      }
    }
  
  return val;
  }



//! accumulate the elements of a subview_col
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const subview_col<eT>& S)
  {
  arma_extra_debug_sigprint();
  
  return (S.n_elem > 0) ? arrayops::accumulate( S.colptr(0), S.n_rows ) : eT(0);
  }



//! @}
