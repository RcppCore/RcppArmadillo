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
  
        ea_type P      = A.get_ea();
  const u32     n_elem = A.get_n_elem();
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  u32 i,j;
  
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



//! explicit handling of Hamming norm (also known as zero norm)
template<typename T1>
arma_inline
arma_warn_unused
u32
accu(const mtOp<u32,T1,op_rel_noteq>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  const u32 n_elem = A.get_n_elem();
  const eT  val    = X.aux;
  
  u32 n_nonzero = 0;
  for(u32 i=0; i<n_elem; ++i)
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
  
        ea_type P      = A.get_ea();
  const u32     n_elem = A.get_n_elem();
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  u32 i,j;
  
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



//! accumulate the elements of a diagview
template<typename eT>
arma_pure
arma_warn_unused
inline
eT
accu(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();  
  
  const u32 n_elem = X.n_elem;
  
  eT val = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
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
  
  if(S.n_elem > 0)
    {
    const u32 S_n_rows = S.n_rows;
    const u32 S_n_cols = S.n_cols;
    
    eT val = eT(0);
    
    for(u32 col=0; col<S_n_cols; ++col)
      {
      val += arrayops::accumulate( S.colptr(col), S_n_rows );
      }
    
    return val;
    }
  else
    {
    return eT(0);
    }
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
  
  const u32 row            = S.aux_row1;
  const u32 start_col      = S.aux_col1;
  const u32 end_col_plus_1 = start_col + S.n_cols;
  
  // S.n_cols might be equal to zero,
  // hence the loop below has a "less than" condition
  
  eT val = eT(0);
  
  u32 i,j;
  
  for(i=start_col, j=start_col+1; j < end_col_plus_1; i+=2, j+=2)
    {
    val += X.at(row,i);
    val += X.at(row,j);
    }
  
  if(i < end_col_plus_1)
    {
    val += X.at(row,i);
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
  
  return arrayops::accumulate( S.colptr(0), S.n_rows );
  }



//! @}
