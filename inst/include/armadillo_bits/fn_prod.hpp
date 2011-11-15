// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_prod
//! @{


//! \brief
//! Delayed product of elements of a matrix along a specified dimension (either rows or columns).
//! The result is stored in a dense matrix that has either one column or one row.
//! For dim = 0, find the sum of each column (i.e. traverse across rows)
//! For dim = 1, find the sum of each row (i.e. traverse across columns)
//! The default is dim = 0.
//! NOTE: this function works differently than in Matlab/Octave.

template<typename T1>
arma_inline
const Op<T1, op_prod>
prod(const Base<typename T1::elem_type,T1>& X, const uword dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_prod>(X.get_ref(), dim, 0);
  }



//! \brief
//! Immediate 'product of all values' operation for a row vector
template<typename eT>
inline
arma_warn_unused
eT
prod(const Row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return arrayops::product(X.memptr(), X.n_elem);
  }



//! \brief
//! Immediate 'product of all values' operation for a column vector
template<typename eT>
inline
arma_warn_unused
eT
prod(const Col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return arrayops::product(X.memptr(), X.n_elem);
  }



//! \brief
//! Immediate 'product of all values' operation,
//! invoked, for example, by: prod(prod(A))

template<typename T1>
inline
typename T1::elem_type
prod(const Op<T1, op_prod>& in)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("prod(): two consecutive prod() calls detected");
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  return arrayops::product( X.memptr(), X.n_elem );
  }



template<typename T1>
inline
const Op<Op<T1, op_prod>, op_prod>
prod(const Op<T1, op_prod>& in, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op<Op<T1, op_prod>, op_prod>(in, dim, 0);
  }



//! product of all values of a subview_row
template<typename eT>
inline
arma_warn_unused
eT
prod(const subview_row<eT>& S)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& X = S.m;
  
  const uword n_elem         = S.n_elem;
  const uword row            = S.aux_row1;
  const uword start_col      = S.aux_col1;
  const uword end_col_plus_1 = start_col + S.n_cols;
  
  eT val = eT(1);
  
  if(n_elem > 0)
    {
    for(uword col=start_col; col<end_col_plus_1; ++col)
      {
      val *= X.at(row,col);
      }
    }
  
  return val;
  }



//! product of all values of a subview_col
template<typename eT>
inline
arma_warn_unused
eT
prod(const subview_col<eT>& S)
  {
  arma_extra_debug_sigprint();
  
  return (S.n_elem > 0) ? arrayops::product( S.colptr(0), S.n_rows ) : eT(1);
  }



//! product of all values of a diagview
template<typename eT>
arma_warn_unused
inline
eT
prod(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword X_n_elem = X.n_elem;
  
  eT val = eT(1);
  
  for(uword i=0; i<X_n_elem; ++i)
    {
    val *= X[i];
    }
  
  return val;
  }



template<typename eT, typename T1>
inline
arma_warn_unused
eT
prod(const subview_elem1<eT,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  const Col<eT> X(A);
  
  return prod(X);
  }



//! @}
