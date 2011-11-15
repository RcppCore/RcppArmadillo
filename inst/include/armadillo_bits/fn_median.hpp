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


//! \addtogroup fn_median
//! @{


template<typename T1>
arma_inline
const Op<T1, op_median>
median(const Base<typename T1::elem_type,T1>& X, const uword dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_median>(X.get_ref(), dim, 0);
  }



//! Immediate 'find the median value of a row vector' operation
template<typename eT>
inline
arma_warn_unused
eT
median(const Row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_elem = A.n_elem;
  
  arma_debug_check( (A_n_elem == 0), "median(): given object has no elements" );
  
  return op_median::direct_median(A.mem, A_n_elem);
  }



//! Immediate 'find the median value of a column vector' operation
template<typename eT>
inline
arma_warn_unused
eT
median(const Col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_elem = A.n_elem;
  
  arma_debug_check( (A_n_elem == 0), "median(): given object has no elements" );
  
  return op_median::direct_median(A.mem, A_n_elem);
  }



//! Immediate 'find the median value of a row vector' operation (complex number version)
template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const Row< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_elem = A.n_elem;
  
  arma_debug_check( (A_n_elem == 0), "median(): given object has no elements" );
  
  uword index1;
  uword index2;
  op_median::direct_cx_median_index(index1, index2, A.mem, A_n_elem);
  
  return (index1 == index2) ? A.mem[index1] : op_median::robust_mean( A.mem[index1], A.mem[index2] );
  }



//! Immediate 'find the median value of a column vector' operation (complex number version)
template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const Col< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_elem = A.n_elem;
  
  arma_debug_check( (A_n_elem == 0), "median(): given object has no elements" );
  
  uword index1;
  uword index2;
  op_median::direct_cx_median_index(index1, index2, A.mem, A_n_elem);
  
  return (index1 == index2) ? A.mem[index1] : op_median::robust_mean( A.mem[index1], A.mem[index2] );
  }



//! find the median value of a subview_row
template<typename eT>
inline
arma_warn_unused
eT
median(const subview_row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given object has no elements" );
  
  return op_median::direct_median(A);
  }



//! find the median value of a subview_col
template<typename eT>
inline
arma_warn_unused
eT
median(const subview_col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given object has no elements" );
  
  return op_median::direct_median(A.colptr(0), A.n_rows);
  }



//! find the median value of a subview_row (complex number version)
template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const subview_row< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given object has no elements" );
  
  uword index1;
  uword index2;
  op_median::direct_cx_median_index(index1, index2, A);
  
  return (index1 == index2) ? A[index1] : op_median::robust_mean(A[index1], A[index2]);
  }



//! find the median value of a subview_col (complex number version)
template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const subview_col< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given object has no elements" );
  
  uword index1;
  uword index2;
  op_median::direct_cx_median_index(index1, index2, A);
  
  return (index1 == index2) ? A[index1] : op_median::robust_mean(A[index1], A[index2]);
  }



template<typename eT>
inline
arma_warn_unused
eT
median(const diagview<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given object has no elements" );
  
  return op_median::direct_median(A);
  }



template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const diagview< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given object has no elements" );
  
  uword index1;
  uword index2;
  op_median::direct_cx_median_index(index1, index2, A);
  
  return (index1 == index2) ? A[index1] : op_median::robust_mean(A[index1], A[index2]);
  }



template<typename eT, typename T1>
inline
arma_warn_unused
eT
median(const subview_elem1<eT,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  const Col<eT> X(A);
  
  return median(X);
  }



//! @}
