// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
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
median(const Base<typename T1::elem_type,T1>& X, const u32 dim = 0)
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
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  return op_median::direct_median(A.mem, A.n_elem);
  }



//! Immediate 'find the median value of a column vector' operation
template<typename eT>
inline
arma_warn_unused
eT
median(const Col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  return op_median::direct_median(A.mem, A.n_elem);
  }



//! Immediate 'find the median value of a row vector' operation (complex number version)
template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const Row< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  u32 index1;
  u32 index2;
  op_median::direct_cx_median_index(index1, index2, A.mem, A.n_elem);
  
  return (A.mem[index1] + A.mem[index2]) / T(2);
  }



//! Immediate 'find the median value of a column vector' operation (complex number version)
template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const Col< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  u32 index1;
  u32 index2;
  op_median::direct_cx_median_index(index1, index2, A.mem, A.n_elem);
  
  return (A.mem[index1] + A.mem[index2]) / T(2);
  }



//! find the median value of a subview_row
template<typename eT>
inline
arma_warn_unused
eT
median(const subview_row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
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
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  return op_median::direct_median(A);
  }



//! find the median value of a subview_row (complex number version)
template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const subview_row< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  u32 index1;
  u32 index2;
  op_median::direct_cx_median_index(index1, index2, A);
  
  return (A[index1] + A[index2]) / T(2);
  }



//! find the median value of a subview_col (complex number version)
template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const subview_col< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  u32 index1;
  u32 index2;
  op_median::direct_cx_median_index(index1, index2, A);
  
  return (A[index1] + A[index2]) / T(2);
  }



template<typename eT>
inline
arma_warn_unused
eT
median(const diagview<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  return op_median::direct_median(A);
  }



template<typename T>
inline
arma_warn_unused
std::complex<T>
median(const diagview< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem == 0), "median(): given vector has no elements" );
  
  u32 index1;
  u32 index2;
  op_median::direct_cx_median_index(index1, index2, A);
  
  return (A[index1] + A[index2]) / T(2);
  }



//! @}
