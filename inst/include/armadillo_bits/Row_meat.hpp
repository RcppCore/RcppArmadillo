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


//! \addtogroup Row
//! @{

template<typename eT>
inline
Row<eT>::Row()
  : Mat<eT>()
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
Row<eT>::Row(const u32 in_n_elem)
  : Mat<eT>(1,in_n_elem)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
Row<eT>::Row(const u32 in_n_rows, const u32 in_n_cols)
  : Mat<eT>(in_n_rows, in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
Row<eT>::Row(const char* text)
  : Mat<eT>(text)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }
  


template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(text);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
Row<eT>::Row(const std::string& text)
  : Mat<eT>(text)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }
  


template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(text);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
Row<eT>::Row(const Row<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const Row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  return *this;
  }



template<typename eT>
inline
Row<eT>::Row(const Mat<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator*=(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



//! construct a row vector from a given auxiliary array
template<typename eT>
inline
Row<eT>::Row(eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const bool copy_aux_mem)
  : Mat<eT>(aux_mem, aux_n_rows, aux_n_cols, copy_aux_mem)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



//! construct a row vector from a given auxiliary array
template<typename eT>
inline
Row<eT>::Row(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols)
  : Mat<eT>(aux_mem, aux_n_rows, aux_n_cols)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



//! construct a row vector from a given auxiliary array
template<typename eT>
inline
Row<eT>::Row(eT* aux_mem, const u32 aux_length, const bool copy_aux_mem)
  : Mat<eT>(aux_mem, 1, aux_length, copy_aux_mem)
  {
  arma_extra_debug_sigprint();
  
//   Mat<eT>::set_size(1, aux_length);
//   arma_check( (Mat<eT>::n_elem != aux_length), "Row(): don't know how to handle the given array" );
// 
//   syslib::copy_elem( Mat<eT>::memptr(), aux_mem, Mat<eT>::n_elem );
  }



//! construct a row vector from a given auxiliary array
template<typename eT>
inline
Row<eT>::Row(const eT* aux_mem, const u32 aux_length)
  : Mat<eT>(aux_mem, 1, aux_length)
  {
  arma_extra_debug_sigprint();
  
//   Mat<eT>::set_size(1, aux_length);
//   arma_check( (Mat<eT>::n_elem != aux_length), "Row(): don't know how to handle the given array" );
// 
//   syslib::copy_elem( Mat<eT>::memptr(), aux_mem, Mat<eT>::n_elem );
  }



template<typename eT>
template<typename T1, typename T2>
inline
Row<eT>::Row
  (
  const Base<typename Row<eT>::pod_type, T1>& A,
  const Base<typename Row<eT>::pod_type, T2>& B
  )
  : Mat<eT>(A,B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
Row<eT>::Row(const subview<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator*=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
Row<eT>::Row(const subview_cube<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator*=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



//! construct a row vector from given a diagview
template<typename eT>
inline
Row<eT>::Row(const diagview<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



//! construct a row vector from given a diagview
template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  //std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator*=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
arma_inline
eT&
Row<eT>::col(const u32 col_num)
  {
  arma_debug_check( (col_num >= Mat<eT>::n_cols), "Row::col(): out of bounds" );
  
  return access::rw(Mat<eT>::mem[col_num]);
  }



template<typename eT>
arma_inline
eT
Row<eT>::col(const u32 col_num)
  const
  {
  arma_debug_check( (col_num >= Mat<eT>::n_cols), "Row::col(): out of bounds" );
  
  return Mat<eT>::mem[col_num];
  }



template<typename eT>
arma_inline
subview_row<eT>
Row<eT>::cols(const u32 in_col1, const u32 in_col2)
  {
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "Row::cols(): indices out of bounds or incorrectly used");
  
  return subview_row<eT>(*this, 0, in_col1, in_col2);
  }



template<typename eT>
arma_inline
const subview_row<eT>
Row<eT>::cols(const u32 in_col1, const u32 in_col2)
  const
  {
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "Row::cols(): indices out of bounds or incorrectly used");
  
  return subview_row<eT>(*this, 0, in_col1, in_col2);
  }



template<typename eT>
template<typename T1, typename op_type>
inline
Row<eT>::Row(const Op<T1, op_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Row<eT>&
Row<eT>::operator=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Row<eT>&
Row<eT>::operator*=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
Row<eT>::Row(const eOp<T1, eop_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Row<eT>&
Row<eT>::operator=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Row<eT>&
Row<eT>::operator*=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
Row<eT>::Row(const mtOp<eT, T1, op_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Row<eT>&
Row<eT>::operator=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Row<eT>&
Row<eT>::operator*=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Row<eT>::Row(const Glue<T1, T2, glue_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Row<eT>&
Row<eT>::operator=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Row<eT>&
Row<eT>::operator*=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
Row<eT>::Row(const eGlue<T1, T2, eglue_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Row<eT>&
Row<eT>::operator=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Row<eT>&
Row<eT>::operator*=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Row<eT>::Row(const mtGlue<eT, T1, T2, glue_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Row<eT>&
Row<eT>::operator=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Row<eT>&
Row<eT>::operator*=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
void
Row<eT>::set_size(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::set_size(1,in_n_elem);
  }



template<typename eT>
inline
void
Row<eT>::set_size(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  // min is used in case in_n_rows is zero
  Mat<eT>::set_size( (std::min)( u32(1), in_n_rows), in_n_cols );
  
  arma_debug_check( (in_n_rows > 1), "Row::set_size(): incompatible dimensions" );
  }



template<typename eT>
template<typename eT2>
inline
void
Row<eT>::copy_size(const Mat<eT2>& x)
  {
  arma_extra_debug_sigprint();
  
  // min is used in case x.n_rows is zero
  Mat<eT>::set_size( (std::min)( u32(1), x.n_rows), x.n_cols );
  
  arma_debug_check( (x.n_rows > 1), "Row::copy_size(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Row<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros();
  }



template<typename eT>
inline
void
Row<eT>::zeros(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros(1, in_n_elem);
  }



template<typename eT>
inline
void
Row<eT>::zeros(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  // min is used in case in_n_rows is zero
  Mat<eT>::zeros( (std::min)( u32(1), in_n_rows), in_n_cols );
  
  arma_debug_check( (in_n_rows > 1), "Row<eT>::zeros(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Row<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::ones();
  }



template<typename eT>
inline
void
Row<eT>::ones(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::ones(1, in_n_elem);
  }



template<typename eT>
inline
void
Row<eT>::ones(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  // min is used in case in_n_rows is zero
  Mat<eT>::ones( (std::min)( u32(1), in_n_rows), in_n_cols );
  
  arma_debug_check( (in_n_rows > 1), "Row<eT>::ones(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Row<eT>::load(const std::string name, const file_type type, const bool print_status)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::load(name, type, print_status);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Row<eT>::load(std::istream& is, const file_type type, const bool print_status)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::load(is, type, print_status);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Row<eT>::quiet_load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::quiet_load(name, type);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Row<eT>::quiet_load(std::istream& is, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::quiet_load(is, type);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
typename Row<eT>::row_iterator
Row<eT>::begin_row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return Mat<eT>::memptr();
  }



template<typename eT>
inline
typename Row<eT>::const_row_iterator
Row<eT>::begin_row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return Mat<eT>::memptr();
  }



template<typename eT>
inline
typename Row<eT>::row_iterator
Row<eT>::end_row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "end_row(): index out of bounds");
  
  return Mat<eT>::memptr() + Mat<eT>::n_cols;
  }



template<typename eT>
inline
typename Row<eT>::const_row_iterator
Row<eT>::end_row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "end_row(): index out of bounds");
  
  return Mat<eT>::memptr() + Mat<eT>::n_cols;
  }



//! @}
