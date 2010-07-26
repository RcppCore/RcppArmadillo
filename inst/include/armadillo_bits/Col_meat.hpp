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


//! \addtogroup Col
//! @{


//! construct an empty column vector
template<typename eT>
inline
Col<eT>::Col()
  : Mat<eT>()
  {
  arma_extra_debug_sigprint();
  }



//! construct a column vector with the specified number of n_elem
template<typename eT>
inline
Col<eT>::Col(const u32 in_n_elem)
  : Mat<eT>(in_n_elem, 1)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
Col<eT>::Col(const u32 in_n_rows, const u32 in_n_cols)
  : Mat<eT>(in_n_rows, in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from specified text
template<typename eT>
inline
Col<eT>::Col(const char* text)
  : Mat<eT>(text)
  {
  arma_extra_debug_sigprint();
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from specified text
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(text);
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from specified text
template<typename eT>
inline
Col<eT>::Col(const std::string& text)
  : Mat<eT>(text)
  {
  arma_extra_debug_sigprint();
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from specified text
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(text);
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from a given column vector
template<typename eT>
inline
Col<eT>::Col(const Col<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  }



//! construct a column vector from a given column vector
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const Col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  return *this;
  }



//! construct a column vector from a given matrix; the matrix must have exactly one column
template<typename eT>
inline
Col<eT>::Col(const Mat<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from a given matrix; the matrix must have exactly one column
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator*=(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from a given auxiliary array of eTs
template<typename eT>
inline
Col<eT>::Col(eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const bool copy_aux_mem)
  : Mat<eT>(aux_mem, aux_n_rows, aux_n_cols, copy_aux_mem)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from a given auxiliary array of eTs
template<typename eT>
inline
Col<eT>::Col(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols)
  : Mat<eT>(aux_mem, aux_n_rows, aux_n_cols)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from a given auxiliary array of eTs
template<typename eT>
inline
Col<eT>::Col(eT* aux_mem, const u32 aux_length, const bool copy_aux_mem)
  : Mat<eT>(aux_mem, aux_length, 1, copy_aux_mem)
  {
  arma_extra_debug_sigprint();
  
//   set_size(aux_length, 1);
// 
//   arma_check( (Mat<eT>::n_elem != aux_length), "Col::Col(): don't know how to handle the given array" );
// 
//   syslib::copy_elem( Mat<eT>::memptr(), aux_mem, Mat<eT>::n_elem );
  }



//! construct a column vector from a given auxiliary array of eTs
template<typename eT>
inline
Col<eT>::Col(const eT* aux_mem, const u32 aux_length)
  : Mat<eT>(aux_mem, aux_length, 1)
  {
  arma_extra_debug_sigprint();
  
//   set_size(aux_length, 1);
// 
//   arma_check( (Mat<eT>::n_elem != aux_length), "Col::Col(): don't know how to handle the given array" );
// 
//   syslib::copy_elem( Mat<eT>::memptr(), aux_mem, Mat<eT>::n_elem );
  }



template<typename eT>
template<typename T1, typename T2>
inline
Col<eT>::Col
  (
  const Base<typename Col<eT>::pod_type, T1>& A,
  const Base<typename Col<eT>::pod_type, T2>& B
  )
  : Mat<eT>(A,B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from given a submatrix; the submatrix must have exactly one column
template<typename eT>
inline
Col<eT>::Col(const subview<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from given a submatrix; the submatrix must have exactly one column
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator*=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from given a subcube; the subcube must have exactly one column
template<typename eT>
inline
Col<eT>::Col(const subview_cube<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from given a subcube; the subcube must have exactly one column
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator*=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from given a diagview
template<typename eT>
inline
Col<eT>::Col(const diagview<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from given a diagview
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator*=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
arma_inline
eT&
Col<eT>::row(const u32 row_num)
  {
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "Col::row(): out of bounds" );
  
  return access::rw(Mat<eT>::mem[row_num]);
  }



template<typename eT>
arma_inline
eT
Col<eT>::row(const u32 row_num)
  const
  {
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "Col::row(): out of bounds" );
  
  return Mat<eT>::mem[row_num];
  }



template<typename eT>
arma_inline
subview_col<eT>
Col<eT>::rows(const u32 in_row1, const u32 in_row2)
  {
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= Mat<eT>::n_rows) ), "Col::rows(): indices out of bounds or incorrectly used");
  
  return subview_col<eT>(*this, 0, in_row1, in_row2);
  }



template<typename eT>
arma_inline
const subview_col<eT>
Col<eT>::rows(const u32 in_row1, const u32 in_row2)
  const
  {
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= Mat<eT>::n_rows) ), "Col::rows(): indices out of bounds or incorrectly used");
  
  return subview_col<eT>(*this, 0, in_row1, in_row2);
  }



//! construct a column vector from Op, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename op_type>
inline
Col<eT>::Col(const Op<T1, op_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from Op, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename op_type>
inline
const Col<eT>&
Col<eT>::operator=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Col<eT>&
Col<eT>::operator*=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
Col<eT>::Col(const eOp<T1, eop_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Col<eT>&
Col<eT>::operator=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Col<eT>&
Col<eT>::operator*=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
Col<eT>::Col(const mtOp<eT, T1, op_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Col<eT>&
Col<eT>::operator=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Col<eT>&
Col<eT>::operator*=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from Glue, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Col<eT>::Col(const Glue<T1, T2, glue_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from Glue, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Col<eT>&
Col<eT>::operator=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Col<eT>&
Col<eT>::operator*=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
Col<eT>::Col(const eGlue<T1, T2, eglue_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Col<eT>&
Col<eT>::operator=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Col<eT>&
Col<eT>::operator*=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Col<eT>::Col(const mtGlue<eT, T1, T2, glue_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Col<eT>&
Col<eT>::operator=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Col<eT>&
Col<eT>::operator*=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! change the number of rows
template<typename eT>
inline
void
Col<eT>::set_size(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::set_size(in_n_elem,1);
  }



//! change the number of n_rows  (this function re-implements mat::set_size() in order to check the number of columns)
template<typename eT>
inline
void
Col<eT>::set_size(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();

  // min() is used in case in_n_cols is zero
  Mat<eT>::set_size( in_n_rows, (std::min)( u32(1), in_n_cols ) );
  
  arma_debug_check( (in_n_cols > 1), "Col::set_size(): incompatible dimensions" );
  }



//! change the number of n_rows  (this function re-implements mat::copy_size() in order to check the number of columns)
template<typename eT>
template<typename eT2>
inline
void
Col<eT>::copy_size(const Mat<eT2>& x)
  {
  arma_extra_debug_sigprint();
  
  // min() is used in case x.n_cols is zero
  Mat<eT>::set_size( x.n_rows, (std::min)( u32(1), x.n_cols ) );
  
  arma_debug_check( (x.n_cols > 1), "Col::copy_size(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Col<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros();
  }



template<typename eT>
inline
void
Col<eT>::zeros(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros(in_n_elem, 1);
  }



template<typename eT>
inline
void
Col<eT>::zeros(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  // min() is used in case in_n_cols is zero
  Mat<eT>::zeros( in_n_rows, (std::min)( u32(1), in_n_cols ) );
  
  arma_debug_check( (in_n_cols > 1), "Col::zeros(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Col<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::ones();
  }



template<typename eT>
inline
void
Col<eT>::ones(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::ones(in_n_elem, 1);
  }



template<typename eT>
inline
void
Col<eT>::ones(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  // min() is used in case in_n_cols is zero
  Mat<eT>::ones( in_n_rows, (std::min)( u32(1), in_n_cols ) );
  
  arma_debug_check( (in_n_cols > 1), "Col::ones(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Col<eT>::load(const std::string name, const file_type type, const bool print_status)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::load(name, type, print_status);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Col<eT>::load(std::istream& is, const file_type type, const bool print_status)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::load(is, type, print_status);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Col<eT>::quiet_load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::quiet_load(name, type);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Col<eT>::quiet_load(std::istream& is, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::quiet_load(is, type);
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



template<typename eT>
inline
typename Col<eT>::row_iterator
Col<eT>::begin_row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return Mat<eT>::memptr() + row_num;
  }



template<typename eT>
inline
typename Col<eT>::const_row_iterator
Col<eT>::begin_row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return Mat<eT>::memptr() + row_num;
  }



template<typename eT>
inline
typename Col<eT>::row_iterator
Col<eT>::end_row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "end_row(): index out of bounds");
  
  return Mat<eT>::memptr() + row_num + 1;
  }



template<typename eT>
inline
typename Col<eT>::const_row_iterator
Col<eT>::end_row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "end_row(): index out of bounds");
  
  return Mat<eT>::memptr() + row_num + 1;
  }



//! @}
