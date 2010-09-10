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
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



//! construct a column vector with the specified number of n_elem
template<typename eT>
inline
Col<eT>::Col(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  Mat<eT>::init(in_n_elem, 1);
  }



template<typename eT>
inline
Col<eT>::Col(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  Mat<eT>::init(in_n_rows, in_n_cols);
  }



//! construct a column vector from specified text
template<typename eT>
inline
Col<eT>::Col(const char* text)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(text);
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



//! construct a column vector from specified text
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(text);
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  return *this;
  }



//! construct a column vector from specified text
template<typename eT>
inline
Col<eT>::Col(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(text);
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



//! construct a column vector from specified text
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(text);
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  return *this;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(val);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
Col<eT>::Col(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  Mat<eT>::operator=(X.get_ref());
  }



template<typename eT>
template<typename T1>
inline
const Col<eT>&
Col<eT>::operator=(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X.get_ref());
  
  return *this;
  }



//! construct a column vector from a given auxiliary array of eTs
template<typename eT>
inline
Col<eT>::Col(eT* aux_mem, const u32 aux_length, const bool copy_aux_mem, const bool strict)
  : Mat<eT>(aux_mem, aux_length, 1, copy_aux_mem, strict)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



//! construct a column vector from a given auxiliary array of eTs
template<typename eT>
inline
Col<eT>::Col(const eT* aux_mem, const u32 aux_length)
  : Mat<eT>(aux_mem, aux_length, 1)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



template<typename eT>
template<typename T1, typename T2>
inline
Col<eT>::Col
  (
  const Base<typename Col<eT>::pod_type, T1>& A,
  const Base<typename Col<eT>::pod_type, T2>& B
  )
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  Mat<eT>::init(A,B);
  }



//! construct a column vector from given a subcube; the subcube must have exactly one column
template<typename eT>
inline
Col<eT>::Col(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  Mat<eT>::operator=(X);
  }



//! construct a column vector from given a subcube; the subcube must have exactly one column
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  return *this;
  }



template<typename eT>
inline
mat_injector< Col<eT> >
Col<eT>::operator<<(const eT val)
  {
  return mat_injector< Col<eT> >(*this, val);
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



template<typename eT>
template<u32 fixed_n_elem>
arma_inline
void
Col<eT>::fixed<fixed_n_elem>::mem_setup()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(fixed_n_elem > 0)
    {
    access::rw(Mat<eT>::n_rows)    = fixed_n_elem;
    access::rw(Mat<eT>::n_cols)    = 1;
    access::rw(Mat<eT>::n_elem)    = fixed_n_elem;
    access::rw(Mat<eT>::vec_state) = 1;
    access::rw(Mat<eT>::mem_state) = 3;
    access::rw(Mat<eT>::mem)       = (fixed_n_elem > Mat_prealloc::mem_n_elem) ? mem_local_extra : Mat<eT>::mem_local;
    }
  else
    {
    access::rw(Mat<eT>::n_rows)    = 0;
    access::rw(Mat<eT>::n_cols)    = 0;
    access::rw(Mat<eT>::n_elem)    = 0;
    access::rw(Mat<eT>::vec_state) = 1;
    access::rw(Mat<eT>::mem_state) = 3;
    access::rw(Mat<eT>::mem)       = 0;
    }
  }



#ifdef ARMA_EXTRA_COL_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_COL_MEAT)
#endif



//! @}
