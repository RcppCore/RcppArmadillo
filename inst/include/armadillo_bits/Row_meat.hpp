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
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  }



template<typename eT>
inline
Row<eT>::Row(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::init(1, in_n_elem);
  }



template<typename eT>
inline
Row<eT>::Row(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::init(in_n_rows, in_n_cols);
  }



template<typename eT>
inline
Row<eT>::Row(const char* text)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(text);
  }
  


template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(text);
  
  return *this;
  }



template<typename eT>
inline
Row<eT>::Row(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(text);
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(text);
  
  return *this;
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(val);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
Row<eT>::Row(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(X.get_ref());
  }



template<typename eT>
template<typename T1>
inline
const Row<eT>&
Row<eT>::operator=(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X.get_ref());
  
  return *this;
  }



//! construct a row vector from a given auxiliary array
template<typename eT>
inline
Row<eT>::Row(eT* aux_mem, const u32 aux_length, const bool copy_aux_mem, const bool strict)
  : Mat<eT>(aux_mem, 1, aux_length, copy_aux_mem, strict)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  }



//! construct a row vector from a given auxiliary array
template<typename eT>
inline
Row<eT>::Row(const eT* aux_mem, const u32 aux_length)
  : Mat<eT>(aux_mem, 1, aux_length)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  }



template<typename eT>
template<typename T1, typename T2>
inline
Row<eT>::Row
  (
  const Base<typename Row<eT>::pod_type, T1>& A,
  const Base<typename Row<eT>::pod_type, T2>& B
  )
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::init(A,B);
  }



template<typename eT>
inline
Row<eT>::Row(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(X);
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  return *this;
  }



template<typename eT>
inline
mat_injector< Row<eT> >
Row<eT>::operator<<(const eT val)
  {
  return mat_injector< Row<eT> >(*this, val);
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



template<typename eT>
template<u32 fixed_n_elem>
arma_inline
void
Row<eT>::fixed<fixed_n_elem>::mem_setup()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(fixed_n_elem > 0)
    {
    access::rw(Mat<eT>::n_rows)    = 1;
    access::rw(Mat<eT>::n_cols)    = fixed_n_elem;
    access::rw(Mat<eT>::n_elem)    = fixed_n_elem;
    access::rw(Mat<eT>::vec_state) = 2;
    access::rw(Mat<eT>::mem_state) = 3;
    access::rw(Mat<eT>::mem)       = (fixed_n_elem > Mat_prealloc::mem_n_elem) ? mem_local_extra : Mat<eT>::mem_local;
    }
  else
    {
    access::rw(Mat<eT>::n_rows)    = 0;
    access::rw(Mat<eT>::n_cols)    = 0;
    access::rw(Mat<eT>::n_elem)    = 0;
    access::rw(Mat<eT>::vec_state) = 2;
    access::rw(Mat<eT>::mem_state) = 3;
    access::rw(Mat<eT>::mem)       = 0;
    }
  }



#ifdef ARMA_EXTRA_ROW_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_ROW_MEAT)
#endif



//! @}
