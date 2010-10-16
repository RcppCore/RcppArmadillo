// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
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



//! remove specified row
template<typename eT>
inline
void
Col<eT>::shed_row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= Mat<eT>::n_rows, "Col::shed_row(): out of bounds");
  
  shed_rows(row_num, row_num);
  }



//! remove specified rows
template<typename eT>
inline
void
Col<eT>::shed_rows(const u32 in_row1, const u32 in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= Mat<eT>::n_rows),
    "Col::shed_rows(): indices out of bounds or incorrectly used"
    );
  
  const u32 n_keep_front = in_row1;
  const u32 n_keep_back  = Mat<eT>::n_rows - (in_row2 + 1);
  
  Col<eT> X(n_keep_front + n_keep_back);
  
        eT* X_mem = X.memptr();
  const eT* t_mem = (*this).memptr();
  
  if(n_keep_front > 0)
    {
    syslib::copy_elem( X_mem, t_mem, n_keep_front );
    }
  
  if(n_keep_back > 0)
    {
    syslib::copy_elem( &(X_mem[n_keep_front]), &(t_mem[in_row2+1]), n_keep_back);
    }
  
  steal_mem(X);
  }



//! insert N rows at the specified row position,
//! optionally setting the elements of the inserted rows to zero
template<typename eT>
inline
void
Col<eT>::insert_rows(const u32 row_num, const u32 N, const bool set_to_zero)
  {
  arma_extra_debug_sigprint();
  
  const u32 t_n_rows = Mat<eT>::n_rows;
  
  const u32 A_n_rows = row_num;
  const u32 B_n_rows = t_n_rows - row_num;
  
  // insertion at row_num == n_rows is in effect an append operation
  arma_debug_check( (row_num > t_n_rows), "Col::insert_rows(): out of bounds");
  
  if(N > 0)
    {
    Col<eT> out(t_n_rows + N);
    
          eT* out_mem = out.memptr();
    const eT*   t_mem = (*this).memptr();
    
    if(A_n_rows > 0)
      {
      syslib::copy_elem( out_mem, t_mem, A_n_rows );
      }
    
    if(B_n_rows > 0)
      {
      syslib::copy_elem( &(out_mem[row_num + N]), &(t_mem[row_num]), B_n_rows );
      }
    
    if(set_to_zero == true)
      {
      arrayops::inplace_set( &(out_mem[row_num]), eT(0), N );
      }
    
    steal_mem(out);
    }
  }



//! insert the given object at the specified row position; 
//! the given object must have one column
template<typename eT>
template<typename T1>
inline
void
Col<eT>::insert_rows(const u32 row_num, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::insert_rows(row_num, X);
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
