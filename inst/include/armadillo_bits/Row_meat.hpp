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



//! remove specified columns
template<typename eT>
inline
void
Row<eT>::shed_col(const u32 col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= Mat<eT>::n_cols, "Row::shed_col(): out of bounds");
  
  shed_cols(col_num, col_num);
  }



//! remove specified columns
template<typename eT>
inline
void
Row<eT>::shed_cols(const u32 in_col1, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols),
    "Row::shed_cols(): indices out of bounds or incorrectly used"
    );
  
  const u32 n_keep_front = in_col1;
  const u32 n_keep_back  = Mat<eT>::n_cols - (in_col2 + 1);
  
  Row<eT> X(n_keep_front + n_keep_back);
  
        eT* X_mem = X.memptr();
  const eT* t_mem = (*this).memptr();
  
  if(n_keep_front > 0)
    {
    syslib::copy_elem( X_mem, t_mem, n_keep_front );
    }
  
  if(n_keep_back > 0)
    {
    syslib::copy_elem( &(X_mem[n_keep_front]), &(t_mem[in_col2+1]), n_keep_back);
    }
  
  steal_mem(X);
  }



//! insert N cols at the specified col position,
//! optionally setting the elements of the inserted cols to zero
template<typename eT>
inline
void
Row<eT>::insert_cols(const u32 col_num, const u32 N, const bool set_to_zero)
  {
  arma_extra_debug_sigprint();
  
  const u32 t_n_cols = Mat<eT>::n_cols;
  
  const u32 A_n_cols = col_num;
  const u32 B_n_cols = t_n_cols - col_num;
  
  // insertion at col_num == n_cols is in effect an append operation
  arma_debug_check( (col_num > t_n_cols), "Row::insert_cols(): out of bounds");
  
  if(N > 0)
    {
    Row<eT> out(t_n_cols + N);
    
          eT* out_mem = out.memptr();
    const eT*   t_mem = (*this).memptr();
    
    if(A_n_cols > 0)
      {
      syslib::copy_elem( out_mem, t_mem, A_n_cols );
      }
    
    if(B_n_cols > 0)
      {
      syslib::copy_elem( &(out_mem[col_num + N]), &(t_mem[col_num]), B_n_cols );
      }
    
    if(set_to_zero == true)
      {
      arrayops::inplace_set( &(out_mem[col_num]), eT(0), N );
      }
    
    steal_mem(out);
    }
  }



//! insert the given object at the specified col position; 
//! the given object must have one row
template<typename eT>
template<typename T1>
inline
void
Row<eT>::insert_cols(const u32 col_num, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::insert_cols(col_num, X);
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
