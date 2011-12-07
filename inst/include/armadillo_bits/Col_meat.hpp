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


//! \addtogroup Col
//! @{


//! construct an empty column vector
template<typename eT>
inline
Col<eT>::Col()
  : Mat<eT>(0, 1)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



template<typename eT>
inline
Col<eT>::Col(const Col<eT>& X)
  : Mat<eT>(X.n_elem, 1)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  arrayops::copy((*this).memptr(), X.memptr(), X.n_elem);
  }



//! construct a column vector with the specified number of n_elem
template<typename eT>
inline
Col<eT>::Col(const uword in_n_elem)
  : Mat<eT>(in_n_elem, 1)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



template<typename eT>
inline
Col<eT>::Col(const uword in_n_rows, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  Mat<eT>::init_warm(in_n_rows, in_n_cols);
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



#if defined(ARMA_USE_CXX11)

template<typename eT>
inline
Col<eT>::Col(const std::initializer_list<eT>& list)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(list);
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const std::initializer_list<eT>& list)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(list);
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  return *this;
  }

#endif



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
Col<eT>::Col(eT* aux_mem, const uword aux_length, const bool copy_aux_mem, const bool strict)
  : Mat<eT>(aux_mem, aux_length, 1, copy_aux_mem, strict)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  }



//! construct a column vector from a given auxiliary array of eTs
template<typename eT>
inline
Col<eT>::Col(const eT* aux_mem, const uword aux_length)
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



template<typename eT>
template<typename T1>
inline
Col<eT>::Col(const BaseCube<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  Mat<eT>::operator=(X);
  }



template<typename eT>
template<typename T1>
inline
const Col<eT>&
Col<eT>::operator=(const BaseCube<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  return *this;
  }



template<typename eT>
inline
Col<eT>::Col(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 1;
  
  Mat<eT>::operator=(X);
  }



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
Col<eT>::row(const uword row_num)
  {
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "Col::row(): out of bounds" );
  
  return access::rw(Mat<eT>::mem[row_num]);
  }



template<typename eT>
arma_inline
eT
Col<eT>::row(const uword row_num) const
  {
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "Col::row(): out of bounds" );
  
  return Mat<eT>::mem[row_num];
  }



template<typename eT>
arma_inline
subview_col<eT>
Col<eT>::rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= Mat<eT>::n_rows) ), "Col::rows(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return subview_col<eT>(*this, 0, in_row1, subview_n_rows);
  }



template<typename eT>
arma_inline
const subview_col<eT>
Col<eT>::rows(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= Mat<eT>::n_rows) ), "Col::rows(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return subview_col<eT>(*this, 0, in_row1, subview_n_rows);
  }



template<typename eT>
arma_inline
subview_col<eT>
Col<eT>::subvec(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= Mat<eT>::n_rows) ), "Col::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return subview_col<eT>(*this, 0, in_row1, subview_n_rows);
  }



template<typename eT>
arma_inline
const subview_col<eT>
Col<eT>::subvec(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= Mat<eT>::n_rows) ), "Col::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return subview_col<eT>(*this, 0, in_row1, subview_n_rows);
  }



template<typename eT>
arma_inline
subview_col<eT>
Col<eT>::subvec(const span& row_span)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;

  const uword local_n_rows = Mat<eT>::n_rows;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword subvec_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;

  arma_debug_check( ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) ), "Col::subvec(): indices out of bounds or incorrectly used");
  
  return subview_col<eT>(*this, 0, in_row1, subvec_n_rows);
  }



template<typename eT>
arma_inline
const subview_col<eT>
Col<eT>::subvec(const span& row_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;

  const uword local_n_rows = Mat<eT>::n_rows;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword subvec_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;

  arma_debug_check( ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) ), "Col::subvec(): indices out of bounds or incorrectly used");
  
  return subview_col<eT>(*this, 0, in_row1, subvec_n_rows);
  }



// template<typename eT>
// arma_inline
// subview_col<eT>
// Col<eT>::operator()(const span& row_span)
//   {
//   arma_extra_debug_sigprint();
//   
//   return subvec(row_span);
//   }
// 
// 
// 
// template<typename eT>
// arma_inline
// const subview_col<eT>
// Col<eT>::operator()(const span& row_span) const
//   {
//   arma_extra_debug_sigprint();
//   
//   return subvec(row_span);
//   }



//! remove specified row
template<typename eT>
inline
void
Col<eT>::shed_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= Mat<eT>::n_rows, "Col::shed_row(): out of bounds");
  
  shed_rows(row_num, row_num);
  }



//! remove specified rows
template<typename eT>
inline
void
Col<eT>::shed_rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= Mat<eT>::n_rows),
    "Col::shed_rows(): indices out of bounds or incorrectly used"
    );
  
  const uword n_keep_front = in_row1;
  const uword n_keep_back  = Mat<eT>::n_rows - (in_row2 + 1);
  
  Col<eT> X(n_keep_front + n_keep_back);
  
        eT* X_mem = X.memptr();
  const eT* t_mem = (*this).memptr();
  
  if(n_keep_front > 0)
    {
    arrayops::copy( X_mem, t_mem, n_keep_front );
    }
  
  if(n_keep_back > 0)
    {
    arrayops::copy( &(X_mem[n_keep_front]), &(t_mem[in_row2+1]), n_keep_back);
    }
  
  Mat<eT>::steal_mem(X);
  }



//! insert N rows at the specified row position,
//! optionally setting the elements of the inserted rows to zero
template<typename eT>
inline
void
Col<eT>::insert_rows(const uword row_num, const uword N, const bool set_to_zero)
  {
  arma_extra_debug_sigprint();
  
  const uword t_n_rows = Mat<eT>::n_rows;
  
  const uword A_n_rows = row_num;
  const uword B_n_rows = t_n_rows - row_num;
  
  // insertion at row_num == n_rows is in effect an append operation
  arma_debug_check( (row_num > t_n_rows), "Col::insert_rows(): out of bounds");
  
  if(N > 0)
    {
    Col<eT> out(t_n_rows + N);
    
          eT* out_mem = out.memptr();
    const eT*   t_mem = (*this).memptr();
    
    if(A_n_rows > 0)
      {
      arrayops::copy( out_mem, t_mem, A_n_rows );
      }
    
    if(B_n_rows > 0)
      {
      arrayops::copy( &(out_mem[row_num + N]), &(t_mem[row_num]), B_n_rows );
      }
    
    if(set_to_zero == true)
      {
      arrayops::inplace_set( &(out_mem[row_num]), eT(0), N );
      }
    
    Mat<eT>::steal_mem(out);
    }
  }



//! insert the given object at the specified row position; 
//! the given object must have one column
template<typename eT>
template<typename T1>
inline
void
Col<eT>::insert_rows(const uword row_num, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::insert_rows(row_num, X);
  }



template<typename eT>
inline
typename Col<eT>::row_iterator
Col<eT>::begin_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return Mat<eT>::memptr() + row_num;
  }



template<typename eT>
inline
typename Col<eT>::const_row_iterator
Col<eT>::begin_row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return Mat<eT>::memptr() + row_num;
  }



template<typename eT>
inline
typename Col<eT>::row_iterator
Col<eT>::end_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "end_row(): index out of bounds");
  
  return Mat<eT>::memptr() + row_num + 1;
  }



template<typename eT>
inline
typename Col<eT>::const_row_iterator
Col<eT>::end_row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "end_row(): index out of bounds");
  
  return Mat<eT>::memptr() + row_num + 1;
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
void
Col<eT>::fixed<fixed_n_elem>::mem_setup()
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::n_rows)    = fixed_n_elem;
  access::rw(Mat<eT>::n_cols)    = 1;
  access::rw(Mat<eT>::n_elem)    = fixed_n_elem;
  access::rw(Mat<eT>::vec_state) = 1;
  access::rw(Mat<eT>::mem_state) = 3;
  access::rw(Mat<eT>::mem)       = (use_extra) ? mem_local_extra : Mat<eT>::mem_local;
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
void
Col<eT>::fixed<fixed_n_elem>::change_to_row()
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::n_cols) = fixed_n_elem;
  access::rw(Mat<eT>::n_rows) = 1;
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
Col<eT>::fixed<fixed_n_elem>::fixed()
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
Col<eT>::fixed<fixed_n_elem>::fixed(const fixed<fixed_n_elem>& X)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  eT* dest = (use_extra) ? mem_local_extra : Mat<eT>::mem_local;
  
  arrayops::copy( dest, X.mem, fixed_n_elem );
  }



template<typename eT>
template<uword fixed_n_elem>
inline
Col<eT>::fixed<fixed_n_elem>::fixed(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  Col<eT>::operator=(X);
  }



template<typename eT>
template<uword fixed_n_elem>
template<typename T1>
inline
Col<eT>::fixed<fixed_n_elem>::fixed(const Base<eT,T1>& A)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  Col<eT>::operator=(A.get_ref());
  }



template<typename eT>
template<uword fixed_n_elem>
template<typename T1, typename T2>
inline
Col<eT>::fixed<fixed_n_elem>::fixed(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  Col<eT>::init(A,B);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
Col<eT>::fixed<fixed_n_elem>::fixed(eT* aux_mem, const bool copy_aux_mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  access::rw(Mat<eT>::n_rows)    = fixed_n_elem;
  access::rw(Mat<eT>::n_cols)    = 1;
  access::rw(Mat<eT>::n_elem)    = fixed_n_elem;
  access::rw(Mat<eT>::vec_state) = 1;
  access::rw(Mat<eT>::mem_state) = 3;
  
  if(copy_aux_mem == true)
    {
    eT* dest = (use_extra) ? mem_local_extra : Mat<eT>::mem_local;
    
    access::rw(Mat<eT>::mem) = dest;
    
    arrayops::copy( dest, aux_mem, fixed_n_elem );
    }
  else
    {
    access::rw(Mat<eT>::mem) = aux_mem;
    }
  }



template<typename eT>
template<uword fixed_n_elem>
inline
Col<eT>::fixed<fixed_n_elem>::fixed(const eT* aux_mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  arrayops::copy( const_cast<eT*>(Mat<eT>::mem), aux_mem, fixed_n_elem );
  }



//! NOTE: this function relies on
//! Col::operator=(text), to change vec_state as well as swapping n_rows and n_cols,
//! and Mat::init(), to check that the given vector will not have a different size than fixed_n_elem.
template<typename eT>
template<uword fixed_n_elem>
inline
Col<eT>::fixed<fixed_n_elem>::fixed(const char* text)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  change_to_row();
  
  Col<eT>::operator=(text);
  }



//! NOTE: this function relies on
//! Col::operator=(text), to change vec_state as well as swapping n_rows and n_cols,
//! and Mat::init(), to check that the given vector will not have a different size than fixed_n_elem.
template<typename eT>
template<uword fixed_n_elem>
inline
Col<eT>::fixed<fixed_n_elem>::fixed(const std::string& text)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  change_to_row();
  
  Col<eT>::operator=(text);
  }



template<typename eT>
template<uword fixed_n_elem>
template<typename T1>
const Col<eT>&
Col<eT>::fixed<fixed_n_elem>::operator=(const Base<eT,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  Col<eT>::operator=(A.get_ref());
  
  return *this;
  }



template<typename eT>
template<uword fixed_n_elem>
const Col<eT>&
Col<eT>::fixed<fixed_n_elem>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Col<eT>::operator=(val);
  
  return *this;
  }



template<typename eT>
template<uword fixed_n_elem>
const Col<eT>&
Col<eT>::fixed<fixed_n_elem>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  change_to_row();
  
  Col<eT>::operator=(text);
  
  return *this; 
  }



template<typename eT>
template<uword fixed_n_elem>
const Col<eT>&
Col<eT>::fixed<fixed_n_elem>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  change_to_row();
  
  Col<eT>::operator=(text);
  
  return *this; 
  }



template<typename eT>
template<uword fixed_n_elem>
const Col<eT>&
Col<eT>::fixed<fixed_n_elem>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Col<eT>::operator=(X);
  
  return *this; 
  }



template<typename eT>
template<uword fixed_n_elem>
inline
subview_row<eT>
Col<eT>::fixed<fixed_n_elem>::operator()(const uword row_num, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_num, col_span);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
const subview_row<eT>
Col<eT>::fixed<fixed_n_elem>::operator()(const uword row_num, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_num, col_span);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
subview_col<eT>
Col<eT>::fixed<fixed_n_elem>::operator()(const span& row_span, const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_span, col_num);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
const subview_col<eT>
Col<eT>::fixed<fixed_n_elem>::operator()(const span& row_span, const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_span, col_num);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
subview<eT>
Col<eT>::fixed<fixed_n_elem>::operator()(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_span, col_span);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
const subview<eT>
Col<eT>::fixed<fixed_n_elem>::operator()(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_span, col_span);
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Col<eT>::fixed<fixed_n_elem>::operator[] (const uword i)
  {
  return access::rw( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Col<eT>::fixed<fixed_n_elem>::operator[] (const uword i) const
  {
  return ( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Col<eT>::fixed<fixed_n_elem>::at(const uword i)
  {
  return access::rw( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Col<eT>::fixed<fixed_n_elem>::at(const uword i) const
  {
  return ( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Col<eT>::fixed<fixed_n_elem>::operator() (const uword i)
  {
  arma_debug_check( (i >= fixed_n_elem), "Col::fixed::operator(): out of bounds");
  
  return access::rw( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Col<eT>::fixed<fixed_n_elem>::operator() (const uword i) const
  {
  arma_debug_check( (i >= fixed_n_elem), "Col::fixed::operator(): out of bounds");
  
  return ( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Col<eT>::fixed<fixed_n_elem>::at(const uword in_row, const uword in_col)
  {
  return access::rw( Mat<eT>::mem[in_row] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Col<eT>::fixed<fixed_n_elem>::at(const uword in_row, const uword in_col) const
  {
  return ( Mat<eT>::mem[in_row] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Col<eT>::fixed<fixed_n_elem>::operator() (const uword in_row, const uword in_col)
  {
  arma_debug_check( ((in_row >= fixed_n_elem) || (in_col >= 1)), "Col::fixed::operator(): out of bounds" );
  
  return access::rw( Mat<eT>::mem[in_row] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Col<eT>::fixed<fixed_n_elem>::operator() (const uword in_row, const uword in_col) const
  {
  arma_debug_check( ((in_row >= fixed_n_elem) || (in_col >= 1)), "Col::fixed::operator(): out of bounds" );
  
  return ( Mat<eT>::mem[in_row] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_hot
inline
const Col<eT>&
Col<eT>::fixed<fixed_n_elem>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set( const_cast<eT*>(Mat<eT>::mem), val, fixed_n_elem );
  
  return *this;
  }



template<typename eT>
template<uword fixed_n_elem>
arma_hot
inline
const Col<eT>&
Col<eT>::fixed<fixed_n_elem>::zeros()
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set( const_cast<eT*>(Mat<eT>::mem), eT(0), fixed_n_elem );
  
  return *this;
  }



template<typename eT>
template<uword fixed_n_elem>
arma_hot
inline
const Col<eT>&
Col<eT>::fixed<fixed_n_elem>::ones()
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set( const_cast<eT*>(Mat<eT>::mem), eT(1), fixed_n_elem );
  
  return *this;
  }



#ifdef ARMA_EXTRA_COL_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_COL_MEAT)
#endif



//! @}
