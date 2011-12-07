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


//! \addtogroup Row
//! @{


//! construct an empty row vector
template<typename eT>
inline
Row<eT>::Row()
  : Mat<eT>(1, 0)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  }



template<typename eT>
inline
Row<eT>::Row(const Row<eT>& X)
  : Mat<eT>(1, X.n_elem)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  arrayops::copy((*this).memptr(), X.memptr(), X.n_elem);
  }



//! construct a row vector with the specified number of n_elem
template<typename eT>
inline
Row<eT>::Row(const uword in_n_elem)
  : Mat<eT>(1, in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  }



template<typename eT>
inline
Row<eT>::Row(const uword in_n_rows, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::init_warm(in_n_rows, in_n_cols);
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



#if defined(ARMA_USE_CXX11)

template<typename eT>
inline
Row<eT>::Row(const std::initializer_list<eT>& list)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(list);
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const std::initializer_list<eT>& list)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(list);
  
  return *this;
  }

#endif



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
Row<eT>::Row(eT* aux_mem, const uword aux_length, const bool copy_aux_mem, const bool strict)
  : Mat<eT>(aux_mem, 1, aux_length, copy_aux_mem, strict)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  }



//! construct a row vector from a given auxiliary array
template<typename eT>
inline
Row<eT>::Row(const eT* aux_mem, const uword aux_length)
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
template<typename T1>
inline
Row<eT>::Row(const BaseCube<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::vec_state) = 2;
  
  Mat<eT>::operator=(X);
  }



template<typename eT>
template<typename T1>
inline
const Row<eT>&
Row<eT>::operator=(const BaseCube<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  return *this;
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
Row<eT>::col(const uword col_num)
  {
  arma_debug_check( (col_num >= Mat<eT>::n_cols), "Row::col(): out of bounds" );
  
  return access::rw(Mat<eT>::mem[col_num]);
  }



template<typename eT>
arma_inline
eT
Row<eT>::col(const uword col_num) const
  {
  arma_debug_check( (col_num >= Mat<eT>::n_cols), "Row::col(): out of bounds" );
  
  return Mat<eT>::mem[col_num];
  }



template<typename eT>
arma_inline
subview_row<eT>
Row<eT>::cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "Row::cols(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return subview_row<eT>(*this, 0, in_col1, subview_n_cols);
  }



template<typename eT>
arma_inline
const subview_row<eT>
Row<eT>::cols(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "Row::cols(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return subview_row<eT>(*this, 0, in_col1, subview_n_cols);
  }



template<typename eT>
arma_inline
subview_row<eT>
Row<eT>::subvec(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "Row::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return subview_row<eT>(*this, 0, in_col1, subview_n_cols);
  }



template<typename eT>
arma_inline
const subview_row<eT>
Row<eT>::subvec(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "Row::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return subview_row<eT>(*this, 0, in_col1, subview_n_cols);
  }



template<typename eT>
arma_inline
subview_row<eT>
Row<eT>::subvec(const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const bool col_all = col_span.whole;

  const uword local_n_cols = Mat<eT>::n_cols;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword subvec_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;

  arma_debug_check( ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) ), "Row::subvec(): indices out of bounds or incorrectly used");
  
  return subview_row<eT>(*this, 0, in_col1, subvec_n_cols);
  }



template<typename eT>
arma_inline
const subview_row<eT>
Row<eT>::subvec(const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool col_all = col_span.whole;

  const uword local_n_cols = Mat<eT>::n_cols;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword subvec_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;

  arma_debug_check( ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) ), "Row::subvec(): indices out of bounds or incorrectly used");
  
  return subview_row<eT>(*this, 0, in_col1, subvec_n_cols);
  }



// template<typename eT>
// arma_inline
// subview_row<eT>
// Row<eT>::operator()(const span& col_span)
//   {
//   arma_extra_debug_sigprint();
//   
//   return subvec(col_span);
//   }
// 
// 
// 
// template<typename eT>
// arma_inline
// const subview_row<eT>
// Row<eT>::operator()(const span& col_span) const
//   {
//   arma_extra_debug_sigprint();
//   
//   return subvec(col_span);
//   }



//! remove specified columns
template<typename eT>
inline
void
Row<eT>::shed_col(const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= Mat<eT>::n_cols, "Row::shed_col(): out of bounds");
  
  shed_cols(col_num, col_num);
  }



//! remove specified columns
template<typename eT>
inline
void
Row<eT>::shed_cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols),
    "Row::shed_cols(): indices out of bounds or incorrectly used"
    );
  
  const uword n_keep_front = in_col1;
  const uword n_keep_back  = Mat<eT>::n_cols - (in_col2 + 1);
  
  Row<eT> X(n_keep_front + n_keep_back);
  
        eT* X_mem = X.memptr();
  const eT* t_mem = (*this).memptr();
  
  if(n_keep_front > 0)
    {
    arrayops::copy( X_mem, t_mem, n_keep_front );
    }
  
  if(n_keep_back > 0)
    {
    arrayops::copy( &(X_mem[n_keep_front]), &(t_mem[in_col2+1]), n_keep_back);
    }
  
  Mat<eT>::steal_mem(X);
  }



//! insert N cols at the specified col position,
//! optionally setting the elements of the inserted cols to zero
template<typename eT>
inline
void
Row<eT>::insert_cols(const uword col_num, const uword N, const bool set_to_zero)
  {
  arma_extra_debug_sigprint();
  
  const uword t_n_cols = Mat<eT>::n_cols;
  
  const uword A_n_cols = col_num;
  const uword B_n_cols = t_n_cols - col_num;
  
  // insertion at col_num == n_cols is in effect an append operation
  arma_debug_check( (col_num > t_n_cols), "Row::insert_cols(): out of bounds");
  
  if(N > 0)
    {
    Row<eT> out(t_n_cols + N);
    
          eT* out_mem = out.memptr();
    const eT*   t_mem = (*this).memptr();
    
    if(A_n_cols > 0)
      {
      arrayops::copy( out_mem, t_mem, A_n_cols );
      }
    
    if(B_n_cols > 0)
      {
      arrayops::copy( &(out_mem[col_num + N]), &(t_mem[col_num]), B_n_cols );
      }
    
    if(set_to_zero == true)
      {
      arrayops::inplace_set( &(out_mem[col_num]), eT(0), N );
      }
    
    Mat<eT>::steal_mem(out);
    }
  }



//! insert the given object at the specified col position; 
//! the given object must have one row
template<typename eT>
template<typename T1>
inline
void
Row<eT>::insert_cols(const uword col_num, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::insert_cols(col_num, X);
  }



template<typename eT>
inline
typename Row<eT>::row_iterator
Row<eT>::begin_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return Mat<eT>::memptr();
  }



template<typename eT>
inline
typename Row<eT>::const_row_iterator
Row<eT>::begin_row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return Mat<eT>::memptr();
  }



template<typename eT>
inline
typename Row<eT>::row_iterator
Row<eT>::end_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "end_row(): index out of bounds");
  
  return Mat<eT>::memptr() + Mat<eT>::n_cols;
  }



template<typename eT>
inline
typename Row<eT>::const_row_iterator
Row<eT>::end_row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= Mat<eT>::n_rows), "end_row(): index out of bounds");
  
  return Mat<eT>::memptr() + Mat<eT>::n_cols;
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
void
Row<eT>::fixed<fixed_n_elem>::mem_setup()
  {
  arma_extra_debug_sigprint();
  
  access::rw(Mat<eT>::n_rows)    = 1;
  access::rw(Mat<eT>::n_cols)    = fixed_n_elem;
  access::rw(Mat<eT>::n_elem)    = fixed_n_elem;
  access::rw(Mat<eT>::vec_state) = 2;
  access::rw(Mat<eT>::mem_state) = 3;
  access::rw(Mat<eT>::mem)       = (use_extra) ? mem_local_extra : Mat<eT>::mem_local;
  }



template<typename eT>
template<uword fixed_n_elem>
inline
Row<eT>::fixed<fixed_n_elem>::fixed()
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
Row<eT>::fixed<fixed_n_elem>::fixed(const fixed<fixed_n_elem>& X)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  eT* dest = (use_extra) ? mem_local_extra : Mat<eT>::mem_local;
  
  arrayops::copy( dest, X.mem, fixed_n_elem );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
Row<eT>::fixed<fixed_n_elem>::fixed(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  Row<eT>::operator=(X);
  }



template<typename eT>
template<uword fixed_n_elem>
template<typename T1>
arma_inline
Row<eT>::fixed<fixed_n_elem>::fixed(const Base<eT,T1>& A)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  Row<eT>::operator=(A.get_ref());
  }



template<typename eT>
template<uword fixed_n_elem>
template<typename T1, typename T2>
arma_inline
Row<eT>::fixed<fixed_n_elem>::fixed(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  Row<eT>::init(A,B);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
Row<eT>::fixed<fixed_n_elem>::fixed(eT* aux_mem, const bool copy_aux_mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  access::rw(Mat<eT>::n_rows)    = 1;
  access::rw(Mat<eT>::n_cols)    = fixed_n_elem;
  access::rw(Mat<eT>::n_elem)    = fixed_n_elem;
  access::rw(Mat<eT>::vec_state) = 2;
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
Row<eT>::fixed<fixed_n_elem>::fixed(const eT* aux_mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  arrayops::copy( const_cast<eT*>(Mat<eT>::mem), aux_mem, fixed_n_elem );
  }



template<typename eT>
template<uword fixed_n_elem>
inline
Row<eT>::fixed<fixed_n_elem>::fixed(const char* text)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  Row<eT>::operator=(text);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
Row<eT>::fixed<fixed_n_elem>::fixed(const std::string& text)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem_setup();
  
  Row<eT>::operator=(text);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
subview_row<eT>
Row<eT>::fixed<fixed_n_elem>::operator()(const uword row_num, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_num, col_span);
  }



template<typename eT>
template<uword fixed_n_elem>
template<typename T1>
const Row<eT>&
Row<eT>::fixed<fixed_n_elem>::operator=(const Base<eT,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  Row<eT>::operator=(A.get_ref());
  
  return *this;
  }



template<typename eT>
template<uword fixed_n_elem>
const Row<eT>&
Row<eT>::fixed<fixed_n_elem>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Row<eT>::operator=(val);
  
  return *this;
  }



template<typename eT>
template<uword fixed_n_elem>
const Row<eT>&
Row<eT>::fixed<fixed_n_elem>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  Row<eT>::operator=(text);
  
  return *this; 
  }



template<typename eT>
template<uword fixed_n_elem>
const Row<eT>&
Row<eT>::fixed<fixed_n_elem>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  Row<eT>::operator=(text);
  
  return *this; 
  }



template<typename eT>
template<uword fixed_n_elem>
const Row<eT>&
Row<eT>::fixed<fixed_n_elem>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Row<eT>::operator=(X);
  
  return *this; 
  }



template<typename eT>
template<uword fixed_n_elem>
inline
const subview_row<eT>
Row<eT>::fixed<fixed_n_elem>::operator()(const uword row_num, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_num, col_span);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
subview_col<eT>
Row<eT>::fixed<fixed_n_elem>::operator()(const span& row_span, const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_span, col_num);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
const subview_col<eT>
Row<eT>::fixed<fixed_n_elem>::operator()(const span& row_span, const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_span, col_num);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
subview<eT>
Row<eT>::fixed<fixed_n_elem>::operator()(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_span, col_span);
  }



template<typename eT>
template<uword fixed_n_elem>
inline
const subview<eT>
Row<eT>::fixed<fixed_n_elem>::operator()(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return Mat<eT>::operator()(row_span, col_span);
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Row<eT>::fixed<fixed_n_elem>::operator[] (const uword i)
  {
  return access::rw( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Row<eT>::fixed<fixed_n_elem>::operator[] (const uword i) const
  {
  return ( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Row<eT>::fixed<fixed_n_elem>::at(const uword i)
  {
  return access::rw( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Row<eT>::fixed<fixed_n_elem>::at(const uword i) const
  {
  return ( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Row<eT>::fixed<fixed_n_elem>::operator() (const uword i)
  {
  arma_debug_check( (i >= fixed_n_elem), "Row::fixed::operator(): out of bounds");
  
  return access::rw( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Row<eT>::fixed<fixed_n_elem>::operator() (const uword i) const
  {
  arma_debug_check( (i >= fixed_n_elem), "Row::fixed::operator(): out of bounds");
  
  return ( Mat<eT>::mem[i] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Row<eT>::fixed<fixed_n_elem>::at(const uword in_row, const uword in_col)
  {
  return access::rw( Mat<eT>::mem[in_col] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Row<eT>::fixed<fixed_n_elem>::at(const uword in_row, const uword in_col) const
  {
  return ( Mat<eT>::mem[in_col] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT&
Row<eT>::fixed<fixed_n_elem>::operator() (const uword in_row, const uword in_col)
  {
  arma_debug_check( ((in_row >= 1) || (in_col >= fixed_n_elem)), "Row::fixed::operator(): out of bounds" );
  
  return access::rw( Mat<eT>::mem[in_col] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_inline
arma_warn_unused
eT
Row<eT>::fixed<fixed_n_elem>::operator() (const uword in_row, const uword in_col) const
  {
  arma_debug_check( ((in_row >= 1) || (in_col >= fixed_n_elem)), "Row::fixed::operator(): out of bounds" );
  
  return ( Mat<eT>::mem[in_col] );
  }



template<typename eT>
template<uword fixed_n_elem>
arma_hot
inline
const Row<eT>&
Row<eT>::fixed<fixed_n_elem>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set( const_cast<eT*>(Mat<eT>::mem), val, fixed_n_elem );
  
  return *this;
  }



template<typename eT>
template<uword fixed_n_elem>
arma_hot
inline
const Row<eT>&
Row<eT>::fixed<fixed_n_elem>::zeros()
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set( const_cast<eT*>(Mat<eT>::mem), eT(0), fixed_n_elem );
  
  return *this;
  }



template<typename eT>
template<uword fixed_n_elem>
arma_hot
inline
const Row<eT>&
Row<eT>::fixed<fixed_n_elem>::ones()
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set( const_cast<eT*>(Mat<eT>::mem), eT(1), fixed_n_elem );
  
  return *this;
  }



#ifdef ARMA_EXTRA_ROW_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_ROW_MEAT)
#endif



//! @}
