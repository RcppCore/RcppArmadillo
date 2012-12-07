// Copyright (C) 2011-2012 Ryan Curtin
// Copyright (C) 2011 Matthew Amidon
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup SpRow
//! @{



template<typename eT>
inline
SpRow<eT>::SpRow()
  : SpMat<eT>(1, 0)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 2;
  }



template<typename eT>
inline
SpRow<eT>::SpRow(const uword in_n_elem)
  : SpMat<eT>(1, in_n_elem)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 2;
  }



template<typename eT>
inline
SpRow<eT>::SpRow(const uword in_n_rows, const uword in_n_cols)
  : SpMat<eT>(in_n_rows, in_n_cols)
  {
  arma_extra_debug_sigprint();

  arma_debug_check((in_n_rows != 1), "SpRow::SpRow(): must have only one row");

  access::rw(SpMat<eT>::vec_state) = 2;
  }



template<typename eT>
inline
SpRow<eT>::SpRow(const char* text)
  : SpMat<eT>(text)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 2;

  arma_debug_check((SpMat<eT>::n_rows != 1), "SpRow::SpRow(): must have only one row");
  }
  


template<typename eT>
inline
const SpRow<eT>&
SpRow<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 2;
  
  SpMat<eT>::operator=(text);
  
  return *this;
  }



template<typename eT>
inline
SpRow<eT>::SpRow(const std::string& text)
  : SpMat<eT>(text)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 2;
  
  arma_debug_check((SpMat<eT>::n_rows != 1), "SpRow::SpRow(): must have only one row");
  }



template<typename eT>
inline
const SpRow<eT>&
SpRow<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT>::operator=(text);
  
  return *this;
  }



template<typename eT>
inline
const SpRow<eT>&
SpRow<eT>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT>::operator=(val);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpRow<eT>::SpRow(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 2;
  
  SpMat<eT>::operator=(X.get_ref());
  }



template<typename eT>
template<typename T1>
inline
const SpRow<eT>&
SpRow<eT>::operator=(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT>::operator=(X.get_ref());
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpRow<eT>::SpRow(const SpBase<eT,T1>& X)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 2;
  
  SpMat<eT>::operator=(X.get_ref());
  }



template<typename eT>
template<typename T1>
inline
const SpRow<eT>&
SpRow<eT>::operator=(const SpBase<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT>::operator=(X.get_ref());
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2>
inline
SpRow<eT>::SpRow
  (
  const SpBase<typename SpRow<eT>::pod_type, T1>& A,
  const SpBase<typename SpRow<eT>::pod_type, T2>& B
  )
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 2;
  
  SpMat<eT>::init(A,B);
  }



template<typename eT>
inline
SpValProxy< SpMat<eT> >
SpRow<eT>::col(const uword col_num)
  {
  arma_debug_check( (col_num >= SpMat<eT>::n_cols), "SpRow::col(): out of bounds" );
  
  return SpMat<eT>::at(0, col_num);
  }



template<typename eT>
inline
eT
SpRow<eT>::col(const uword col_num) const
  {
  arma_debug_check( (col_num >= SpMat<eT>::n_cols), "SpRow::col(): out of bounds" );
  
  return SpMat<eT>::at(0, col_num);
  }


/*
template<typename eT>
arma_inline
subview_row<eT>
SpRow<eT>::cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "SpRow::cols(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return subview_row<eT>(*this, 0, in_col1, subview_n_cols);
  }



template<typename eT>
arma_inline
const subview_row<eT>
SpRow<eT>::cols(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "SpRow::cols(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return subview_row<eT>(*this, 0, in_col1, subview_n_cols);
  }



template<typename eT>
arma_inline
subview_row<eT>
SpRow<eT>::subvec(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "SpRow::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return subview_row<eT>(*this, 0, in_col1, subview_n_cols);
  }



template<typename eT>
arma_inline
const subview_row<eT>
SpRow<eT>::subvec(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= Mat<eT>::n_cols) ), "SpRow::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return subview_row<eT>(*this, 0, in_col1, subview_n_cols);
  }



template<typename eT>
arma_inline
subview_row<eT>
SpRow<eT>::subvec(const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const bool col_all = col_span.whole;

  const uword local_n_cols = Mat<eT>::n_cols;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword subvec_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;

  arma_debug_check( ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) ), "SpRow::subvec(): indices out of bounds or incorrectly used");
  
  return subview_row<eT>(*this, 0, in_col1, subvec_n_cols);
  }



template<typename eT>
arma_inline
const subview_row<eT>
SpRow<eT>::subvec(const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool col_all = col_span.whole;

  const uword local_n_cols = Mat<eT>::n_cols;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword subvec_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;

  arma_debug_check( ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) ), "SpRow::subvec(): indices out of bounds or incorrectly used");
  
  return subview_row<eT>(*this, 0, in_col1, subvec_n_cols);
  }
*/


// template<typename eT>
// arma_inline
// subview_row<eT>
// SpRow<eT>::operator()(const span& col_span)
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
// SpRow<eT>::operator()(const span& col_span) const
//   {
//   arma_extra_debug_sigprint();
//   
//   return subvec(col_span);
//   }



//! remove specified columns
template<typename eT>
inline
void
SpRow<eT>::shed_col(const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= SpMat<eT>::n_cols, "SpRow::shed_col(): out of bounds");
  
  shed_cols(col_num, col_num);
  }



//! remove specified columns
template<typename eT>
inline
void
SpRow<eT>::shed_cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= SpMat<eT>::n_cols),
    "SpRow::shed_cols(): indices out of bounds or incorrectly used"
    );
  
  const uword diff = (in_col2 - in_col1 + 1);

  // This is doubleplus easy because we have all the column pointers stored.
  const uword start = SpMat<eT>::col_ptrs[in_col1];
  const uword end   = SpMat<eT>::col_ptrs[in_col2 + 1];

  if (start != end)
    {
    const uword elem_diff = end - start;

    eT*    new_values      = memory::acquire_chunked<eT>   (SpMat<eT>::n_nonzero - elem_diff);
    uword* new_row_indices = memory::acquire_chunked<uword>(SpMat<eT>::n_nonzero - elem_diff);

    // Copy first set of elements, if necessary.
    if (start > 0)
      {
      arrayops::copy(new_values, SpMat<eT>::values, start);
      arrayops::copy(new_row_indices, SpMat<eT>::row_indices, start);
      }

    // Copy last set of elements, if necessary.
    if (end != SpMat<eT>::n_nonzero)
      {
      arrayops::copy(new_values + start, SpMat<eT>::values + end, (SpMat<eT>::n_nonzero - end));
      arrayops::copy(new_row_indices + start, SpMat<eT>::row_indices + end, (SpMat<eT>::n_nonzero - end));
      }

    memory::release(SpMat<eT>::values);
    memory::release(SpMat<eT>::row_indices);

    access::rw(SpMat<eT>::values) = new_values;
    access::rw(SpMat<eT>::row_indices) = new_row_indices;

    access::rw(SpMat<eT>::n_nonzero) -= elem_diff;
    }

  // Update column pointers.
  uword* new_col_ptrs = memory::acquire<uword>(SpMat<eT>::n_cols - diff + 1);

  // Copy first part of column pointers.
  if (in_col1 > 0)
    {
    arrayops::copy(new_col_ptrs, SpMat<eT>::col_ptrs, in_col1);
    }

  // Copy last part of column pointers (and adjust their values as necessary).
  if (in_col2 < SpMat<eT>::n_cols - 1)
    {
    arrayops::copy(new_col_ptrs + in_col1, SpMat<eT>::col_ptrs + in_col2 + 1, SpMat<eT>::n_cols - in_col2);
    // Modify their values.
    arrayops::inplace_minus(new_col_ptrs + in_col1, (end - start), SpMat<eT>::n_cols - in_col2);
    }

  memory::release(SpMat<eT>::col_ptrs);

  access::rw(SpMat<eT>::col_ptrs) = new_col_ptrs;

  access::rw(SpMat<eT>::n_cols) -= diff;
  access::rw(SpMat<eT>::n_elem) -= diff;
  }



// //! insert N cols at the specified col position,
// //! optionally setting the elements of the inserted cols to zero
// template<typename eT>
// inline
// void
// SpRow<eT>::insert_cols(const uword col_num, const uword N, const bool set_to_zero)
//   {
//   arma_extra_debug_sigprint();
// 
//   // insertion at col_num == n_cols is in effect an append operation
//   arma_debug_check( (col_num > SpMat<eT>::n_cols), "SpRow::insert_cols(): out of bounds");
// 
//   arma_debug_check( (set_to_zero == false), "SpRow::insert_cols(): cannot set elements to nonzero values");
// 
//   uword newVal = (col_num == 0) ? 0 : SpMat<eT>::col_ptrs[col_num];
//   SpMat<eT>::col_ptrs.insert(col_num, N, newVal);
//   uword* new_col_ptrs = memory::acquire<uword>(SpMat<eT>::n_cols + N);
// 
//   arrayops::copy(new_col_ptrs, SpMat<eT>::col_ptrs, col_num);
// 
//   uword fill_value = (col_num == 0) ? 0 : SpMat<eT>::col_ptrs[col_num - 1];
//   arrayops::inplace_set(new_col_ptrs + col_num, fill_value, N);
// 
//   arrayops::copy(new_col_ptrs + col_num + N, SpMat<eT>::col_ptrs + col_num, SpMat<eT>::n_cols - col_num);
// 
//   access::rw(SpMat<eT>::n_cols) += N;
//   access::rw(SpMat<eT>::n_elem) += N;
//   }
// 
// 
// 
// //! insert the given object at the specified col position; 
// //! the given object must have one row
// template<typename eT>
// template<typename T1>
// inline
// void
// SpRow<eT>::insert_cols(const uword col_num, const Base<eT,T1>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   SpMat<eT>::insert_cols(col_num, X);
//   }



template<typename eT>
inline
typename SpRow<eT>::row_iterator
SpRow<eT>::begin_row(const uword row_num)
  {
  arma_extra_debug_sigprint();

  // Since this is a row, row_num can only be 0.  But the option is provided for
  // compatibility.
  arma_debug_check((row_num >= 1), "SpRow::row(): invalid row index");
  
  return SpMat<eT>::begin();
  }



template<typename eT>
inline
typename SpRow<eT>::const_row_iterator
SpRow<eT>::begin_row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  // Since this is a row, row_num can only be 0.  But the option is provided for
  // compatibility.
  arma_debug_check((row_num >= 1), "SpRow::row(): invalid row index");
  
  return SpMat<eT>::begin();
  }



template<typename eT>
inline
typename SpRow<eT>::row_iterator
SpRow<eT>::end_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  // Since this is a row, row_num can only be 0.  But the option is provided for
  // compatibility.
  arma_debug_check((row_num >= 1), "SpRow::row(): invalid row index");
  
  return SpMat<eT>::end();
  }



template<typename eT>
inline
typename SpRow<eT>::const_row_iterator
SpRow<eT>::end_row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  // Since this is a row, row_num can only be 0.  But the option is provided for
  // compatibility.
  arma_debug_check((row_num >= 1), "SpRow::row(): invalid row index");
  
  return SpMat<eT>::end();
  }



  
#ifdef ARMA_EXTRA_SPROW_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_SPROW_MEAT)
#endif



//! @}
