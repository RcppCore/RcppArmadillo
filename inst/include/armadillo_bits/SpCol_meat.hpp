// Copyright (C) 2011-2012 Ryan Curtin <ryan@igglybob.com>
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


//! \addtogroup SpCol
//! @{


//! construct an empty column vector
template<typename eT>
inline
SpCol<eT>::SpCol()
  : SpMat<eT>(0, 1)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 1;
  }



//! construct a column vector with the specified number of n_elem
template<typename eT>
inline
SpCol<eT>::SpCol(const uword in_n_elem)
  : SpMat<eT>(in_n_elem, 1)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 1;
  }



template<typename eT>
inline
SpCol<eT>::SpCol(const uword in_n_rows, const uword in_n_cols)
  : SpMat<eT>(in_n_rows, in_n_cols)
  {
  arma_extra_debug_sigprint();

  arma_debug_check((in_n_cols != 1), "SpCol::SpCol(): must have only one column");

  access::rw(SpMat<eT>::vec_state) = 1;
  }



//! construct a column vector from specified text
template<typename eT>
inline
SpCol<eT>::SpCol(const char* text)
  : SpMat<eT>(text)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 1;

  arma_debug_check((SpMat<eT>::n_cols != 1), "SpCol::SpCol(): must have only one column");
  }



//! construct a column vector from specified text
template<typename eT>
inline
const SpCol<eT>&
SpCol<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();

  SpMat<eT>::init(std::string(text));

  access::rw(SpMat<eT>::vec_state) = 1;

  return *this;
  }



//! construct a column vector from specified text
template<typename eT>
inline
SpCol<eT>::SpCol(const std::string& text)
  : SpMat<eT>(text)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 1;

  arma_debug_check((SpMat<eT>::n_cols != 1), "SpCol::SpCol(): must have only one column");
  }



//! construct a column vector from specified text
template<typename eT>
inline
const SpCol<eT>&
SpCol<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();

  SpMat<eT>::init(std::string(text));

  return *this;
  }



template<typename eT>
inline
const SpCol<eT>&
SpCol<eT>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();

  SpMat<eT>::operator=(val);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpCol<eT>::SpCol(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 1;

  SpMat<eT>::operator=(X.get_ref());
  }



template<typename eT>
template<typename T1>
inline
const SpCol<eT>&
SpCol<eT>::operator=(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 1;

  SpMat<eT>::operator=(X.get_ref());

  return *this;
  }



template<typename eT>
template<typename T1, typename T2>
inline
SpCol<eT>::SpCol
  (
  const Base<typename SpCol<eT>::pod_type, T1>& A,
  const Base<typename SpCol<eT>::pod_type, T2>& B
  )
  {
  arma_extra_debug_sigprint();

  access::rw(SpMat<eT>::vec_state) = 1;

  SpMat<eT>::init(A,B);
  }



template<typename eT>
arma_inline
SpValProxy<SpMat<eT> >&
SpCol<eT>::row(const uword row_num)
  {
  arma_debug_check( (row_num >= SpMat<eT>::n_rows), "SpCol::row(): out of bounds" );

  return SpMat<eT>::at(row_num, 0);
  }



template<typename eT>
arma_inline
eT
SpCol<eT>::row(const uword row_num) const
  {
  arma_debug_check( (row_num >= SpMat<eT>::n_rows), "SpCol::row(): out of bounds" );
  
  return SpMat<eT>::at(row_num, 0);
  }


/*
template<typename eT>
arma_inline
subview_col<eT>
SpCol<eT>::rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= SpMat<eT>::n_rows) ), "Col::rows(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return subview_col<eT>(*this, 0, in_row1, subview_n_rows);
  }



template<typename eT>
arma_inline
const subview_col<eT>
SpCol<eT>::rows(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= SpMat<eT>::n_rows) ), "Col::rows(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return subview_col<eT>(*this, 0, in_row1, subview_n_rows);
  }



template<typename eT>
arma_inline
subview_col<eT>
SpCol<eT>::subvec(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= SpMat<eT>::n_rows) ), "Col::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return subview_col<eT>(*this, 0, in_row1, subview_n_rows);
  }



template<typename eT>
arma_inline
const subview_col<eT>
SpCol<eT>::subvec(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= SpMat<eT>::n_rows) ), "Col::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return subview_col<eT>(*this, 0, in_row1, subview_n_rows);
  }



template<typename eT>
arma_inline
subview_col<eT>
SpCol<eT>::subvec(const span& row_span)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;

  const uword local_n_rows = SpMat<eT>::n_rows;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword subvec_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;

  arma_debug_check( ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) ), "Col::subvec(): indices out of bounds or incorrectly used");
  
  return subview_col<eT>(*this, 0, in_row1, subvec_n_rows);
  }



template<typename eT>
arma_inline
const subview_col<eT>
SpCol<eT>::subvec(const span& row_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;

  const uword local_n_rows = SpMat<eT>::n_rows;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword subvec_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;

  arma_debug_check( ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) ), "Col::subvec(): indices out of bounds or incorrectly used");
  
  return subview_col<eT>(*this, 0, in_row1, subvec_n_rows);
  }
*/


//! remove specified row
template<typename eT>
inline
void
SpCol<eT>::shed_row(const uword row_num)
  {
  arma_extra_debug_sigprint();

  arma_debug_check( row_num >= SpMat<eT>::n_rows, "Col::shed_row(): out of bounds");
  
  shed_rows(row_num, row_num);
  }



//! remove specified rows
template<typename eT>
inline
void
SpCol<eT>::shed_rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();

  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= SpMat<eT>::n_rows),
    "Col::shed_rows(): indices out of bounds or incorrectly used"
    );
  
  const uword diff = (in_row2 - in_row1 + 1);

  // This is easy because everything is in one column.
  uword start = 0, end = 0;
  bool start_found = false, end_found = false;
  for(uword i = 0; i < SpMat<eT>::n_nonzero; ++i)
    {
    // Start position found?
    if (SpMat<eT>::row_indices[i] >= in_row1 && !start_found)
      {
      start = i;
      start_found = true;
      }

    // End position found?
    if (SpMat<eT>::row_indices[i] > in_row2)
      {
      end = i;
      end_found = true;
      break;
      }
    }

  if (!end_found)
    {
    end = SpMat<eT>::n_nonzero;
    }

  // Now we can make the copy.
  if (start != end)
    {
    const uword elem_diff = end - start;

    eT*    new_values      = memory::acquire_chunked<eT>   (SpMat<eT>::n_nonzero - elem_diff);
    uword* new_row_indices = memory::acquire_chunked<uword>(SpMat<eT>::n_nonzero - elem_diff);

    // Copy before the section we are dropping (if it exists).
    if (start > 0)
      {
      arrayops::copy(new_values, SpMat<eT>::values, start);
      arrayops::copy(new_row_indices, SpMat<eT>::row_indices, start);
      }

    // Copy after the section we are dropping (if it exists).
    if (end != SpMat<eT>::n_nonzero)
      {
      arrayops::copy(new_values + start, SpMat<eT>::values + end, (SpMat<eT>::n_nonzero - end));
      arrayops::copy(new_row_indices + start, SpMat<eT>::row_indices + end, (SpMat<eT>::n_nonzero - end));
      arrayops::inplace_minus(new_row_indices + start, diff, (SpMat<eT>::n_nonzero - end));
      }

    memory::release(SpMat<eT>::values);
    memory::release(SpMat<eT>::row_indices);

    access::rw(SpMat<eT>::values) = new_values;
    access::rw(SpMat<eT>::row_indices) = new_row_indices;

    access::rw(SpMat<eT>::n_nonzero) -= elem_diff;
    access::rw(SpMat<eT>::col_ptrs[1]) -= elem_diff;
    }

  access::rw(SpMat<eT>::n_rows) -= diff;
  access::rw(SpMat<eT>::n_elem) -= diff;

  }



//! insert N rows at the specified row position,
//! optionally setting the elements of the inserted rows to zero
template<typename eT>
inline
void
SpCol<eT>::insert_rows(const uword row_num, const uword N, const bool set_to_zero)
  {
  arma_extra_debug_sigprint();

  arma_debug_check(set_to_zero == false, "SpCol::insert_rows(): cannot set nonzero values");

  arma_debug_check((row_num > SpMat<eT>::n_rows), "SpCol::insert_rows(): out of bounds");

  for(uword row = 0; row < SpMat<eT>::n_rows; ++row)
    {
    if (SpMat<eT>::row_indices[row] >= row_num)
      {
      access::rw(SpMat<eT>::row_indices[row]) += N;
      }
    }

  access::rw(SpMat<eT>::n_rows) += N;
  access::rw(SpMat<eT>::n_elem) += N;
  }



//! insert the given object at the specified row position;
//! the given object must have one column
template<typename eT>
template<typename T1>
inline
void
SpCol<eT>::insert_rows(const uword row_num, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();

  SpMat<eT>::insert_rows(row_num, X);
  }



template<typename eT>
inline
typename SpCol<eT>::row_iterator
SpCol<eT>::begin_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= SpMat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return row_iterator(*this, row_num, 0);
  }



template<typename eT>
inline
typename SpCol<eT>::const_row_iterator
SpCol<eT>::begin_row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= SpMat<eT>::n_rows), "begin_row(): index out of bounds");
  
  return const_row_iterator(*this, row_num, 0);
  }



template<typename eT>
inline
typename SpCol<eT>::row_iterator
SpCol<eT>::end_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= SpMat<eT>::n_rows), "end_row(): index out of bounds");
  
  return row_iterator(*this, row_num + 1, 0);
  }



template<typename eT>
inline
typename SpCol<eT>::const_row_iterator
SpCol<eT>::end_row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= SpMat<eT>::n_rows), "end_row(): index out of bounds");
  
  return const_row_iterator(*this, row_num + 1, 0);
  }


#ifdef ARMA_EXTRA_COL_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_COL_MEAT)
#endif



//! @}
