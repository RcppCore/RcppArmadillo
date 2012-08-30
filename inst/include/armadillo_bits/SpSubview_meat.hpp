// Copyright (C) 2011-2012 Ryan Curtin <ryan@igglybob.com>
// Copyright (C) 2011 Matthew Amidon
// Copyright (C) 2012 Conrad Sanderson
//
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose.  You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



template<typename eT>
arma_inline
SpSubview<eT>::SpSubview(const SpMat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols)
  : m(in_m)
  , aux_row1(in_row1)
  , aux_col1(in_col1)
  , n_rows(in_n_rows)
  , n_cols(in_n_cols)
  , n_elem(in_n_rows * in_n_cols)
  , n_nonzero(0)
  {
  arma_extra_debug_sigprint();

  // There must be a O(1) way to do this
  uword end     = m.col_ptrs[in_col1 + in_n_cols];
  uword end_row = in_row1 + in_n_rows;
  uword count   = 0;

  for(uword i = m.col_ptrs[in_col1]; i < end; ++i)
    {
    if(m.row_indices[i] >= in_row1 && m.row_indices[i] < end_row)
      {
      ++count;
      }
    }

  access::rw(n_nonzero) = count;
  }



template<typename eT>
arma_inline
SpSubview<eT>::SpSubview(SpMat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols)
  : m(in_m)
  , aux_row1(in_row1)
  , aux_col1(in_col1)
  , n_rows(in_n_rows)
  , n_cols(in_n_cols)
  , n_elem(in_n_rows * in_n_cols)
  , n_nonzero(0)
  {
  arma_extra_debug_sigprint();

  // There must be a O(1) way to do this
  uword end     = m.col_ptrs[in_col1 + in_n_cols];
  uword end_row = in_row1 + in_n_rows;
  uword count   = 0;

  for(uword i = m.col_ptrs[in_col1]; i < end; ++i)
    {
    if(m.row_indices[i] >= in_row1 && m.row_indices[i] < end_row)
      {
      ++count;
      }
    }

  access::rw(n_nonzero) = count;
  }



template<typename eT>
inline
SpSubview<eT>::~SpSubview()
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator+=(const eT val)
  {
  arma_extra_debug_sigprint();

  if(val == 0)
    {
    return *this;
    }

  const uword start_row = aux_row1;
  const uword end_row   = aux_row1 + n_rows;
  const uword start_col = aux_col1;
  const uword end_col   = aux_col1 + n_cols;

  const uword old_n_nonzero = m.n_nonzero;

  // iterate over our part of the sparse matrix
  for(uword col = start_col; col < end_col; ++col)
    {
    for(uword row = start_row; row < end_row; ++row)
      {
      access::rw(m).at(row, col) += val;
      }
    }

  const uword new_n_nonzero = m.n_nonzero;

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();

  if(val == 0)
    {
    return *this;
    }

  const uword start_row = aux_row1;
  const uword end_row   = aux_row1 + n_rows;
  const uword start_col = aux_col1;
  const uword end_col   = aux_col1 + n_cols;

  const uword old_n_nonzero = m.n_nonzero;

  for(uword col = start_col; col < end_col; ++col)
    {
    for(uword row = start_row; row < end_row; ++row)
      {
      access::rw(m).at(row, col) -= val;
      }
    }

  const uword new_n_nonzero = m.n_nonzero;

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();

  if(val == eT(0))
    {
    // Turn it all into zeros.
    for(typename SpMat<eT>::iterator it(access::rw(m), aux_row1, aux_col1); it.col() < (aux_col1 + n_cols); ++it)
      {
      if((it.row() >= aux_row1) && (it.row() < (aux_row1 + n_cols)))
        {
        (*it) = eT(0); // zero out the value.
        }
      }

    return *this;
    }

  const uword start_row = aux_row1;
  const uword end_row   = aux_row1 + n_rows;
  const uword start_col = aux_col1;
  const uword end_col   = aux_col1 + n_cols;

  for(uword c = start_col; c < end_col; ++c)
    {
    for(uword r = m.col_ptrs[c]; r < m.col_ptrs[c + 1]; ++r)
      {
      if(m.row_indices[r] >= start_row && m.row_indices[r] < end_row)
        {
        access::rw(m.values[r]) *= val;
        }
      }
    }

  return *this;
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();

  const uword start_col = aux_col1;
  const uword end_col   = aux_col1 + n_cols;
  const uword start_row = aux_row1;
  const uword end_row   = aux_row1 + n_rows;

  for(uword c = start_col; c < end_col; ++c)
    {
    for(uword r = start_row; r < end_row; ++r)
      {
      access::rw(m).at(r, c) /= val;
      }
    }

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const Proxy<T1> P(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "insert into sparse submatrix");

  const uword old_n_nonzero = m.n_nonzero;

  for(uword c = 0; c < n_cols; ++c)
    {
    for(uword r = 0; r < n_rows; ++r)
      {
      at(r, c) = P.at(r, c);
      }
    }

  const uword new_n_nonzero = m.n_nonzero;

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator+=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const Proxy<T1> P(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "addition into sparse submatrix");

  const uword old_n_nonzero = m.n_nonzero;

  for(uword c = 0; c < n_cols; ++c)
    {
    for(uword r = 0; r < n_rows; ++r)
      {
      at(r, c) += P.at(r, c);
      }
    }

  const uword new_n_nonzero = m.n_nonzero;

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator-=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const Proxy<T1> P(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "subtraction into sparse submatrix");

  const uword old_n_nonzero = m.n_nonzero;

  for(uword c = 0; c < n_cols; ++c)
    {
    for(uword r = 0; r < n_rows; ++r)
      {
      at(r, c) -= P.at(r, c);
      }
    }

  const uword new_n_nonzero = m.n_nonzero;

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator*=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const Proxy<T1> P(x.get_ref());

  // Must be exactly the same size for this (we can't modify our own size).
  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "matrix multiplication into sparse submatrix");

  SpMat<eT> tmp(*this);
  Mat<eT> other_tmp(x.get_ref());
  tmp *= other_tmp;
  operator=(tmp);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator%=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const Proxy<T1> P(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "element-wise multiplication into sparse submatrix");

  const uword old_n_nonzero = m.n_nonzero;

  for(typename SpMat<eT>::iterator it(access::rw(m), aux_row1, aux_col1); it.col() < (aux_col1 + n_cols); ++it)
    {
    if((it.row() >= aux_row1) && (it.row() < (aux_row1 + n_rows)))
      {
      (*it) *= P.at(it.row() - aux_row1, it.col() - aux_col1);
      }
    }

  const uword new_n_nonzero = m.n_nonzero;

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator/=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const Proxy<T1> P(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "element-wise division into sparse submatrix");

  const uword old_n_nonzero = m.n_nonzero;

  for(uword c = 0; c < n_cols; ++c)
    {
    for(uword r = 0; r < n_rows; ++r)
      {
      at(r, c) /= P.at(r, c);
      }
    }

  const uword new_n_nonzero = m.n_nonzero;

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator=(const SpSubview<eT>& x)
  {
  arma_extra_debug_sigprint();

  arma_debug_assert_same_size(n_rows, n_cols, x.n_rows, x.n_cols, "insertion into sparse submatrix");

  const_iterator cit = x.begin();
  iterator it = begin();

  while((cit.pos() < x.n_nonzero) || (it.pos() < n_nonzero))
    {
    if((cit.row() == it.row()) && (cit.col() == it.col()))
      {
      (*it) = (*cit);
      ++it;
      ++cit;
      }
    else
      {
      if((cit.col() > it.col()) || ((cit.col() == it.col()) && (cit.row() > it.row())))
        {
        // cit is "ahead"
        (*it) = eT(0); // erase element
        it.internal_pos--; // update iterator so it still works
        ++it;
        }
      else
        {
        // it is "ahead"
        at(cit.row(), cit.col()) = (*cit);
        it.internal_pos++; // update iterator so it still works
        ++cit;
        }
      }
    }

  access::rw(n_nonzero) = x.n_nonzero;

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T1> p(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "insertion into sparse submatrix");

  typename SpProxy<T1>::const_iterator_type cit = p.begin();
  iterator it(*this);

  while((cit.pos() < p.get_n_nonzero()) || (it.pos() < n_nonzero))
    {
    if(cit == it) // at the same location
      {
      (*it) = (*cit);
      ++it;
      ++cit;
      }
    else
      {
      if((cit.col() > it.col()) || ((cit.col() == it.col()) && (cit.row() > it.row())))
        {
        // cit is "ahead"
        (*it) = eT(0); // erase element
        it.internal_pos--; // update iterator so it still works
        ++it;
        }
      else
        {
        // it is "ahead"
        at(cit.row(), cit.col()) = (*cit);
        it.internal_pos++; // update iterator so it still works
        ++cit;
        }
      }
    }

  access::rw(n_nonzero) = p.get_n_nonzero();

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator+=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T1> p(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "addition into sparse submatrix");

  typename SpProxy<T1>::const_iterator_type cit = p.begin();
  iterator it(*this);

  const uword old_n_nonzero = p.get_n_nonzero();

  while(it.pos() < n_nonzero)
    {
    if(it == cit)
      {
      const eT val = (*it) + (*cit);
      (*it) = val;
      if(val == 0)
        {
        it.internal_pos--; // to keep it valid
        }
      ++it;
      ++cit;
      }
    else
      {
      if((cit.col() > it.col()) || ((cit.col() == it.col()) && (cit.row() > it.row())))
        {
        // cit is "ahead"; we don't need to do anything
        ++it;
        }
      else
        {
        // it is "ahead"; we need to insert a new value
        at(cit.row(), cit.col()) = (*cit);
        it.internal_pos++; // update iterator so it still works
        ++cit;
        }
      }
    }

  const uword new_n_nonzero = p.get_n_nonzero();

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator-=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T1> p(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "subtraction into sparse submatrix");

  typename SpProxy<T1>::const_iterator_type cit = p.begin();
  iterator it(*this);

  const uword old_n_nonzero = p.get_n_nonzero();

  while(it.pos() < n_nonzero)
    {
    if(cit == it)
      {
      const eT val = (*it) - (*cit);
      (*it) = val;
      if(val == 0)
        {
        it.internal_pos--; // to keep it valid
        }
      ++it;
      ++cit;
      }
    else
      {
      if((cit.col() > it.col()) || ((cit.col() == it.col()) && (cit.row() > it.row())))
        {
        // cit is "ahead"; we don't need to do anything
        ++it;
        }
      else
        {
        // it is "ahead"; we need to insert a new value
        at(cit.row(), cit.col()) = -(*cit);
        it.internal_pos++; // update iterator so it still works
        ++cit;
        }
      }
    }

  const uword new_n_nonzero = p.get_n_nonzero();

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator*=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T1> p(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "matrix multiplication into sparse submatrix");

  // Because we have to use a temporary anyway, it doesn't make sense to
  // reimplement this here.
  return operator=((*this) * x.get_ref());
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator%=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T1> p(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "element-wise multiplication into sparse submatrix");

  typename SpProxy<T1>::const_iterator_type cit = p.begin();
  iterator it(*this);

  const uword old_n_nonzero = p.get_n_nonzero();

  while(it.pos() < n_nonzero)
    {
    if(it == cit)
      {
      (*it) *= (*cit);
      ++it;
      ++cit;
      }
    else
      {
      if((cit.col() > it.col()) || ((cit.col() == it.col()) && (cit.row() > it.row())))
        {
        // cit is "ahead"; this value becomes zero.
        (*it) = eT(0);
        it.internal_pos--; // keep it accurate
        ++it;
        }
      else
        {
        // it is "ahead"; we just need to catch cit up
        ++cit;
        }
      }
    }

  const uword new_n_nonzero = p.get_n_nonzero();

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator/=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();

  SpProxy<T1> p(x.get_ref());

  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "element-wise division into sparse submatrix");

  const uword old_n_nonzero = p.get_n_nonzero();

  // If you are using this function, you are probably misguided.
  for(uword i = 0; i < n_elem; ++i)
    {
    at(i) /= p[i];
    }

  const uword new_n_nonzero = p.get_n_nonzero();

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

  return *this;
  }



template<typename eT>
inline
void
SpSubview<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();

  (*this).operator=(val);
  }



template<typename eT>
inline
void
SpSubview<eT>::zeros()
  {
  arma_extra_debug_sigprint();

  // we can be a little smarter here
  for(typename SpMat<eT>::iterator it(access::rw(m), aux_row1, aux_col1);
      (it.col < (aux_col1 + n_cols - 1)) || (it.col == (aux_col1 + n_cols - 1) && it.row < (aux_row1 + n_rows));
      ++it)
    {
    // column will always be valid; no need to check that
    if((it.row >= aux_row1) && (it.row < (aux_row1 + n_rows)))
      {
      (*it) = 0;
      }
    }

  access::rw(n_nonzero) = 0;
  }



template<typename eT>
inline
void
SpSubview<eT>::ones()
  {
  arma_extra_debug_sigprint();

  (*this).fill(eT(1));
  }



template<typename eT>
inline
void
SpSubview<eT>::eye()
  {
  arma_extra_debug_sigprint();

  // clear other things
  (*this).zeros();

  // now the diagonal ones
  const uword end_index = std::min(n_rows, n_cols);

  for(uword ind = 0; ind < end_index; ++ind)
    {
    m.at(ind + aux_row1, ind + aux_col1) = eT(1);
    }

  access::rw(n_nonzero) = end_index;
  }



template<typename eT>
inline
SpValProxy<SpSubview<eT> >
SpSubview<eT>::operator[](const uword i)
  {
  const uword row = i % n_rows;
  const uword col = i / n_rows;

  return (*this).at(row, col);
  }



template<typename eT>
inline
eT
SpSubview<eT>::operator[](const uword i) const
  {
  const uword row = i % n_rows;
  const uword col = i / n_rows;

  return (*this).at(row, col);
  }



template<typename eT>
inline
SpValProxy<SpSubview<eT> >
SpSubview<eT>::operator()(const uword i)
  {
  arma_debug_check( (i >= n_elem), "SpSubview::operator(): index out of bounds");

  const uword row = i % n_rows;
  const uword col = i / n_cols;

  return (*this).at(row, col);
  }



template<typename eT>
inline
eT
SpSubview<eT>::operator()(const uword i) const
  {
  arma_debug_check( (i >= n_elem), "SpSubview::operator(): index out of bounds");

  const uword row = i % n_rows;
  const uword col = i / n_cols;

  return (*this).at(row, col);
  }



template<typename eT>
inline
SpValProxy<SpSubview<eT> >
SpSubview<eT>::operator()(const uword in_row, const uword in_col)
  {
  arma_debug_check( (in_row >= n_rows) || (in_col >= n_cols), "SpSubview::operator(): index out of bounds");

  return (*this).at(in_row, in_col);
  }



template<typename eT>
inline
eT
SpSubview<eT>::operator()(const uword in_row, const uword in_col) const
  {
  arma_debug_check( (in_row >= n_rows) || (in_col >= n_cols), "SpSubview::operator(): index out of bounds");

  return (*this).at(in_row, in_col);
  }



template<typename eT>
inline
SpValProxy<SpSubview<eT> >
SpSubview<eT>::at(const uword i)
  {
  const uword row = i % n_rows;
  const uword col = i / n_cols;

  return (*this).at(row, col);
  }



template<typename eT>
inline
eT
SpSubview<eT>::at(const uword i) const
  {
  const uword row = i % n_rows;
  const uword col = i / n_cols;

  return (*this).at(row, col);
  }



template<typename eT>
inline
SpValProxy<SpSubview<eT> >
SpSubview<eT>::at(const uword in_row, const uword in_col)
  {
  const uword colptr      = m.col_ptrs[in_col + aux_col1];
  const uword next_colptr = m.col_ptrs[in_col + aux_col1 + 1];

  // Step through the row indices to see if our element exists.
  for(uword i = colptr; i < next_colptr; ++i)
    {
    // First check that we have not stepped past it.
    if((in_row + aux_row1) < m.row_indices[i])
      {
      return SpValProxy<SpSubview<eT> >(in_row, in_col, *this); // Proxy for a zero value.
      }

    // Now check if we are at the correct place.
    if((in_row + aux_row1) == m.row_indices[i]) // If we are, return a reference to the value.
      {
      return SpValProxy<SpSubview<eT> >(in_row, in_col, *this, &access::rw(m.values[i]));
      }
    }

  // We did not find it, so it does not exist.
  return SpValProxy<SpSubview<eT> >(in_row, in_col, *this);
  }



template<typename eT>
inline
eT
SpSubview<eT>::at(const uword in_row, const uword in_col) const
  {
  return m.at(aux_row1 + in_row, aux_col1 + in_col);
  }



template<typename eT>
inline
bool
SpSubview<eT>::check_overlap(const SpSubview<eT>& x) const
  {
  const subview<eT>& t = *this;

  if(&t.m != &x.m)
    {
    return false;
    }
  else
    {
    if( (t.n_elem == 0) || (x.n_elem == 0) )
      {
      return false;
      }
    else
      {
      const uword t_row_start  = t.aux_row1;
      const uword t_row_end_p1 = t_row_start + t.n_rows;

      const uword t_col_start  = t.aux_col1;
      const uword t_col_end_p1 = t_col_start + t.n_cols;

      const uword x_row_start  = x.aux_row1;
      const uword x_row_end_p1 = x_row_start + x.n_rows;

      const uword x_col_start  = x.aux_col1;
      const uword x_col_end_p1 = x_col_start + x.n_cols;

      const bool outside_rows = ( (x_row_start >= t_row_end_p1) || (t_row_start >= x_row_end_p1) );
      const bool outside_cols = ( (x_col_start >= t_col_end_p1) || (t_col_start >= x_col_end_p1) );

      return ( (outside_rows == false) && (outside_cols == false) );
      }
    }
  }



template<typename eT>
inline
bool
SpSubview<eT>::is_vec() const
  {
  return ( (n_rows == 1) || (n_cols == 1) );
  }



template<typename eT>
inline
void
SpSubview<eT>::swap_rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();

  arma_debug_check((in_row1 >= n_rows) || (in_row2 >= n_rows), "SpSubview::swap_rows(): invalid row index");

  const uword start_col = aux_col1;
  const uword end_col   = aux_col1 + n_cols;

  for(uword c = start_col; c < end_col; ++c)
    {
    eT val = m.at(in_row1 + aux_row1, c);
    access::rw(m).at(in_row2 + aux_row1, c) = m.at(in_row1 + aux_row1, c);
    access::rw(m).at(in_row1 + aux_row1, c) = val;
    }
  }



template<typename eT>
inline
void
SpSubview<eT>::swap_cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();

  arma_debug_check((in_col1 >= n_cols) || (in_col2 >= n_cols), "SpSubview::swap_cols(): invalid column index");

  const uword start_row = aux_row1;
  const uword end_row   = aux_row1 + n_rows;

  for(uword r = start_row; r < end_row; ++r)
    {
    eT val = m.at(r, in_col1 + aux_col1);
    access::rw(m).at(r, in_col1 + aux_col1) = m.at(r, in_col2 + aux_col1);
    access::rw(m).at(r, in_col2 + aux_col1) = val;
    }
  }



template<typename eT>
inline
typename SpSubview<eT>::iterator
SpSubview<eT>::begin()
  {
  return iterator(*this);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::begin() const
  {
  return const_iterator(*this);
  }



template<typename eT>
inline
typename SpSubview<eT>::iterator
SpSubview<eT>::begin_col(const uword col_num)
  {
  return iterator(*this, 0, col_num);
  }


template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::begin_col(const uword col_num) const
  {
  return const_iterator(*this, 0, col_num);
  }



template<typename eT>
inline
typename SpSubview<eT>::row_iterator
SpSubview<eT>::begin_row()
  {
  return row_iterator(*this);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::begin_row() const
  {
  return const_row_iterator(*this);
  }



template<typename eT>
inline
typename SpSubview<eT>::iterator
SpSubview<eT>::end()
  {
  return iterator(*this, n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::end() const
  {
  return const_iterator(*this, n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::row_iterator
SpSubview<eT>::end_row()
  {
  return row_iterator(*this, n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::end_row() const
  {
  return const_row_iterator(*this, n_nonzero);
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
eT&
SpSubview<eT>::add_element(const uword in_row, const uword in_col, const eT in_val)
  {
  arma_extra_debug_sigprint();

  // This may not actually add an element.
  const uword old_n_nonzero = m.n_nonzero;
  eT& retval = access::rw(m).add_element(in_row + aux_row1, in_col + aux_col1, in_val);
  // Update n_nonzero (if necessary).
  access::rw(n_nonzero) += (m.n_nonzero - old_n_nonzero);

  return retval;
  }



template<typename eT>
inline
void
SpSubview<eT>::delete_element(const uword in_row, const uword in_col)
  {
  arma_extra_debug_sigprint();

  // This may not actually delete an element.
  const uword old_n_nonzero = m.n_nonzero;
  access::rw(m).delete_element(in_row + aux_row1, in_col + aux_col1);
  access::rw(n_nonzero) -= (old_n_nonzero - m.n_nonzero);
  }



/**
 * Sparse subview col
 *
template<typename eT>
inline
SpSubview_col<eT>::SpSubview_col(const Mat<eT>& in_m, const uword in_col)
  {
  arma_extra_debug_sigprint();
  }

template<typename eT>
inline
SpSubview_col<eT>::SpSubview_col(Mat<eT>& in_m, const uword in_col)
  {
  arma_extra_debug_sigprint();
  }

template<typename eT>
inline
SpSubview_col<eT>::SpSubview_col(const Mat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows)
  {
  arma_extra_debug_sigprint();
  }

template<typename eT>
inline
SpSubview_col<eT>::SpSubview_col(Mat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows)
  {
  arma_extra_debug_sigprint();
  }
*/

/**
 * Sparse subview row
 *
template<typename eT>
inline
SpSubview_row<eT>::SpSubview_row(const Mat<eT>& in_m, const uword in_row)
  {
  arma_extra_debug_sigprint();
  }

template<typename eT>
inline
SpSubview_row<eT>::SpSubview_row(Mat<eT>& in_m, const uword in_row)
  {
  arma_extra_debug_sigprint();
  }

template<typename eT>
inline
SpSubview_row<eT>::SpSubview_row(const Mat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  }

template<typename eT>
inline
SpSubview_row<eT>::SpSubview_row(Mat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  }
*/
