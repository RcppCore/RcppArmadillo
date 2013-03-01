// Copyright (C) 2011-2012 Ryan Curtin
// Copyright (C) 2011 Matthew Amidon
// Copyright (C) 2012 Conrad Sanderson
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup SpSubview
//! @{

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
  uword lend     = m.col_ptrs[in_col1 + in_n_cols];
  uword lend_row = in_row1 + in_n_rows;
  uword count   = 0;

  for(uword i = m.col_ptrs[in_col1]; i < lend; ++i)
    {
    if(m.row_indices[i] >= in_row1 && m.row_indices[i] < lend_row)
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
  uword lend     = m.col_ptrs[in_col1 + in_n_cols];
  uword lend_row = in_row1 + in_n_rows;
  uword count    = 0;

  for(uword i = m.col_ptrs[in_col1]; i < lend; ++i)
    {
    if(m.row_indices[i] >= in_row1 && m.row_indices[i] < lend_row)
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

  if(val == eT(0))
    {
    return *this;
    }

  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;

  const uword old_n_nonzero = m.n_nonzero;

  // iterate over our part of the sparse matrix
  for(uword lcol = lstart_col; lcol < lend_col; ++lcol)
  for(uword lrow = lstart_row; lrow < lend_row; ++lrow)
    {
    access::rw(m).at(lrow, lcol) += val;
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

  if(val == eT(0))
    {
    return *this;
    }

  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;

  const uword old_n_nonzero = m.n_nonzero;

  for(uword lcol = lstart_col; lcol < lend_col; ++lcol)
  for(uword lrow = lstart_row; lrow < lend_row; ++lrow)
    {
    access::rw(m).at(lrow, lcol) -= val;
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
    for(iterator it(*this); it != end(); ++it)
      {
      (*it) = eT(0); // zero out the value.
      it.internal_pos--;
      }

    return *this;
    }

  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;

  for(uword c = lstart_col; c < lend_col; ++c)
    {
    for(uword r = m.col_ptrs[c]; r < m.col_ptrs[c + 1]; ++r)
      {
      if(m.row_indices[r] >= lstart_row && m.row_indices[r] < lend_row)
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

  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;
  
  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;

  const uword old_n_nonzero = m.n_nonzero;

  for(uword c = lstart_col; c < lend_col; ++c)
  for(uword r = lstart_row; r < lend_row; ++r)
    {
    access::rw(m).at(r, c) /= val;
    }

  const uword new_n_nonzero = m.n_nonzero;

  access::rw(n_nonzero) += (new_n_nonzero - old_n_nonzero);

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

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "addition");

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

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "subtraction");

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
  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "matrix multiplication");

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

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "element-wise multiplication");

  const uword old_n_nonzero = m.n_nonzero;

  for(iterator it(*this); it != end(); ++it)
    {
    (*it) *= P.at(it.row(), it.col());
    if(P.at(it.row(), it.col()) == eT(0))
      {
      it.internal_pos--;
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

  arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "element-wise division");

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
  
  const bool alias = ( &m == &(x.m) );
  
  if(alias == false)
    {
    const_iterator cit = x.begin();
    iterator        it = begin();
    
    while((cit != x.end()) || (it != end()))
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
    }
  else
    {
    const SpMat<eT> tmp(x);
    
    (*this).operator=(tmp);
    }
  
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
  
  if(p.is_alias(m) == false)
    {
    typename SpProxy<T1>::const_iterator_type cit = p.begin();
    iterator it(*this);

    while((cit != p.end()) || (it != end()))
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
    }
  else
    {
    const SpMat<eT> tmp(p.Q);
    
    (*this).operator=(tmp);
    }
  
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
  
  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "addition");
  
  if(p.is_alias(m) == false)
    {
    typename SpProxy<T1>::const_iterator_type cit = p.begin();

    while(cit != p.end())
      {
      at(cit.row(), cit.col()) += (*cit);
      ++cit;
      }
    }
  else
    {
    const SpMat<eT> tmp(p.Q);
    
    (*this).operator+=(tmp);
    }
  
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
  
  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "subtraction");
  
  if(p.is_alias(m) == false)
    {
    typename SpProxy<T1>::const_iterator_type cit = p.begin();
    
    while(cit != p.end())
      {
      at(cit.row(), cit.col()) -= (*cit);
      ++cit;
      }
    }
  else
    {
    const SpMat<eT> tmp(p.Q);
    
    (*this).operator-=(tmp);
    }
  
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

  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "matrix multiplication");

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
  
  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "element-wise multiplication");
  
  if(p.is_alias(m) == false)
    {
    typename SpProxy<T1>::const_iterator_type cit = p.begin();
    iterator it(*this);

    while((it != end()) || (cit != p.end()))
      {
      if((cit.row() == it.row()) && (cit.col() == it.col()))
        {
        (*it) *= (*cit);
        ++it;
        ++cit;
        }
      else
        {
        if((cit.col() > it.col()) || ((cit.col() == it.col()) && (cit.row() > it.row())))
          {
          // cit is "ahead"
          (*it) = eT(0); // erase element -- x has a zero here
          it.internal_pos--; // update iterator so it still works
          ++it;
          }
        else
          {
          // it is "ahead"
          ++cit;
          }
        }
      }
    }
  else
    {
    const SpMat<eT> tmp(p.Q);
    
    (*this).operator%=(tmp);
    }
  
  return *this;
  }



//! If you are using this function, you are probably misguided.
template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator/=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  SpProxy<T1> p(x.get_ref());
  
  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "element-wise division");
  
  if(p.is_alias(m) == false)
    {
    for(uword lcol = 0; lcol < n_cols; ++lcol)
    for(uword lrow = 0; lrow < n_rows; ++lrow)
      {
      at(lrow,lcol) /= p.at(lrow,lcol);
      }
    }
  else
    {
    const SpMat<eT> tmp(p.Q);
    
    (*this).operator/=(tmp);
    }
  
  return *this;
  }



template<typename eT>
inline
void
SpSubview<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val != eT(0))
    {
    // TODO: implement a faster version; the code below is slow
    
    const uword lstart_row = aux_row1;
    const uword lend_row   = aux_row1 + n_rows;
    
    const uword lstart_col = aux_col1;
    const uword lend_col   = aux_col1 + n_cols;

    const uword orig_nonzero = m.n_nonzero;
    
    // iterate over our part of the sparse matrix
    for(uword lcol = lstart_col; lcol < lend_col; ++lcol)
    for(uword lrow = lstart_row; lrow < lend_row; ++lrow)
      {
      access::rw(m).at(lrow, lcol) = val;
      }

    access::rw(n_nonzero) += (m.n_nonzero - orig_nonzero);

    }
  else
    {
    (*this).zeros();
    }
  }



template<typename eT>
inline
void
SpSubview<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  // we can be a little smarter here
  iterator it(*this);

  while(it != end())
    {
    (*it) = eT(0);
    it.internal_pos--; // hack to update iterator without requiring a new one
    ++it;
    }
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
  
  const uword orig_nonzero = m.n_nonzero;
  
  // now the diagonal ones
  const uword end_index = std::min(n_rows, n_cols);
  
  for(uword ind = 0; ind < end_index; ++ind)
    {
    access::rw(m).at(ind + aux_row1, ind + aux_col1) = eT(1);
    }
  
  access::rw(n_nonzero) += (m.n_nonzero - orig_nonzero);
  }



template<typename eT>
arma_hot
inline
SpValProxy<SpSubview<eT> >
SpSubview<eT>::operator[](const uword i)
  {
  const uword lrow = i % n_rows;
  const uword lcol = i / n_rows;

  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::operator[](const uword i) const
  {
  const uword lrow = i % n_rows;
  const uword lcol = i / n_rows;

  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
SpValProxy< SpSubview<eT> >
SpSubview<eT>::operator()(const uword i)
  {
  arma_debug_check( (i >= n_elem), "SpSubview::operator(): index out of bounds");

  const uword lrow = i % n_rows;
  const uword lcol = i / n_rows;

  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::operator()(const uword i) const
  {
  arma_debug_check( (i >= n_elem), "SpSubview::operator(): index out of bounds");

  const uword lrow = i % n_rows;
  const uword lcol = i / n_rows;

  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
SpValProxy< SpSubview<eT> >
SpSubview<eT>::operator()(const uword in_row, const uword in_col)
  {
  arma_debug_check( (in_row >= n_rows) || (in_col >= n_cols), "SpSubview::operator(): index out of bounds");

  return (*this).at(in_row, in_col);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::operator()(const uword in_row, const uword in_col) const
  {
  arma_debug_check( (in_row >= n_rows) || (in_col >= n_cols), "SpSubview::operator(): index out of bounds");

  return (*this).at(in_row, in_col);
  }



template<typename eT>
arma_hot
inline
SpValProxy< SpSubview<eT> >
SpSubview<eT>::at(const uword i)
  {
  const uword lrow = i % n_rows;
  const uword lcol = i / n_cols;

  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::at(const uword i) const
  {
  const uword lrow = i % n_rows;
  const uword lcol = i / n_cols;

  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
SpValProxy< SpSubview<eT> >
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
arma_hot
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
SpSubview<eT>
SpSubview<eT>::row(const uword row_num)
  {
  arma_extra_debug_sigprint();

  arma_debug_check(row_num >= n_rows, "SpSubview::row(): out of bounds");

  return submat(row_num, 0, row_num, n_cols - 1);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::row(const uword row_num) const
  {
  arma_extra_debug_sigprint();

  arma_debug_check(row_num >= n_rows, "SpSubview::row(): out of bounds");

  return submat(row_num, 0, row_num, n_cols - 1);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::col(const uword col_num)
  {
  arma_extra_debug_sigprint();

  arma_debug_check(col_num >= n_cols, "SpSubview::col(): out of bounds");

  return submat(0, col_num, n_rows - 1, col_num);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::col(const uword col_num) const
  {
  arma_extra_debug_sigprint();

  arma_debug_check(col_num >= n_cols, "SpSubview::col(): out of bounds");

  return submat(0, col_num, n_rows - 1, col_num);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();

  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "SpSubview::rows(): indices out of bounds or incorrectly used"
    );

  return submat(in_row1, 0, in_row2, n_cols - 1);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::rows(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();

  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "SpSubview::rows(): indices out of bounds or incorrectly used"
    );

  return submat(in_row1, 0, in_row2, n_cols - 1);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();

  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "SpSubview::cols(): indices out of bounds or incorrectly used"
    );

  return submat(0, in_col1, n_rows - 1, in_col2);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::cols(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();

  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "SpSubview::cols(): indices out of bounds or incorrectly used"
    );

  return submat(0, in_col1, n_rows - 1, in_col2);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2)
  {
  arma_extra_debug_sigprint();

  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "SpSubview::submat(): indices out of bounds or incorrectly used"
    );

  return access::rw(m).submat(in_row1 + aux_row1, in_col1 + aux_col1, in_row2 + aux_row1, in_col2 + aux_col1);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const
  {
  arma_extra_debug_sigprint();

  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "SpSubview::submat(): indices out of bounds or incorrectly used"
    );

  return m.submat(in_row1 + aux_row1, in_col1 + aux_col1, in_row2 + aux_row1, in_col2 + aux_col1);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::submat(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();

  const bool row_all = row_span.whole;
  const bool col_all = row_span.whole;

  const uword in_row1 = row_all ? 0      : row_span.a;
  const uword in_row2 = row_all ? n_rows : row_span.b;

  const uword in_col1 = col_all ? 0      : col_span.a;
  const uword in_col2 = col_all ? n_cols : col_span.b;

  arma_debug_check
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= n_rows)))
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= n_cols))),
    "SpSubview::submat(): indices out of bounds or incorrectly used"
    );

  return submat(in_row1, in_col1, in_row2, in_col2);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::submat(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();

  const bool row_all = row_span.whole;
  const bool col_all = row_span.whole;

  const uword in_row1 = row_all ? 0          : row_span.a;
  const uword in_row2 = row_all ? n_rows - 1 : row_span.b;

  const uword in_col1 = col_all ? 0          : col_span.a;
  const uword in_col2 = col_all ? n_cols - 1 : col_span.b;

  arma_debug_check
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= n_rows)))
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= n_cols))),
    "SpSubview::submat(): indices out of bounds or incorrectly used"
    );

  return submat(in_row1, in_col1, in_row2, in_col2);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::operator()(const uword row_num, const span& col_span)
  {
  arma_extra_debug_sigprint();

  return submat(span(row_num, row_num), col_span);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::operator()(const uword row_num, const span& col_span) const
  {
  arma_extra_debug_sigprint();

  return submat(span(row_num, row_num), col_span);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::operator()(const span& row_span, const uword col_num)
  {
  arma_extra_debug_sigprint();

  return submat(row_span, span(col_num, col_num));
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::operator()(const span& row_span, const uword col_num) const
  {
  arma_extra_debug_sigprint();

  return submat(row_span, span(col_num, col_num));
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::operator()(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();

  return submat(row_span, col_span);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::operator()(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();

  return submat(row_span, col_span);
  }



template<typename eT>
inline
void
SpSubview<eT>::swap_rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();

  arma_debug_check((in_row1 >= n_rows) || (in_row2 >= n_rows), "SpSubview::swap_rows(): invalid row index");

  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;

  for(uword c = lstart_col; c < lend_col; ++c)
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

  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;

  for(uword r = lstart_row; r < lend_row; ++r)
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
SpSubview<eT>::begin_row(const uword row_num)
  {
  return row_iterator(*this, row_num, 0);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::begin_row(const uword row_num) const
  {
  return const_row_iterator(*this, row_num, 0);
  }



template<typename eT>
inline
typename SpSubview<eT>::iterator
SpSubview<eT>::end()
  {
  return iterator(*this, 0, n_cols, n_nonzero, m.n_nonzero - n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::end() const
  {
  return const_iterator(*this, 0, n_cols, n_nonzero, m.n_nonzero - n_nonzero);
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
typename SpSubview<eT>::row_iterator
SpSubview<eT>::end_row(const uword row_num)
  {
  return row_iterator(*this, row_num + 1, 0);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::end_row(const uword row_num) const
  {
  return const_row_iterator(*this, row_num + 1, 0);
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


//! @}
