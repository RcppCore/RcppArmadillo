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


//! \addtogroup diagview
//! @{


template<typename eT>
inline
diagview<eT>::~diagview()
  {
  arma_extra_debug_sigprint();
  }


template<typename eT>
arma_inline
diagview<eT>::diagview(const Mat<eT>& in_m, const uword in_row_offset, const uword in_col_offset, const uword in_len)
  : m(in_m)
  , m_ptr(0)
  , row_offset(in_row_offset)
  , col_offset(in_col_offset)
  , n_rows(in_len)
  , n_elem(in_len)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
diagview<eT>::diagview(Mat<eT>& in_m, const uword in_row_offset, const uword in_col_offset, const uword in_len)
  : m(in_m)
  , m_ptr(&in_m)
  , row_offset(in_row_offset)
  , col_offset(in_col_offset)
  , n_rows(in_len)
  , n_elem(in_len)
  {
  arma_extra_debug_sigprint();
  }



//! set a diagonal of our matrix using a diagonal from a foreign matrix
template<typename eT>
inline
void
diagview<eT>::operator= (const diagview<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>& t = *this;
  
  arma_debug_check( (t.n_elem != x.n_elem), "diagview: diagonals have incompatible lengths");
  
        Mat<eT>& t_m = *(t.m_ptr);
  const Mat<eT>& x_m = x.m;
  
  if(&t_m != &x_m)
    {
    const uword t_n_elem     = t.n_elem;
    const uword t_row_offset = t.row_offset;
    const uword t_col_offset = t.col_offset;
    
    const uword x_row_offset = x.row_offset;
    const uword x_col_offset = x.col_offset;
    
    uword i,j;
    for(i=0, j=1; j < t_n_elem; i+=2, j+=2)
      {
      const eT tmp_i = x_m.at(i + x_row_offset, i + x_col_offset);
      const eT tmp_j = x_m.at(j + x_row_offset, j + x_col_offset);
      
      t_m.at(i + t_row_offset, i + t_col_offset) = tmp_i;
      t_m.at(j + t_row_offset, j + t_col_offset) = tmp_j;
      }
    
    if(i < t_n_elem)
      {
      t_m.at(i + t_row_offset, i + t_col_offset) = x_m.at(i + x_row_offset, i + x_col_offset);
      }
    }
  else
    {
    const Mat<eT> tmp = x;
    
    (*this).operator=(tmp);
    }
  }



template<typename eT>
inline
void
diagview<eT>::operator+=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& t_m = (*m_ptr);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i<t_n_elem; ++i)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) += val;
    }
  }



template<typename eT>
inline
void
diagview<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& t_m = (*m_ptr);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i<t_n_elem; ++i)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) -= val;
    }
  }



template<typename eT>
inline
void
diagview<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& t_m = (*m_ptr);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i<t_n_elem; ++i)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) *= val;
    }
  }



template<typename eT>
inline
void
diagview<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& t_m = (*m_ptr);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i<t_n_elem; ++i)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) /= val;
    }
  }



//! set a diagonal of our matrix using data from a foreign object
template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator= (const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(o.get_ref());
  const Mat<eT>& x = tmp.M;
  
  diagview<eT>& t = *this;
  
  arma_debug_check
    (
    ( (t.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& t_m = *(t.m_ptr);
  
  const uword t_n_elem     = t.n_elem;
  const uword t_row_offset = t.row_offset;
  const uword t_col_offset = t.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < t_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    t_m.at( i + t_row_offset,  i + t_col_offset) = tmp_i;
    t_m.at( j + t_row_offset,  j + t_col_offset) = tmp_j;
    }
  
  if(i < t_n_elem)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) = x_mem[i];
    }
  }



template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator+=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(o.get_ref());
  const Mat<eT>& x = tmp.M;
  
  diagview<eT>& t = *this;
  
  arma_debug_check
    (
    ( (t.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& t_m = *(t.m_ptr);
  
  const uword t_n_elem     = t.n_elem;
  const uword t_row_offset = t.row_offset;
  const uword t_col_offset = t.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < t_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    t_m.at( i + t_row_offset,  i + t_col_offset) += tmp_i;
    t_m.at( j + t_row_offset,  j + t_col_offset) += tmp_j;
    }
  
  if(i < t_n_elem)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) += x_mem[i];
    }
  }



template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator-=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(o.get_ref());
  const Mat<eT>& x = tmp.M;
  
  diagview<eT>& t = *this;
  
  arma_debug_check
    (
    ( (t.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& t_m = *(t.m_ptr);
  
  const uword t_n_elem     = t.n_elem;
  const uword t_row_offset = t.row_offset;
  const uword t_col_offset = t.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < t_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    t_m.at( i + t_row_offset,  i + t_col_offset) -= tmp_i;
    t_m.at( j + t_row_offset,  j + t_col_offset) -= tmp_j;
    }
  
  if(i < t_n_elem)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) -= x_mem[i];
    }
  }



template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator%=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(o.get_ref());
  const Mat<eT>& x = tmp.M;
  
  diagview<eT>& t = *this;
  
  arma_debug_check
    (
    ( (t.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& t_m = *(t.m_ptr);
  
  const uword t_n_elem     = t.n_elem;
  const uword t_row_offset = t.row_offset;
  const uword t_col_offset = t.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < t_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    t_m.at( i + t_row_offset,  i + t_col_offset) *= tmp_i;
    t_m.at( j + t_row_offset,  j + t_col_offset) *= tmp_j;
    }
  
  if(i < t_n_elem)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) *= x_mem[i];
    }
  }



template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator/=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(o.get_ref());
  const Mat<eT>& x = tmp.M;
  
  diagview<eT>& t = *this;
  
  arma_debug_check
    (
    ( (t.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& t_m = *(t.m_ptr);
  
  const uword t_n_elem     = t.n_elem;
  const uword t_row_offset = t.row_offset;
  const uword t_col_offset = t.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < t_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    t_m.at( i + t_row_offset,  i + t_col_offset) /= tmp_i;
    t_m.at( j + t_row_offset,  j + t_col_offset) /= tmp_j;
    }
  
  if(i < t_n_elem)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) /= x_mem[i];
    }
  }



//! extract a diagonal and store it as a column vector
template<typename eT>
inline
void
diagview<eT>::extract(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: we're assuming that the matrix has already been set to the correct size and there is no aliasing;
  // size setting and alias checking is done by either the Mat contructor or operator=()
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] = tmp_i;
    out_mem[j] = tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] = in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



//! X += Y.diag()
template<typename eT>
inline
void
diagview<eT>::plus_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, in.n_cols, "addition");
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] += tmp_i;
    out_mem[j] += tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] += in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



//! X -= Y.diag()
template<typename eT>
inline
void
diagview<eT>::minus_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, in.n_cols, "subtraction");
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] -= tmp_i;
    out_mem[j] -= tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] -= in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



//! X %= Y.diag()
template<typename eT>
inline
void
diagview<eT>::schur_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, in.n_cols, "element-wise multiplication");
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] *= tmp_i;
    out_mem[j] *= tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] *= in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



//! X /= Y.diag()
template<typename eT>
inline
void
diagview<eT>::div_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, in.n_cols, "element-wise division");
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] /= tmp_i;
    out_mem[j] /= tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] /= in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator[](const uword i)
  {
  return (*m_ptr).at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator[](const uword i) const
  {
  return m.at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::at(const uword i)
  {
  return (*m_ptr).at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::at(const uword i) const
  {
  return m.at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator()(const uword i)
  {
  arma_debug_check( (i >= n_elem), "diagview::operator(): out of bounds" );
  
  return (*m_ptr).at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator()(const uword i) const
  {
  arma_debug_check( (i >= n_elem), "diagview::operator(): out of bounds" );
  
  return m.at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::at(const uword row, const uword col)
  {
  return (*m_ptr).at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::at(const uword row, const uword col) const
  {
  return m.at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator()(const uword row, const uword col)
  {
  arma_debug_check( ((row >= n_elem) || (col > 0)), "diagview::operator(): out of bounds" );
  
  return (*m_ptr).at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator()(const uword row, const uword col) const
  {
  arma_debug_check( ((row >= n_elem) || (col > 0)), "diagview::operator(): out of bounds" );
  
  return m.at(row+row_offset, row+col_offset);
  }



template<typename eT>
inline
void
diagview<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& x = (*m_ptr);
  
  for(uword i=0; i<n_elem; ++i)
    {
    x.at(i+row_offset, i+col_offset) = val;
    }
  }



template<typename eT>
inline
void
diagview<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(0));
  }



template<typename eT>
inline
void
diagview<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(1));
  }



//! @}
