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
diagview<eT>::diagview(const Mat<eT>& in_m, const u32 in_row_offset, const u32 in_col_offset, const u32 in_len)
  : m(in_m)
  , m_ptr(0)
  , row_offset(in_row_offset)
  , col_offset(in_col_offset)
  , n_rows(in_len)
  , n_cols( (in_len > 0) ? 1 : 0 )
  , n_elem(in_len)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
diagview<eT>::diagview(Mat<eT>& in_m, const u32 in_row_offset, const u32 in_col_offset, const u32 in_len)
  : m(in_m)
  , m_ptr(&in_m)
  , row_offset(in_row_offset)
  , col_offset(in_col_offset)
  , n_rows(in_len)
  , n_cols( (in_len > 0) ? 1 : 0 )
  , n_elem(in_len)
  {
  arma_extra_debug_sigprint();
  }



//! set a diagonal of our matrix using data from a foreign object
template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator= (const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  arma_check( (m_ptr == 0), "diagview::operator=(): matrix is read only");
  
  const unwrap<T1> tmp(o.get_ref());
  const Mat<eT>& x = tmp.M;
  
  diagview& t = *this;
  
  arma_debug_check( !x.is_vec(), "diagview::operator=(): need a vector");
  arma_debug_check( (t.n_elem != x.n_elem), "diagview::operator=(): diagonal and given vector have incompatible lengths");
  
  Mat<eT>& t_m = *(t.m_ptr);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    t_m.at(i+row_offset, i+col_offset) = x.mem[i];
    }
  
  }



//! set a diagonal of our matrix using a diagonal from a foreign matrix
template<typename eT>
inline
void
diagview<eT>::operator= (const diagview<eT>& x)
  {
  arma_extra_debug_sigprint();
  arma_check( (m_ptr == 0), "diagview::operator=(): matrix is read only");
  
  diagview<eT>& t = *this;
  
  arma_debug_check( (t.n_elem != x.n_elem), "diagview::operator=(): diagonals have incompatible lengths");
  
        Mat<eT>& t_m = *(t.m_ptr);
  const Mat<eT>& x_m = x.m;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    t_m.at(i+t.row_offset, i+t.col_offset) = x_m.at(i+x.row_offset, i+x.col_offset);
    }
    
  }



//! \brief
//! extract a diagonal and store it as a column vector
template<typename eT>
inline
void
diagview<eT>::extract(Mat<eT>& actual_out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& in_m = in.m;
  const bool alias = (&actual_out == &in_m);
  
  Mat<eT>* tmp = (alias) ? new Mat<eT> : 0;
  Mat<eT>& out = (alias) ? (*tmp)      : actual_out;
  
  const u32 in_n_elem     = in.n_elem;
  const u32 in_row_offset = in.row_offset;
  const u32 in_col_offset = in.col_offset;
  
  out.set_size(in_n_elem,1);
  eT* out_mem = out.memptr();
  
  for(u32 i=0; i<in_n_elem; ++i)
    {
    out_mem[i] = in_m.at(i+in_row_offset, i+in_col_offset);
    }
  
  
  if(alias)
    {
    actual_out = out;
    delete tmp;
    }
  
  }



//! X += Y.diagview(...)
template<typename eT>
inline
void
diagview<eT>::plus_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, 1, "matrix addition");
  
  const Mat<eT>& in_m = in.m;
  
  const u32 in_n_elem     = in.n_elem;
  const u32 in_row_offset = in.row_offset;
  const u32 in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  for(u32 i=0; i<in_n_elem; ++i)
    {
    out_mem[i] += in_m.at(i+in_row_offset, i+in_col_offset);
    }
  }



//! X -= Y.diagview(...)
template<typename eT>
inline
void
diagview<eT>::minus_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, 1, "matrix subtraction");
  
  const Mat<eT>& in_m = in.m;
  
  const u32 in_n_elem     = in.n_elem;
  const u32 in_row_offset = in.row_offset;
  const u32 in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  for(u32 i=0; i<in_n_elem; ++i)
    {
    out_mem[i] -= in_m.at(i+in_row_offset, i+in_col_offset);
    }
  }



//! X %= Y.submat(...)
template<typename eT>
inline
void
diagview<eT>::schur_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, 1, "element-wise matrix multiplication");
  
  const Mat<eT>& in_m = in.m;
  
  const u32 in_n_elem     = in.n_elem;
  const u32 in_row_offset = in.row_offset;
  const u32 in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  for(u32 i=0; i<in_n_elem; ++i)
    {
    out_mem[i] *= in_m.at(i+in_row_offset, i+in_col_offset);
    }
  }



//! X /= Y.diagview(...)
template<typename eT>
inline
void
diagview<eT>::div_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, 1, "element-wise matrix division");
  
  const Mat<eT>& in_m = in.m;
  
  const u32 in_n_elem     = in.n_elem;
  const u32 in_row_offset = in.row_offset;
  const u32 in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  for(u32 i=0; i<in_n_elem; ++i)
    {
    out_mem[i] /= in_m.at(i+in_row_offset, i+in_col_offset);
    }
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator[](const u32 i)
  {
  return (*m_ptr).at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator[](const u32 i) const
  {
  return m.at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::at(const u32 row, const u32 col)
  {
  return (*m_ptr).at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::at(const u32 row, const u32 col) const
  {
  return m.at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator()(const u32 i)
  {
  arma_check( (m_ptr == 0), "diagview::operator(): matrix is read only");
  arma_debug_check( (i >= n_elem), "diagview::operator(): out of bounds" );
  
  return (*m_ptr).at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator()(const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "diagview::operator(): out of bounds" );
  
  return m.at(i+row_offset, i+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator()(const u32 row, const u32 col)
  {
  arma_check( (m_ptr == 0), "diagview::operator(): matrix is read only");
  arma_debug_check( ((row >= n_elem) || (col > 0)), "diagview::operator(): out of bounds" );
  
  return (*m_ptr).at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator()(const u32 row, const u32 col) const
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
  arma_check( (m_ptr == 0), "diagview::fill(): matrix is read only");
  
  Mat<eT>& x = (*m_ptr);
  
  for(u32 i=0; i<n_elem; ++i)
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
  arma_check( (m_ptr == 0), "diagview::zeros(): matrix is read only");
  
  Mat<eT>& x = (*m_ptr);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    x.at(i+row_offset, i+col_offset) = eT(0);
    }
  }



template<typename eT>
inline
void
diagview<eT>::ones()
  {
  arma_extra_debug_sigprint();
  arma_check( (m_ptr == 0), "diagview::ones(): matrix is read only");
  
  Mat<eT>& x = (*m_ptr);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    x.at(i+row_offset, i+col_offset) = eT(1);
    }
  }



//! @}
