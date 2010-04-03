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


//! \addtogroup subview_cube
//! @{


template<typename eT>
inline
subview_cube<eT>::~subview_cube()
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview_cube<eT>::subview_cube
  (
  const Cube<eT>& in_m,
  const u32       in_row1,
  const u32       in_col1,
  const u32       in_slice1,
  const u32       in_row2,
  const u32       in_col2,
  const u32       in_slice2
  )
  : m           (in_m)
  , m_ptr       (0)
  , aux_row1    (in_row1)
  , aux_col1    (in_col1)
  , aux_slice1  (in_slice1)
  , aux_row2    (in_row2)
  , aux_col2    (in_col2)
  , aux_slice2  (in_slice2)
  , n_rows      (1 + in_row2 - in_row1)
  , n_cols      (1 + in_col2 - in_col1)
  , n_elem_slice(n_rows * n_cols)
  , n_slices    (1 + in_slice2 - in_slice1)
  , n_elem      (n_elem_slice * n_slices)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview_cube<eT>::subview_cube
  (
        Cube<eT>& in_m,
  const u32       in_row1,
  const u32       in_col1,
  const u32       in_slice1,
  const u32       in_row2,
  const u32       in_col2,
  const u32       in_slice2
  )
  : m           (in_m)
  , m_ptr       (&in_m)
  , aux_row1    (in_row1)
  , aux_col1    (in_col1)
  , aux_slice1  (in_slice1)
  , aux_row2    (in_row2)
  , aux_col2    (in_col2)
  , aux_slice2  (in_slice2)
  , n_rows      (1 + in_row2 - in_row1)
  , n_cols      (1 + in_col2 - in_col1)
  , n_elem_slice(n_rows * n_cols)
  , n_slices    (1 + in_slice2 - in_slice1)
  , n_elem      (n_elem_slice * n_slices)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
void
subview_cube<eT>::operator+= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 slice = 0; slice < n_slices; ++slice)
    {
    for(u32 col = 0; col < n_cols; ++col)
      {
      eT* coldata = slice_colptr(slice,col);
      
      for(u32 row = 0; row < n_rows; ++row)
        {
        coldata[row] += val;
        }
      
      }
    }

  }



template<typename eT>
inline
void
subview_cube<eT>::operator-= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 slice = 0; slice < n_slices; ++slice)
    {
    for(u32 col = 0; col<n_cols; ++col)
      {
      eT* coldata = slice_colptr(slice,col);
      
      for(u32 row = 0; row<n_rows; ++row)
        {
        coldata[row] -= val;
        }
      
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::operator*= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 slice = 0; slice < n_slices; ++slice)
    {
    for(u32 col = 0; col<n_cols; ++col)
      {
      eT* coldata = slice_colptr(slice,col);
      
      for(u32 row = 0; row<n_rows; ++row)
        {
        coldata[row] *= val;
        }
      
      }
    }  
  }



template<typename eT>
inline
void
subview_cube<eT>::operator/= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 slice = 0; slice < n_slices; ++slice)
    {
    for(u32 col = 0; col<n_cols; ++col)
      {
      eT* coldata = slice_colptr(slice,col);
      
      for(u32 row = 0; row<n_rows; ++row)
        {
        coldata[row] /= val;
        }
      
      }
    }  
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(in.get_ref());
  
  const Cube<eT>&         x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "copy into subcube");
  
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col < t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row < t.n_rows; ++row)
        {
        t_coldata[row] = x_coldata[row];
        }
        
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator+= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(in.get_ref());
  
  const Cube<eT>&         x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "cube addition");
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col < t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row < t.n_rows; ++row)
        {
        t_coldata[row] += x_coldata[row];
        }
      }
    }
  
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator-= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(in.get_ref());
  
  const Cube<eT>&         x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "cube subtraction");
  
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col < t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row < t.n_rows; ++row)
        {
        t_coldata[row] -= x_coldata[row];
        }
      }
    }  
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator%= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(in.get_ref());
  
  const Cube<eT>&         x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "cube schur product");
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col<t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row<t.n_rows; ++row)
        {
        t_coldata[row] *= x_coldata[row];
        }
      }
    }
  
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator/= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(in.get_ref());
  
  const Cube<eT>&         x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "element-wise cube division");
  
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col<t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row<t.n_rows; ++row)
        {
        t_coldata[row] /= x_coldata[row];
        }
      }
    }
  
  }



//! x.subcube(...) = y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::operator= (const subview_cube<eT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Cube<eT>*         tmp_cube         = overlap ? new Cube<eT>(x_in.m) : 0;
  const subview_cube<eT>* tmp_subview_cube = overlap ? new subview_cube<eT>(*tmp_cube, x_in.aux_row1, x_in.aux_col1, x_in.aux_slice1, x_in.aux_row2, x_in.aux_col2, x_in.aux_slice2) : 0;
  const subview_cube<eT>& x                = overlap ? (*tmp_subview_cube) : x_in;
  
  subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "copy into subcube");
  
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col < t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row < t.n_rows; ++row)
        {
        t_coldata[row] = x_coldata[row];
        }
      }
    }
    
  if(overlap)
    {
    delete tmp_subview_cube;
    delete tmp_cube;
    }
  
  }



template<typename eT>
inline
void
subview_cube<eT>::operator+= (const subview_cube<eT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Cube<eT>*         tmp_cube         = overlap ? new Cube<eT>(x_in.m) : 0;
  const subview_cube<eT>* tmp_subview_cube = overlap ? new subview_cube<eT>(*tmp_cube, x_in.aux_row1, x_in.aux_col1, x_in.aux_slice1, x_in.aux_row2, x_in.aux_col2, x_in.aux_slice2) : 0;
  const subview_cube<eT>& x                = overlap ? (*tmp_subview_cube) : x_in;
  
  subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "cube addition");
  
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col < t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row < t.n_rows; ++row)
        {
        t_coldata[row] += x_coldata[row];
        }
      }
    }
    
  if(overlap)
    {
    delete tmp_subview_cube;
    delete tmp_cube;
    }
  
  }



template<typename eT>
inline
void
subview_cube<eT>::operator-= (const subview_cube<eT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Cube<eT>*         tmp_cube         = overlap ? new Cube<eT>(x_in.m) : 0;
  const subview_cube<eT>* tmp_subview_cube = overlap ? new subview_cube<eT>(*tmp_cube, x_in.aux_row1, x_in.aux_col1, x_in.aux_slice1, x_in.aux_row2, x_in.aux_col2, x_in.aux_slice2) : 0;
  const subview_cube<eT>& x                = overlap ? (*tmp_subview_cube) : x_in;
  
  subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "cube subtraction");
  
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col < t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row < t.n_rows; ++row)
        {
        t_coldata[row] -= x_coldata[row];
        }
      }
    }
    
  if(overlap)
    {
    delete tmp_subview_cube;
    delete tmp_cube;
    }
    
  }



template<typename eT>
inline
void
subview_cube<eT>::operator%= (const subview_cube<eT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Cube<eT>*         tmp_cube         = overlap ? new Cube<eT>(x_in.m) : 0;
  const subview_cube<eT>* tmp_subview_cube = overlap ? new subview_cube<eT>(*tmp_cube, x_in.aux_row1, x_in.aux_col1, x_in.aux_slice1, x_in.aux_row2, x_in.aux_col2, x_in.aux_slice2) : 0;
  const subview_cube<eT>& x                = overlap ? (*tmp_subview_cube) : x_in;
  
  subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "element-wise cube multiplication");
  
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col < t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row < t.n_rows; ++row)
        {
        t_coldata[row] *= x_coldata[row];
        }
      }
    }
    
  if(overlap)
    {
    delete tmp_subview_cube;
    delete tmp_cube;
    }
  
  }



template<typename eT>
inline
void
subview_cube<eT>::operator/= (const subview_cube<eT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Cube<eT>*         tmp_cube         = overlap ? new Cube<eT>(x_in.m) : 0;
  const subview_cube<eT>* tmp_subview_cube = overlap ? new subview_cube<eT>(*tmp_cube, x_in.aux_row1, x_in.aux_col1, x_in.aux_slice1, x_in.aux_row2, x_in.aux_col2, x_in.aux_slice2) : 0;
  const subview_cube<eT>& x                = overlap ? (*tmp_subview_cube) : x_in;
  
  subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "element-wise cube division");
  
  
  for(u32 slice = 0; slice < t.n_slices; ++slice)
    {
    for(u32 col = 0; col < t.n_cols; ++col)
      {
            eT* t_coldata = t.slice_colptr(slice,col);
      const eT* x_coldata = x.slice_colptr(slice,col);
      
      for(u32 row = 0; row < t.n_rows; ++row)
        {
        t_coldata[row] /= x_coldata[row];
        }
      }
    }
    
  if(overlap)
    {
    delete tmp_subview_cube;
    delete tmp_cube;
    }
  
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "copy into subcube");
  
  
  for(u32 col = 0; col < t.n_cols; ++col)
    {
          eT* t_coldata = t.slice_colptr(t.aux_slice1, col);
    const eT* x_coldata = x.colptr(col);
    
    for(u32 row = 0; row < t.n_rows; ++row)
      {
      t_coldata[row] = x_coldata[row];
      }
      
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator+= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "cube addition");
  
  for(u32 col = 0; col < t.n_cols; ++col)
    {
          eT* t_coldata = t.slice_colptr(t.aux_slice1, col);
    const eT* x_coldata = x.colptr(col);
    
    for(u32 row = 0; row < t.n_rows; ++row)
      {
      t_coldata[row] += x_coldata[row];
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator-= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "cube subtraction");
  
  for(u32 col = 0; col < t.n_cols; ++col)
    {
          eT* t_coldata = t.slice_colptr(t.aux_slice1, col);
    const eT* x_coldata = x.colptr(col);
    
    for(u32 row = 0; row < t.n_rows; ++row)
      {
      t_coldata[row] -= x_coldata[row];
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator%= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "cube schur product");
  
  for(u32 col = 0; col<t.n_cols; ++col)
    {
          eT* t_coldata = t.slice_colptr(t.aux_slice1, col);
    const eT* x_coldata = x.colptr(col);
    
    for(u32 row = 0; row<t.n_rows; ++row)
      {
      t_coldata[row] *= x_coldata[row];
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator/= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "element-wise cube division");
  
  
  for(u32 col = 0; col<t.n_cols; ++col)
    {
          eT* t_coldata = t.slice_colptr(t.aux_slice1, col);
    const eT* x_coldata = x.colptr(col);
    
    for(u32 row = 0; row<t.n_rows; ++row)
      {
      t_coldata[row] /= x_coldata[row];
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();

  for(u32 slice = 0; slice < n_slices; ++slice)
    {
    for(u32 col = 0; col < n_cols; ++col)
      {
      eT* coldata = slice_colptr(slice,col);
      
      for(u32 row = 0; row < n_rows; ++row)
        {
        coldata[row] = val;
        }
      }
    }
  
  }



template<typename eT>
inline
void
subview_cube<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(0));
  }



template<typename eT>
inline
void
subview_cube<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(1));
  }



template<typename eT>
arma_inline
eT&
subview_cube<eT>::operator[](const u32 i)
  {
  arma_check( (m_ptr == 0), "subview_cube::operator[]: cube is read-only");
  
  const u32 in_slice = i / n_elem_slice;
  const u32 offset   = in_slice * n_elem_slice;
  const u32 j        = i - offset;
  
  const u32 in_col   = j / n_rows;
  const u32 in_row   = j % n_rows;

  const u32 index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview_cube<eT>::operator[](const u32 i) const
  {
  const u32 in_slice = i / n_elem_slice;
  const u32 offset   = in_slice * n_elem_slice;
  const u32 j        = i - offset;
  
  const u32 in_col   = j / n_rows;
  const u32 in_row   = j % n_rows;

  const u32 index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT&
subview_cube<eT>::operator()(const u32 i)
  {
  arma_check( (m_ptr == 0), "subview_cube::operator(): matrix is read-only");
  arma_debug_check( (i >= n_elem), "subview_cube::operator(): index out of bounds");
  
  const u32 in_slice = i / n_elem_slice;
  const u32 offset   = in_slice * n_elem_slice;
  const u32 j        = i - offset;
  
  const u32 in_col   = j / n_rows;
  const u32 in_row   = j % n_rows;

  const u32 index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview_cube<eT>::operator()(const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "subview_cube::operator(): index out of bounds");
  
  const u32 in_slice = i / n_elem_slice;
  const u32 offset   = in_slice * n_elem_slice;
  const u32 j        = i - offset;
  
  const u32 in_col   = j / n_rows;
  const u32 in_row   = j % n_rows;

  const u32 index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT&
subview_cube<eT>::operator()(const u32 in_row, const u32 in_col, const u32 in_slice)
  {
  arma_check( (m_ptr == 0), "subview_cube::operator(): matrix is read-only");
  arma_debug_check( ( (in_row >= n_rows) || (in_col >= n_cols) || (in_slice >= n_slices) ), "subview_cube::operator(): location out of bounds");
  
  const u32 index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview_cube<eT>::operator()(const u32 in_row, const u32 in_col, const u32 in_slice) const
  {
  arma_debug_check( ( (in_row >= n_rows) || (in_col >= n_cols) || (in_slice >= n_slices) ), "subview_cube::operator(): location out of bounds");
  
  const u32 index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT&
subview_cube<eT>::at(const u32 in_row, const u32 in_col, const u32 in_slice)
  {
  arma_check( (m_ptr == 0), "subview_cube::at(): cube is read-only");
  
  const u32 index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview_cube<eT>::at(const u32 in_row, const u32 in_col, const u32 in_slice) const
  {
  const u32 index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT*
subview_cube<eT>::slice_colptr(const u32 in_slice, const u32 in_col)
  {
  arma_check( (m_ptr == 0), "subview_cube::slice_colptr(): cube is read-only");
    
  return & access::rw((*m_ptr).mem[  (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1  ]);
  }



template<typename eT>
arma_inline
const eT*
subview_cube<eT>::slice_colptr(const u32 in_slice, const u32 in_col) const
  {
  return & m.mem[ (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 ];
  }



template<typename eT>
inline
bool
subview_cube<eT>::check_overlap(const subview_cube<eT>& x) const
  {
  const subview_cube<eT>& t = *this;
  
  if(&t.m != &x.m)
    {
    return false;
    }
  else
    {
    const bool row_overlap =
      (
      ( (x.aux_row1 >= t.aux_row1) && (x.aux_row1 <= t.aux_row2) )
      || 
      ( (x.aux_row2 >= t.aux_row1) && (x.aux_row2 <= t.aux_row2) )
      );
    
    const bool col_overlap =
      (
      ( (x.aux_col1 >= t.aux_col1) && (x.aux_col1 <= t.aux_col2) )
      || 
      ( (x.aux_col2 >= t.aux_col1) && (x.aux_col2 <= t.aux_col2) )
      );
    
    const bool slice_overlap =
      (
      ( (x.aux_slice1 >= t.aux_slice1) && (x.aux_slice1 <= t.aux_slice2) )
      || 
      ( (x.aux_slice2 >= t.aux_slice1) && (x.aux_slice2 <= t.aux_slice2) )
      );
    
    return (row_overlap & col_overlap & slice_overlap);
    }
  }



template<typename eT>
inline
bool
subview_cube<eT>::check_overlap(const Mat<eT>& x) const
  {
  const subview_cube<eT>& t = *this;
  
  for(u32 slice = t.aux_slice1; slice <= t.aux_slice2; ++slice)
    {
    const Mat<eT>& y = *(t.m.mat_ptrs[slice]);
  
    if( x.memptr() == y.memptr() )
      {
      return true;
      }
    }
  
  return false;
  }



//! cube X = Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::extract(Cube<eT>& actual_out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  //
  const bool alias = (&actual_out == &in.m);
  
  Cube<eT>* tmp = (alias) ? new Cube<eT> : 0;
  Cube<eT>& out = (alias) ? (*tmp)       : actual_out;
  
  //
  
  const u32 n_rows   = in.n_rows;
  const u32 n_cols   = in.n_cols;
  const u32 n_slices = in.n_slices;
  
  out.set_size(n_rows, n_cols, n_slices);
  
  arma_extra_debug_print(arma_boost::format("out.n_rows = %d   out.n_cols = %d    out.n_slices = %d    in.m.n_rows = %d   in.m.n_cols = %d   in.m.n_slices = %d") % out.n_rows % out.n_cols % out.n_slices % in.m.n_rows % in.m.n_cols % in.m.n_slices);
  
  
  for(u32 slice = 0; slice<n_slices; ++slice)
    {
    for(u32 col = 0; col<n_cols; ++col)
      {
            eT* out_coldata = out.slice_colptr(slice,col);
      const eT*  in_coldata =  in.slice_colptr(slice,col);
      
      for(u32 row = 0; row<n_rows; ++row)
        {
        out_coldata[row] = in_coldata[row];
        }
      
      }
    }
  
  
  if(alias)
    {
    actual_out = out;
    delete tmp;
    }
  
  }



//! mat X = Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::extract(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_slices != 1), "subview_cube::extract(): given subcube doesn't have exactly one slice" );
  
  out.set_size(in.n_rows, in.n_cols);
  
  for(u32 col = 0; col < in.n_cols; ++col)
    {
    const eT* in_coldata  = in.slice_colptr(in.aux_slice1, col);
          eT* out_coldata = out.colptr(col);
    
    for(u32 row = 0; row < in.n_rows; ++row)
      {
      out_coldata[row] = in_coldata[row];
      }
    }
  }



//! cube X += Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::plus_inplace(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "cube addition");
  
  const u32 n_rows   = out.n_rows;
  const u32 n_cols   = out.n_cols;
  const u32 n_slices = out.n_slices;
  
  for(u32 slice = 0; slice<n_slices; ++slice)
    {
    for(u32 col = 0; col<n_cols; ++col)
      {
            eT* out_coldata = out.slice_colptr(slice,col);
      const eT*  in_coldata =  in.slice_colptr(slice,col);
      
      for(u32 row = 0; row<n_rows; ++row)
        {
        out_coldata[row] += in_coldata[row];
        }
      }
    }
  }



//! cube X -= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::minus_inplace(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "cube subtraction");
  
  const u32 n_rows   = out.n_rows;
  const u32 n_cols   = out.n_cols;
  const u32 n_slices = out.n_slices;
  
  for(u32 slice = 0; slice<n_slices; ++slice)
    {
    for(u32 col = 0; col<n_cols; ++col)
      {
            eT* out_coldata = out.slice_colptr(slice,col);
      const eT*  in_coldata =  in.slice_colptr(slice,col);
      
      for(u32 row = 0; row<n_rows; ++row)
        {
        out_coldata[row] -= in_coldata[row];
        }
      }
    }
  }



//! cube X %= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::schur_inplace(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "cube schur product");
  
  const u32 n_rows   = out.n_rows;
  const u32 n_cols   = out.n_cols;
  const u32 n_slices = out.n_slices;
  
  for(u32 slice = 0; slice<n_slices; ++slice)
    {
    for(u32 col = 0; col<n_cols; ++col)
      {
            eT* out_coldata = out.slice_colptr(slice,col);
      const eT*  in_coldata =  in.slice_colptr(slice,col);
      
      for(u32 row = 0; row<n_rows; ++row)
        {
        out_coldata[row] *= in_coldata[row];
        }
      }
    }
  }



//! cube X /= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::div_inplace(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "element-wise cube division");
  
  const u32 n_rows   = out.n_rows;
  const u32 n_cols   = out.n_cols;
  const u32 n_slices = out.n_slices;
  
  for(u32 slice = 0; slice<n_slices; ++slice)
    {
    for(u32 col = 0; col<n_cols; ++col)
      {
            eT* out_coldata = out.slice_colptr(slice,col);
      const eT*  in_coldata =  in.slice_colptr(slice,col);
      
      for(u32 row = 0; row<n_rows; ++row)
        {
        out_coldata[row] /= in_coldata[row];
        }
      }
    }
  }



//! mat X += Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::plus_inplace(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "matrix addition");
  
  for(u32 col = 0; col < in.n_cols; ++col)
    {
    const eT* in_coldata  = in.slice_colptr(in.aux_slice1, col);
          eT* out_coldata = out.colptr(col);
    
    for(u32 row = 0; row < in.n_rows; ++row)
      {
      out_coldata[row] += in_coldata[row];
      }
    }
    
  }



//! mat X -= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::minus_inplace(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "matrix subtraction");
  
  for(u32 col = 0; col < in.n_cols; ++col)
    {
    const eT* in_coldata  = in.slice_colptr(in.aux_slice1, col);
          eT* out_coldata = out.colptr(col);
    
    for(u32 row = 0; row < in.n_rows; ++row)
      {
      out_coldata[row] -= in_coldata[row];
      }
    }
    
  }



//! mat X %= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::schur_inplace(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "matrix schur product");
  
  for(u32 col = 0; col < in.n_cols; ++col)
    {
    const eT* in_coldata  = in.slice_colptr(in.aux_slice1, col);
          eT* out_coldata = out.colptr(col);
    
    for(u32 row = 0; row < in.n_rows; ++row)
      {
      out_coldata[row] *= in_coldata[row];
      }
    }
    
  }



//! mat X /= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::div_inplace(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "matrix element-wise division");
  
  for(u32 col = 0; col < in.n_cols; ++col)
    {
    const eT* in_coldata  = in.slice_colptr(in.aux_slice1, col);
          eT* out_coldata = out.colptr(col);
    
    for(u32 row = 0; row < in.n_rows; ++row)
      {
      out_coldata[row] /= in_coldata[row];
      }
    }
    
  }



//! @}
