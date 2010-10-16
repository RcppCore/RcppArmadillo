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


//! \addtogroup subview
//! @{


template<typename eT>
inline
subview<eT>::~subview()
  {
  arma_extra_debug_sigprint();
  }


template<typename eT>
arma_inline
subview<eT>::subview(const Mat<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_row2,  const u32 in_col2)
  : m(in_m)
  , m_ptr(0)
  , aux_row1(in_row1)
  , aux_col1(in_col1)
  , aux_row2(in_row2)
  , aux_col2(in_col2)
  , n_rows(1 + in_row2 - in_row1)
  , n_cols(1 + in_col2 - in_col1)
  , n_elem(n_rows*n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview<eT>::subview(Mat<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_row2,  const u32 in_col2)
  : m(in_m)
  , m_ptr(&in_m)
  , aux_row1(in_row1)
  , aux_col1(in_col1)
  , aux_row2(in_row2)
  , aux_col2(in_col2)
  , n_rows(1 + in_row2 - in_row1)
  , n_cols(1 + in_col2 - in_col1)
  , n_elem(n_rows*n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
void
subview<eT>::operator+= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  const u32 local_n_cols = n_cols;
  const u32 local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      at(0, col) += val;
      }
    }
  else
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      arrayops::inplace_plus( colptr(col), val, local_n_rows );
      }
    }
  }



template<typename eT>
inline
void
subview<eT>::operator-= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  const u32 local_n_cols = n_cols;
  const u32 local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      at(0, col) -= val;
      }
    }
  else
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      arrayops::inplace_minus( colptr(col), val, local_n_rows );
      }
    }
  }



template<typename eT>
inline
void
subview<eT>::operator*= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  const u32 local_n_cols = n_cols;
  const u32 local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      at(0, col) *= val;
      }
    }
  else
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      arrayops::inplace_mul( colptr(col), val, local_n_rows );
      }
    }
  }



template<typename eT>
inline
void
subview<eT>::operator/= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  const u32 local_n_cols = n_cols;
  const u32 local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      at(0, col) /= val;
      }
    }
  else
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      arrayops::inplace_div( colptr(col), val, local_n_rows );
      }
    }
  }



template<typename eT>
template<typename T1>
arma_inline
void
subview<eT>::operator= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1>   tmp(in.get_ref());
  const Mat<eT>& x = tmp.M;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "insert into submatrix");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    const eT* x_mem = x.memptr();
    
    for(u32 col=0; col<t_n_cols; ++col)
      {
      at(0,col) = x_mem[col];
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      syslib::copy_elem( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview<eT>::operator+= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1>   tmp(in.get_ref());
  const Mat<eT>& x = tmp.M;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "matrix addition");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    const eT* x_mem = x.memptr();
    
    for(u32 col=0; col<t_n_cols; ++col)
      {
      at(0,col) += x_mem[col];
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      arrayops::inplace_plus( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview<eT>::operator-= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1>   tmp(in.get_ref());
  const Mat<eT>& x = tmp.M;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "matrix subtraction");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    const eT* x_mem = x.memptr();
    
    for(u32 col=0; col<t_n_cols; ++col)
      {
      at(0,col) -= x_mem[col];
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      arrayops::inplace_minus( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview<eT>::operator%= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1>   tmp(in.get_ref());
  const Mat<eT>& x = tmp.M;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "matrix schur product");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    const eT* x_mem = x.memptr();
    
    for(u32 col=0; col<t_n_cols; ++col)
      {
      at(0,col) *= x_mem[col];
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      arrayops::inplace_mul( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview<eT>::operator/= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1>   tmp(in.get_ref());
  const Mat<eT>& x = tmp.M;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "element-wise matrix division");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    const eT* x_mem = x.memptr();
    
    for(u32 col=0; col<t_n_cols; ++col)
      {
      at(0,col) /= x_mem[col];
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      arrayops::inplace_div( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
  }



//! x.submat(...) = y.submat(...)
template<typename eT>
inline
void
subview<eT>::operator= (const subview<eT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Mat<eT>*     tmp_mat     = overlap ? new Mat<eT>(x_in.m) : 0;
  const subview<eT>* tmp_subview = overlap ? new subview<eT>(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.aux_row2, x_in.aux_col2) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "insert into submatrix");
  
  const u32 t_n_cols = t.n_cols;
  const u32 t_n_rows = t.n_rows;
  
  if(t_n_rows == 1)
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      t.at(0,col) = x.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      syslib::copy_elem( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
  
  if(overlap)
    {
    delete tmp_subview;
    delete tmp_mat;
    }
  }



template<typename eT>
inline
void
subview<eT>::operator+= (const subview<eT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Mat<eT>*     tmp_mat     = overlap ? new Mat<eT>(x_in.m) : 0;
  const subview<eT>* tmp_subview = overlap ? new subview(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.aux_row2, x_in.aux_col2) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "matrix addition");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      t.at(0,col) += x.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      arrayops::inplace_plus( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
  
  if(overlap)
    {
    delete tmp_subview;
    delete tmp_mat;
    }
  }



template<typename eT>
inline
void
subview<eT>::operator-= (const subview<eT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Mat<eT>*     tmp_mat     = overlap ? new Mat<eT>(x_in.m) : 0;
  const subview<eT>* tmp_subview = overlap ? new subview(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.aux_row2, x_in.aux_col2) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "matrix subtraction");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      t.at(0,col) -= x.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      arrayops::inplace_minus( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
    
  if(overlap)
    {
    delete tmp_subview;
    delete tmp_mat;
    }
  
  }



template<typename eT>
inline
void
subview<eT>::operator%= (const subview& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Mat<eT>*     tmp_mat     = overlap ? new Mat<eT>(x_in.m) : 0;
  const subview<eT>* tmp_subview = overlap ? new subview(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.aux_row2, x_in.aux_col2) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "matrix schur product");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      t.at(0,col) *= x.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      arrayops::inplace_mul( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
  
  if(overlap)
    {
    delete tmp_subview;
    delete tmp_mat;
    }
  
  }



template<typename eT>
inline
void
subview<eT>::operator/= (const subview& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
  
        Mat<eT>*     tmp_mat     = overlap ? new Mat<eT>(x_in.m) : 0;
  const subview<eT>* tmp_subview = overlap ? new subview(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.aux_row2, x_in.aux_col2) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "element-wise matrix division");
  
  const u32 t_n_rows = t.n_rows;
  const u32 t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      t.at(0,col) /= x.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<t_n_cols; ++col)
      {
      arrayops::inplace_div( t.colptr(col), x.colptr(col), t_n_rows );
      }
    }
    
  if(overlap)
    {
    delete tmp_subview;
    delete tmp_mat;
    }
  
  }



template<typename eT>
inline
void
subview<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  const u32 local_n_cols = n_cols;
  const u32 local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      at(0,col) = val;
      }
    }
  else
    {
    for(u32 col=0; col<local_n_cols; ++col)
      {
      arrayops::inplace_set( colptr(col), val, local_n_rows );
      }
    }
  }



template<typename eT>
inline
void
subview<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(0));
  }



template<typename eT>
arma_inline
eT&
subview<eT>::operator[](const u32 i)
  {
  const u32 in_col = i / n_rows;
  const u32 in_row = i % n_rows;
    
  const u32 index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview<eT>::operator[](const u32 i) const
  {
  const u32 in_col = i / n_rows;
  const u32 in_row = i % n_rows;
  
  const u32 index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT&
subview<eT>::operator()(const u32 i)
  {
  arma_debug_check( (i >= n_elem), "subview::operator(): index out of bounds");
    
  const u32 in_col = i / n_rows;
  const u32 in_row = i % n_rows;
  
  const u32 index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview<eT>::operator()(const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "subview::operator(): index out of bounds");
  
  const u32 in_col = i / n_rows;
  const u32 in_row = i % n_rows;
  
  const u32 index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT&
subview<eT>::operator()(const u32 in_row, const u32 in_col)
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "subview::operator(): index out of bounds");
  
  const u32 index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview<eT>::operator()(const u32 in_row, const u32 in_col) const
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "subview::operator(): index out of bounds");
  
  const u32 index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT&
subview<eT>::at(const u32 in_row, const u32 in_col)
  {
  const u32 index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview<eT>::at(const u32 in_row, const u32 in_col) const
  {
  const u32 index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT*
subview<eT>::colptr(const u32 in_col)
  {
  return & access::rw((*m_ptr).mem[ (in_col + aux_col1)*m.n_rows + aux_row1 ]);
  }



template<typename eT>
arma_inline
const eT*
subview<eT>::colptr(const u32 in_col) const
  {
  return & m.mem[ (in_col + aux_col1)*m.n_rows + aux_row1 ];
  }



template<typename eT>
inline
bool
subview<eT>::check_overlap(const subview<eT>& x) const
  {
  const subview<eT>& t = *this;
  
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
    
    const bool overlap = ( (row_overlap == true) && (col_overlap == true) );
    
    return overlap;
    }
  }



template<typename eT>
inline
bool
subview<eT>::is_vec() const
  {
  return ( (n_rows == 1) || (n_cols == 1) );
  }



//! X = Y.submat(...)
template<typename eT>
inline
void
subview<eT>::extract(Mat<eT>& actual_out, const subview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  //
  const bool alias = (&actual_out == &in.m);
  
  Mat<eT>* tmp = (alias) ? new Mat<eT> : 0;
  Mat<eT>& out = (alias) ? (*tmp)      : actual_out;
  
  //
  
  const u32 n_rows = in.n_rows;  // number of rows in the subview
  const u32 n_cols = in.n_cols;  // number of columns in the subview
  
  out.set_size(n_rows, n_cols);
  
  arma_extra_debug_print(arma_boost::format("out.n_rows = %d   out.n_cols = %d    in.m.n_rows = %d  in.m.n_cols = %d") % out.n_rows % out.n_cols % in.m.n_rows % in.m.n_cols );
  

  if(in.is_vec() == true)
    {
    if(n_cols == 1)   // a column vector
      {
      arma_extra_debug_print("subview::extract(): copying col (going across rows)");
      
      // in.colptr(0) the first column of the subview, taking into account any row offset
      syslib::copy_elem( out.memptr(), in.colptr(0), n_rows );
      }
    else   // a row vector
      {
      arma_extra_debug_print("subview::extract(): copying row (going across columns)");
      
      const Mat<eT>& X = in.m;
      
            eT* out_mem   = out.memptr();
      const u32 row       = in.aux_row1;
      const u32 start_col = in.aux_col1;
      
      for(u32 i=0; i<n_cols; ++i)
        {
        out_mem[i] = X.at(row, i+start_col);
        }
      }
    }
  else   // general submatrix
    {
    arma_extra_debug_print("subview::extract(): general submatrix");
    
    for(u32 col = 0; col<n_cols; ++col)   
      {
      syslib::copy_elem( out.colptr(col), in.colptr(col), n_rows );
      }
    }
  
  
  if(alias)
    {
    actual_out = out;
    delete tmp;
    }
  
  }



//! X += Y.submat(...)
template<typename eT>
inline
void
subview<eT>::plus_inplace(Mat<eT>& out, const subview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "matrix addition");
  
  const u32 n_rows = out.n_rows;
  const u32 n_cols = out.n_cols;
  
  if(n_rows == 1)
    {
    eT* out_mem = out.memptr();
    
    for(u32 col=0; col<n_cols; ++col)
      {
      out_mem[col] += in.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<n_cols; ++col)
      {
      arrayops::inplace_plus(out.colptr(col), in.colptr(col), n_rows);
      }
    }
  }



//! X -= Y.submat(...)
template<typename eT>
inline
void
subview<eT>::minus_inplace(Mat<eT>& out, const subview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "matrix subtraction");
  
  const u32 n_rows = out.n_rows;
  const u32 n_cols = out.n_cols;
  
  if(n_rows == 1)
    {
    eT* out_mem = out.memptr();
    
    for(u32 col=0; col<n_cols; ++col)
      {
      out_mem[col] -= in.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<n_cols; ++col)
      {
      arrayops::inplace_minus(out.colptr(col), in.colptr(col), n_rows);
      }
    }
  }



//! X %= Y.submat(...)
template<typename eT>
inline
void
subview<eT>::schur_inplace(Mat<eT>& out, const subview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "matrix schur product");
  
  const u32 n_rows = out.n_rows;
  const u32 n_cols = out.n_cols;
  
  if(n_rows == 1)
    {
    eT* out_mem = out.memptr();
    
    for(u32 col=0; col<n_cols; ++col)
      {
      out_mem[col] *= in.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<n_cols; ++col)
      {
      arrayops::inplace_mul(out.colptr(col), in.colptr(col), n_rows);
      }
    }
  }



//! X /= Y.submat(...)
template<typename eT>
inline
void
subview<eT>::div_inplace(Mat<eT>& out, const subview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "element-wise matrix division");
  
  const u32 n_rows = out.n_rows;
  const u32 n_cols = out.n_cols;
  
  if(n_rows == 1)
    {
    eT* out_mem = out.memptr();
    
    for(u32 col=0; col<n_cols; ++col)
      {
      out_mem[col] /= in.at(0,col);
      }
    }
  else
    {
    for(u32 col=0; col<n_cols; ++col)
      {
      arrayops::inplace_div(out.colptr(col), in.colptr(col), n_rows);
      }
    }
  }



//
//
//



template<typename eT>
arma_inline
subview_col<eT>::subview_col(const Mat<eT>& in_m, const u32 in_col)
  : subview<eT>(in_m, 0, in_col, in_m.n_rows-1, in_col)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview_col<eT>::subview_col(Mat<eT>& in_m, const u32 in_col)
  : subview<eT>(in_m, 0, in_col, in_m.n_rows-1, in_col)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview_col<eT>::subview_col(const Mat<eT>& in_m, const u32 in_col, const u32 in_row1, const u32 in_row2)
  : subview<eT>(in_m, in_row1, in_col, in_row2, in_col)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview_col<eT>::subview_col(Mat<eT>& in_m, const u32 in_col, const u32 in_row1, const u32 in_row2)
  : subview<eT>(in_m, in_row1, in_col, in_row2, in_col)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
void
subview_col<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::operator=(X);
  arma_debug_check( (subview<eT>::n_cols > 1), "subview_col(): incompatible dimensions" );
  }



template<typename eT>
inline
void
subview_col<eT>::operator=(const subview_col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::operator=(X); // interprets 'subview_col' as 'subview'
  arma_debug_check( (subview<eT>::n_cols > 1), "subview_col(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1>
inline
void
subview_col<eT>::operator=(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::operator=(X);
  arma_debug_check( (subview<eT>::n_cols > 1), "subview_col(): incompatible dimensions" );
  }




//
//
//



template<typename eT>
arma_inline
subview_row<eT>::subview_row(const Mat<eT>& in_m, const u32 in_row)
  : subview<eT>(in_m, in_row, 0, in_row, in_m.n_cols-1)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview_row<eT>::subview_row(Mat<eT>& in_m, const u32 in_row)
  : subview<eT>(in_m, in_row, 0, in_row, in_m.n_cols-1)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview_row<eT>::subview_row(const Mat<eT>& in_m, const u32 in_row, const u32 in_col1, const u32 in_col2)
  : subview<eT>(in_m, in_row, in_col1, in_row, in_col2)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
arma_inline
subview_row<eT>::subview_row(Mat<eT>& in_m, const u32 in_row, const u32 in_col1, const u32 in_col2)
  : subview<eT>(in_m, in_row, in_col1, in_row, in_col2)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
void
subview_row<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::operator=(X);
  arma_debug_check( (subview<eT>::n_rows > 1), "subview_row(): incompatible dimensions" );
  }



template<typename eT>
inline
void
subview_row<eT>::operator=(const subview_row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::operator=(X); // interprets 'subview_row' as 'subview'
  arma_debug_check( (subview<eT>::n_rows > 1), "subview_row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1>
inline
void
subview_row<eT>::operator=(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::operator=(X);
  arma_debug_check( (subview<eT>::n_rows > 1), "subview_row(): incompatible dimensions" );
  }



//! @}
