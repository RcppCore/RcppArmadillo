// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// Copyright (C)      2011 James Sanders
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
inline
subview<eT>::subview(const Mat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols)
  : m(in_m)
  , m_ptr(0)
  , aux_row1(in_row1)
  , aux_col1(in_col1)
  , n_rows(in_n_rows)
  , n_cols(in_n_cols)
  , n_elem(in_n_rows*in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview<eT>::subview(Mat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols)
  : m(in_m)
  , m_ptr(&in_m)
  , aux_row1(in_row1)
  , aux_col1(in_col1)
  , n_rows(in_n_rows)
  , n_cols(in_n_cols)
  , n_elem(in_n_rows*in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
void
subview<eT>::operator+= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_cols = n_cols;
  const uword local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    Mat<eT>& X = (*m_ptr);
    
    const uword row           = aux_row1;
    const uword start_col     = aux_col1;
    const uword end_col_plus1 = start_col + local_n_cols;
    
    uword i,j;
    
    for(i=start_col, j=start_col+1; j < end_col_plus1; i+=2, j+=2)
      {
      X.at(row, i) += val;
      X.at(row, j) += val;
      }
    
    if(i < end_col_plus1)
      {
      X.at(row, i) += val;
      }
    }
  else
    {
    for(uword col=0; col<local_n_cols; ++col)
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
  
  const uword local_n_cols = n_cols;
  const uword local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    Mat<eT>& X = (*m_ptr);
    
    const uword row           = aux_row1;
    const uword start_col     = aux_col1;
    const uword end_col_plus1 = start_col + local_n_cols;
    
    uword i,j;
    
    for(i=start_col, j=start_col+1; j < end_col_plus1; i+=2, j+=2)
      {
      X.at(row, i) -= val;
      X.at(row, j) -= val;
      }
    
    if(i < end_col_plus1)
      {
      X.at(row, i) -= val;
      }
    }
  else
    {
    for(uword col=0; col<local_n_cols; ++col)
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
  
  const uword local_n_cols = n_cols;
  const uword local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    Mat<eT>& X = (*m_ptr);
    
    const uword row           = aux_row1;
    const uword start_col     = aux_col1;
    const uword end_col_plus1 = start_col + local_n_cols;
    
    uword i,j;
    
    for(i=start_col, j=start_col+1; j < end_col_plus1; i+=2, j+=2)
      {
      X.at(row, i) *= val;
      X.at(row, j) *= val;
      }
    
    if(i < end_col_plus1)
      {
      X.at(row, i) *= val;
      }
    }
  else
    {
    for(uword col=0; col<local_n_cols; ++col)
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
  
  const uword local_n_cols = n_cols;
  const uword local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    Mat<eT>& X = (*m_ptr);
    
    const uword row           = aux_row1;
    const uword start_col     = aux_col1;
    const uword end_col_plus1 = start_col + local_n_cols;
    
    uword i,j;
    
    for(i=start_col, j=start_col+1; j < end_col_plus1; i+=2, j+=2)
      {
      X.at(row, i) /= val;
      X.at(row, j) /= val;
      }
    
    if(i < end_col_plus1)
      {
      X.at(row, i) /= val;
      }
    }
  else
    {
    for(uword col=0; col<local_n_cols; ++col)
      {
      arrayops::inplace_div( colptr(col), val, local_n_rows );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview<eT>::operator= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> P(in.get_ref());
  
  subview<eT>& t = *this;
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
    
  arma_debug_assert_same_size(t, P, "insert into submatrix");
  
  const bool alias = P.is_alias(t.m);
  
  arma_extra_debug_warn(alias, "aliasing detected");
  
  if( (alias == true) || (is_Mat<typename Proxy<T1>::stored_type>::value == true) )
    {
    const unwrap_check<typename Proxy<T1>::stored_type> tmp(P.Q, t.m);
    const Mat<eT>& x = tmp.M;
    
    if(t_n_rows == 1)
      {
      const eT* x_mem = x.memptr();
      
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        A.at(row, start_col+i) = x_mem[i];
        A.at(row, start_col+j) = x_mem[j];
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) = x_mem[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        arrayops::copy( t.colptr(col), x.colptr(col), t_n_rows );
        }
      }
    }
  else
    {
    if(t_n_rows == 1)
      {
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp1 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        const eT tmp2 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,j) : P[j];
        
        A.at(row, start_col+i) = tmp1;
        A.at(row, start_col+j) = tmp2;
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) = (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        eT* t_col_data = t.colptr(col);
        
        uword i,j;
        for(i=0, j=1; j<t_n_rows; i+=2, j+=2)
          {
          const eT tmp1 = P.at(i,col);
          const eT tmp2 = P.at(j,col);
          
          t_col_data[i] = tmp1;
          t_col_data[j] = tmp2;
          }
        
        if(i < t_n_rows)
          {
          t_col_data[i] = P.at(i,col);
          }
        }
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
  
  const Proxy<T1> P(in.get_ref());
  
  subview<eT>& t = *this;
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
  
  arma_debug_assert_same_size(t, P, "addition");
  
  const bool alias = P.is_alias(t.m);
  
  arma_extra_debug_warn(alias, "aliasing detected");
  
  if( (alias == true) || (is_Mat<typename Proxy<T1>::stored_type>::value == true) )
    {
    const unwrap_check<typename Proxy<T1>::stored_type> tmp(P.Q, t.m);
    const Mat<eT>& x = tmp.M;
    
    if(t_n_rows == 1)
      {
      const eT* x_mem = x.memptr();
      
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        A.at(row, start_col+i) += x_mem[i];
        A.at(row, start_col+j) += x_mem[j];
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) += x_mem[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        arrayops::inplace_plus( t.colptr(col), x.colptr(col), t_n_rows );
        }
      }
    }
  else
    {
    if(t_n_rows == 1)
      {
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp1 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        const eT tmp2 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,j) : P[j];
        
        A.at(row, start_col+i) += tmp1;
        A.at(row, start_col+j) += tmp2;
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) += (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        eT* t_col_data = t.colptr(col);
        
        uword i,j;
        for(i=0, j=1; j<t_n_rows; i+=2, j+=2)
          {
          const eT val1 = P.at(i,col);
          const eT val2 = P.at(j,col);
          
          t_col_data[i] += val1;
          t_col_data[j] += val2;
          }
        
        if(i < t_n_rows)
          {
          t_col_data[i] += P.at(i,col);
          }
        }
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
  
  const Proxy<T1> P(in.get_ref());
  
  subview<eT>& t = *this;
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
  
  arma_debug_assert_same_size(t, P, "subtraction");
  
  const bool alias = P.is_alias(t.m);
  
  if( (alias == true) || (is_Mat<typename Proxy<T1>::stored_type>::value == true) )
    {
    const unwrap_check<typename Proxy<T1>::stored_type> tmp(P.Q, t.m);
    const Mat<eT>& x = tmp.M;
    
    if(t_n_rows == 1)
      {
      const eT* x_mem = x.memptr();
      
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        A.at(row, start_col+i) -= x_mem[i];
        A.at(row, start_col+j) -= x_mem[j];
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) -= x_mem[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        arrayops::inplace_minus( t.colptr(col), x.colptr(col), t_n_rows );
        }
      }
    }
  else
    {
    if(t_n_rows == 1)
      {
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp1 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        const eT tmp2 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,j) : P[j];
        
        A.at(row, start_col+i) -= tmp1;
        A.at(row, start_col+j) -= tmp2;
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) -= (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        eT* t_col_data = t.colptr(col);
        
        uword i,j;
        for(i=0, j=1; j<t_n_rows; i+=2, j+=2)
          {
          const eT val1 = P.at(i,col);
          const eT val2 = P.at(j,col);
          
          t_col_data[i] -= val1;
          t_col_data[j] -= val2;
          }
        
        if(i < t_n_rows)
          {
          t_col_data[i] -= P.at(i,col);
          }
        }
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
  
  const Proxy<T1> P(in.get_ref());
  
  subview<eT>& t = *this;
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
  
  arma_debug_assert_same_size(t, P, "element-wise multiplication");
  
  const bool alias = P.is_alias(t.m);
  
  arma_extra_debug_warn(alias, "aliasing detected");
  
  if( (alias == true) || (is_Mat<typename Proxy<T1>::stored_type>::value == true) )
    {
    const unwrap_check<typename Proxy<T1>::stored_type> tmp(P.Q, t.m);
    const Mat<eT>& x = tmp.M;
    
    if(t_n_rows == 1)
      {
      const eT* x_mem = x.memptr();
      
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        A.at(row, start_col+i) *= x_mem[i];
        A.at(row, start_col+j) *= x_mem[j];
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) *= x_mem[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        arrayops::inplace_mul( t.colptr(col), x.colptr(col), t_n_rows );
        }
      }
    }
  else
    {
    if(t_n_rows == 1)
      {
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp1 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        const eT tmp2 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,j) : P[j];
        
        A.at(row, start_col+i) *= tmp1;
        A.at(row, start_col+j) *= tmp2;
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) *= (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        eT* t_col_data = t.colptr(col);
        
        uword i,j;
        for(i=0, j=1; j<t_n_rows; i+=2, j+=2)
          {
          const eT val1 = P.at(i,col);
          const eT val2 = P.at(j,col);
          
          t_col_data[i] *= val1;
          t_col_data[j] *= val2;
          }
        
        if(i < t_n_rows)
          {
          t_col_data[i] *= P.at(i,col);
          }
        }
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
  
  const Proxy<T1> P(in.get_ref());
  
  subview<eT>& t = *this;
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
  
  arma_debug_assert_same_size(t, P, "element-wise division");
  
  const bool alias = P.is_alias(t.m);
  
  arma_extra_debug_warn(alias, "aliasing detected");
  
  if( (alias == true) || (is_Mat<typename Proxy<T1>::stored_type>::value == true) )
    {
    const unwrap_check<typename Proxy<T1>::stored_type> tmp(P.Q, t.m);
    const Mat<eT>& x = tmp.M;
    
    if(t_n_rows == 1)
      {
      const eT* x_mem = x.memptr();
      
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        A.at(row, start_col+i) /= x_mem[i];
        A.at(row, start_col+j) /= x_mem[j];
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) /= x_mem[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        arrayops::inplace_div( t.colptr(col), x.colptr(col), t_n_rows );
        }
      }
    }
  else
    {
    if(t_n_rows == 1)
      {
      Mat<eT>& A = (*m_ptr);
      
      const uword row       = aux_row1;
      const uword start_col = aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp1 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        const eT tmp2 = (Proxy<T1>::prefer_at_accessor) ? P.at(0,j) : P[j];
        
        A.at(row, start_col+i) /= tmp1;
        A.at(row, start_col+j) /= tmp2;
        }
      
      if(i < t_n_cols)
        {
        A.at(row, start_col+i) /= (Proxy<T1>::prefer_at_accessor) ? P.at(0,i) : P[i];
        }
      }
    else
      {
      for(uword col=0; col<t_n_cols; ++col)
        {
        eT* t_col_data = t.colptr(col);
        
        uword i,j;
        for(i=0, j=1; j<t_n_rows; i+=2, j+=2)
          {
          const eT val1 = P.at(i,col);
          const eT val2 = P.at(j,col);
          
          t_col_data[i] /= val1;
          t_col_data[j] /= val2;
          }
        
        if(i < t_n_rows)
          {
          t_col_data[i] /= P.at(i,col);
          }
        }
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
  const subview<eT>* tmp_subview = overlap ? new subview<eT>(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.n_rows, x_in.n_cols) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "insert into submatrix");
  
  const uword t_n_cols = t.n_cols;
  const uword t_n_rows = t.n_rows;
  
  if(t_n_rows == 1)
    {
          Mat<eT>& A = *(t.m_ptr);
    const Mat<eT>& B = x.m;
    
    const uword row_A = t.aux_row1;
    const uword row_B = x.aux_row1;
    
    const uword start_col_A = t.aux_col1;
    const uword start_col_B = x.aux_col1;
    
    uword i,j;
    
    for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
      {
      const eT tmp1 = B.at(row_B, start_col_B + i);
      const eT tmp2 = B.at(row_B, start_col_B + j);
      
      A.at(row_A, start_col_A + i) = tmp1;
      A.at(row_A, start_col_A + j) = tmp2;
      }
    
    if(i < t_n_cols)
      {
      A.at(row_A, start_col_A + i) = B.at(row_B, start_col_B + i);
      }
    }
  else
    {
    for(uword col=0; col<t_n_cols; ++col)
      {
      arrayops::copy( t.colptr(col), x.colptr(col), t_n_rows );
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
  const subview<eT>* tmp_subview = overlap ? new subview(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.n_rows, x_in.n_cols) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "addition");
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
          Mat<eT>& A = *(t.m_ptr);
    const Mat<eT>& B = x.m;
    
    const uword row_A = t.aux_row1;
    const uword row_B = x.aux_row1;
    
    const uword start_col_A = t.aux_col1;
    const uword start_col_B = x.aux_col1;
    
    uword i,j;
    
    for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
      {
      const eT tmp1 = B.at(row_B, start_col_B + i);
      const eT tmp2 = B.at(row_B, start_col_B + j);
      
      A.at(row_A, start_col_A + i) += tmp1;
      A.at(row_A, start_col_A + j) += tmp2;
      }
    
    if(i < t_n_cols)
      {
      A.at(row_A, start_col_A + i) += B.at(row_B, start_col_B + i);
      }
    }
  else
    {
    for(uword col=0; col<t_n_cols; ++col)
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
  const subview<eT>* tmp_subview = overlap ? new subview(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.n_rows, x_in.n_cols) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "subtraction");
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
          Mat<eT>& A = *(t.m_ptr);
    const Mat<eT>& B = x.m;
    
    const uword row_A = t.aux_row1;
    const uword row_B = x.aux_row1;
    
    const uword start_col_A = t.aux_col1;
    const uword start_col_B = x.aux_col1;
    
    uword i,j;
    
    for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
      {
      const eT tmp1 = B.at(row_B, start_col_B + i);
      const eT tmp2 = B.at(row_B, start_col_B + j);
      
      A.at(row_A, start_col_A + i) -= tmp1;
      A.at(row_A, start_col_A + j) -= tmp2;
      }
    
    if(i < t_n_cols)
      {
      A.at(row_A, start_col_A + i) -= B.at(row_B, start_col_B + i);
      }
    }
  else
    {
    for(uword col=0; col<t_n_cols; ++col)
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
  const subview<eT>* tmp_subview = overlap ? new subview(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.n_rows, x_in.n_cols) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "element-wise multiplication");
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
          Mat<eT>& A = *(t.m_ptr);
    const Mat<eT>& B = x.m;
    
    const uword row_A = t.aux_row1;
    const uword row_B = x.aux_row1;
    
    const uword start_col_A = t.aux_col1;
    const uword start_col_B = x.aux_col1;
    
    uword i,j;
    
    for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
      {
      const eT tmp1 = B.at(row_B, start_col_B + i);
      const eT tmp2 = B.at(row_B, start_col_B + j);
      
      A.at(row_A, start_col_A + i) *= tmp1;
      A.at(row_A, start_col_A + j) *= tmp2;
      }
    
    if(i < t_n_cols)
      {
      A.at(row_A, start_col_A + i) *= B.at(row_B, start_col_B + i);
      }
    }
  else
    {
    for(uword col=0; col<t_n_cols; ++col)
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
  const subview<eT>* tmp_subview = overlap ? new subview(*tmp_mat, x_in.aux_row1, x_in.aux_col1, x_in.n_rows, x_in.n_cols) : 0;
  const subview<eT>&           x = overlap ? (*tmp_subview) : x_in;
  
  subview<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, "element-wise division");
  
  const uword t_n_rows = t.n_rows;
  const uword t_n_cols = t.n_cols;
  
  if(t_n_rows == 1)
    {
          Mat<eT>& A = *(t.m_ptr);
    const Mat<eT>& B = x.m;
    
    const uword row_A = t.aux_row1;
    const uword row_B = x.aux_row1;
    
    const uword start_col_A = t.aux_col1;
    const uword start_col_B = x.aux_col1;
    
    uword i,j;
    
    for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
      {
      const eT tmp1 = B.at(row_B, start_col_B + i);
      const eT tmp2 = B.at(row_B, start_col_B + j);
      
      A.at(row_A, start_col_A + i) /= tmp1;
      A.at(row_A, start_col_A + j) /= tmp2;
      }
    
    if(i < t_n_cols)
      {
      A.at(row_A, start_col_A + i) /= B.at(row_B, start_col_B + i);
      }
    }
  else
    {
    for(uword col=0; col<t_n_cols; ++col)
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
  
  const uword local_n_cols = n_cols;
  const uword local_n_rows = n_rows;
  
  if(local_n_rows == 1)
    {
    Mat<eT>& X = (*m_ptr);
    
    const uword row           = aux_row1;
    const uword start_col     = aux_col1;
    const uword end_col_plus1 = start_col + local_n_cols;
    
    uword i,j;
    
    for(i=start_col, j=start_col+1; j < end_col_plus1; i+=2, j+=2)
      {
      X.at(row, i) = val;
      X.at(row, j) = val;
      }
    
    if(i < end_col_plus1)
      {
      X.at(row, i) = val;
      }
    }
  else
    {
    for(uword col=0; col<local_n_cols; ++col)
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
inline
void
subview<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(1));
  }



template<typename eT>
inline
void
subview<eT>::eye()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(0));
  
  const uword N = (std::min)(n_rows, n_cols);
  
  for(uword i=0; i<N; ++i)
    {
    at(i,i) = eT(1);
    }
  }



template<typename eT>
inline
eT&
subview<eT>::operator[](const uword i)
  {
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
    
  const uword index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
inline
eT
subview<eT>::operator[](const uword i) const
  {
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  const uword index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
inline
eT&
subview<eT>::operator()(const uword i)
  {
  arma_debug_check( (i >= n_elem), "subview::operator(): index out of bounds");
    
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  const uword index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
inline
eT
subview<eT>::operator()(const uword i) const
  {
  arma_debug_check( (i >= n_elem), "subview::operator(): index out of bounds");
  
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  const uword index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
inline
eT&
subview<eT>::operator()(const uword in_row, const uword in_col)
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "subview::operator(): index out of bounds");
  
  const uword index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
inline
eT
subview<eT>::operator()(const uword in_row, const uword in_col) const
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "subview::operator(): index out of bounds");
  
  const uword index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
inline
eT&
subview<eT>::at(const uword in_row, const uword in_col)
  {
  const uword index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return access::rw( (*m_ptr).mem[index] );
  }



template<typename eT>
inline
eT
subview<eT>::at(const uword in_row, const uword in_col) const
  {
  const uword index = (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT*
subview<eT>::colptr(const uword in_col)
  {
  return & access::rw((*m_ptr).mem[ (in_col + aux_col1)*m.n_rows + aux_row1 ]);
  }



template<typename eT>
arma_inline
const eT*
subview<eT>::colptr(const uword in_col) const
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
subview<eT>::is_vec() const
  {
  return ( (n_rows == 1) || (n_cols == 1) );
  }



//! X = Y.submat(...)
template<typename eT>
inline
void
subview<eT>::extract(Mat<eT>& out, const subview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: we're assuming that the matrix has already been set to the correct size and there is no aliasing;
  // size setting and alias checking is done by either the Mat contructor or operator=()
  
  const uword n_rows = in.n_rows;  // number of rows in the subview
  const uword n_cols = in.n_cols;  // number of columns in the subview
  
  arma_extra_debug_print(arma_boost::format("out.n_rows = %d   out.n_cols = %d    in.m.n_rows = %d  in.m.n_cols = %d") % out.n_rows % out.n_cols % in.m.n_rows % in.m.n_cols );
  
  
  if(in.is_vec() == true)
    {
    if(n_cols == 1)   // a column vector
      {
      arma_extra_debug_print("subview::extract(): copying col (going across rows)");
      
      // in.colptr(0) the first column of the subview, taking into account any row offset
      arrayops::copy( out.memptr(), in.colptr(0), n_rows );
      }
    else   // a row vector (possibly empty)
      {
      arma_extra_debug_print("subview::extract(): copying row (going across columns)");
      
      const Mat<eT>& X = in.m;
      
      eT* out_mem = out.memptr();
      
      const uword row       = in.aux_row1;
      const uword start_col = in.aux_col1;
      
      uword i,j;
      
      for(i=0, j=1; j < n_cols; i+=2, j+=2)
        {
        const eT tmp1 = X.at(row, start_col+i);
        const eT tmp2 = X.at(row, start_col+j);
        
        out_mem[i] = tmp1;
        out_mem[j] = tmp2;
        }
      
      if(i < n_cols)
        {
        out_mem[i] = X.at(row, start_col+i);
        }
      }
    }
  else   // general submatrix
    {
    arma_extra_debug_print("subview::extract(): general submatrix");
    
    for(uword col = 0; col<n_cols; ++col)   
      {
      arrayops::copy( out.colptr(col), in.colptr(col), n_rows );
      }
    }
  }



//! X += Y.submat(...)
template<typename eT>
inline
void
subview<eT>::plus_inplace(Mat<eT>& out, const subview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "addition");
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows == 1)
    {
    eT* out_mem = out.memptr();
    
    const Mat<eT>& X = in.m;
    
    const uword row       = in.aux_row1;
    const uword start_col = in.aux_col1;
    
    uword i,j;
    for(i=0, j=1; j < n_cols; i+=2, j+=2)
      {
      const eT tmp1 = X.at(row, start_col+i);
      const eT tmp2 = X.at(row, start_col+j);
        
      out_mem[i] += tmp1;
      out_mem[j] += tmp2;
      }
    
    if(i < n_cols)
      {
      out_mem[i] += X.at(row, start_col+i);
      }
    }
  else
    {
    for(uword col=0; col<n_cols; ++col)
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
  
  arma_debug_assert_same_size(out, in, "subtraction");
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows == 1)
    {
    eT* out_mem = out.memptr();
    
    const Mat<eT>& X = in.m;
    
    const uword row       = in.aux_row1;
    const uword start_col = in.aux_col1;
    
    uword i,j;
    for(i=0, j=1; j < n_cols; i+=2, j+=2)
      {
      const eT tmp1 = X.at(row, start_col+i);
      const eT tmp2 = X.at(row, start_col+j);
        
      out_mem[i] -= tmp1;
      out_mem[j] -= tmp2;
      }
    
    if(i < n_cols)
      {
      out_mem[i] -= X.at(row, start_col+i);
      }
    }
  else
    {
    for(uword col=0; col<n_cols; ++col)
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
  
  arma_debug_assert_same_size(out, in, "element-wise multiplication");
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows == 1)
    {
    eT* out_mem = out.memptr();
    
    const Mat<eT>& X = in.m;
    
    const uword row       = in.aux_row1;
    const uword start_col = in.aux_col1;
    
    uword i,j;
    for(i=0, j=1; j < n_cols; i+=2, j+=2)
      {
      const eT tmp1 = X.at(row, start_col+i);
      const eT tmp2 = X.at(row, start_col+j);
        
      out_mem[i] *= tmp1;
      out_mem[j] *= tmp2;
      }
    
    if(i < n_cols)
      {
      out_mem[i] *= X.at(row, start_col+i);
      }
    }
  else
    {
    for(uword col=0; col<n_cols; ++col)
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
  
  arma_debug_assert_same_size(out, in, "element-wise division");
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows == 1)
    {
    eT* out_mem = out.memptr();
    
    const Mat<eT>& X = in.m;
    
    const uword row       = in.aux_row1;
    const uword start_col = in.aux_col1;
    
    uword i,j;
    for(i=0, j=1; j < n_cols; i+=2, j+=2)
      {
      const eT tmp1 = X.at(row, start_col+i);
      const eT tmp2 = X.at(row, start_col+j);
        
      out_mem[i] /= tmp1;
      out_mem[j] /= tmp2;
      }
    
    if(i < n_cols)
      {
      out_mem[i] /= X.at(row, start_col+i);
      }
    }
  else
    {
    for(uword col=0; col<n_cols; ++col)
      {
      arrayops::inplace_div(out.colptr(col), in.colptr(col), n_rows);
      }
    }
  }



//! creation of subview (row vector)
template<typename eT>
inline
subview_row<eT>
subview<eT>::row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= n_rows, "subview::row(): out of bounds" );
  
  const uword base_row = aux_row1 + row_num;
  
  return subview_row<eT>(*m_ptr, base_row, aux_col1, n_cols);
  }



//! creation of subview (row vector)
template<typename eT>
inline
const subview_row<eT>
subview<eT>::row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= n_rows, "subview::row(): out of bounds" );
  
  const uword base_row = aux_row1 + row_num;
  
  return subview_row<eT>(m, base_row, aux_col1, n_cols);
  }



template<typename eT>
inline
subview_row<eT>
subview<eT>::operator()(const uword row_num, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const bool col_all = col_span.whole;
  
  const uword local_n_cols = n_cols;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword submat_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  const uword base_col1     = aux_col1 + in_col1;  
  const uword base_row      = aux_row1 + row_num;
  
  arma_debug_check
    (
    (row_num >= n_rows)
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "subview::operator(): indices out of bounds or incorrectly used"
    );
  
  return subview_row<eT>(*m_ptr, base_row, base_col1, submat_n_cols);
  }



template<typename eT>
inline
const subview_row<eT>
subview<eT>::operator()(const uword row_num, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool col_all = col_span.whole;
  
  const uword local_n_cols = n_cols;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword submat_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  const uword base_col1     = aux_col1 + in_col1;
  const uword base_row      = aux_row1 + row_num;
  
  arma_debug_check
    (
    (row_num >= n_rows)
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "subview::operator(): indices out of bounds or incorrectly used"
    );
  
  return subview_row<eT>(m, base_row, base_col1, submat_n_cols);
  }



//! creation of subview (column vector)
template<typename eT>
inline
subview_col<eT>
subview<eT>::col(const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "subview::col(): out of bounds");
  
  const uword base_col = aux_col1 + col_num;
  
  return subview_col<eT>(*m_ptr, base_col, aux_row1, n_rows);
  }



//! creation of subview (column vector)
template<typename eT>
inline
const subview_col<eT>
subview<eT>::col(const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "subview::col(): out of bounds");
  
  const uword base_col = aux_col1 + col_num;
  
  return subview_col<eT>(m, base_col, aux_row1, n_rows);
  }



template<typename eT>
inline
subview_col<eT>
subview<eT>::operator()(const span& row_span, const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  
  const uword local_n_rows = n_rows;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword submat_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const uword base_row1       = aux_row1 + in_row1;  
  const uword base_col        = aux_col1 + col_num;
  
  arma_debug_check
    (
    (col_num >= n_cols)
    ||
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ,
    "subview::operator(): indices out of bounds or incorrectly used"
    );
  
  return subview_col<eT>(*m_ptr, base_col, base_row1, submat_n_rows);
  }



template<typename eT>
inline
const subview_col<eT>
subview<eT>::operator()(const span& row_span, const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  
  const uword local_n_rows = n_rows;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword submat_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const uword base_row1       = aux_row1 + in_row1;
  const uword base_col        = aux_col1 + col_num;
  
  arma_debug_check
    (
    (col_num >= n_cols)
    ||
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ,
    "subview::operator(): indices out of bounds or incorrectly used"
    );
  
  return subview_col<eT>(m, base_col, base_row1, submat_n_rows);
  }



//! create a Col object which uses memory from an existing matrix object.
//! this approach is currently not alias safe
//! and does not take into account that the parent matrix object could be deleted.
//! if deleted memory is accessed by the created Col object,
//! it will cause memory corruption and/or a crash
template<typename eT>
inline
Col<eT>
subview<eT>::unsafe_col(const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "subview::unsafe_col(): out of bounds");
  
  return Col<eT>(colptr(col_num), n_rows, false, true);
  }



//! create a Col object which uses memory from an existing matrix object.
//! this approach is currently not alias safe
//! and does not take into account that the parent matrix object could be deleted.
//! if deleted memory is accessed by the created Col object,
//! it will cause memory corruption and/or a crash
template<typename eT>
inline
const Col<eT>
subview<eT>::unsafe_col(const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "subview::unsafe_col(): out of bounds");
  
  return Col<eT>(const_cast<eT*>(colptr(col_num)), n_rows, false, true);
  }



//! creation of subview (submatrix comprised of specified row vectors)
template<typename eT>
inline
subview<eT>
subview<eT>::rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "subview::rows(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  const uword base_row1 = aux_row1 + in_row1;
  
  return subview<eT>(*m_ptr, base_row1, aux_col1, subview_n_rows, n_cols );
  }



//! creation of subview (submatrix comprised of specified row vectors)
template<typename eT>
inline
const subview<eT>
subview<eT>::rows(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "subview::rows(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  const uword base_row1 = aux_row1 + in_row1;
  
  return subview<eT>(m, base_row1, aux_col1, subview_n_rows, n_cols );
  }



//! creation of subview (submatrix comprised of specified column vectors)
template<typename eT>
inline
subview<eT>
subview<eT>::cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "subview::cols(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  const uword base_col1 = aux_col1 + in_col1;
  
  return subview<eT>(*m_ptr, aux_row1, base_col1, n_rows, subview_n_cols);
  }



//! creation of subview (submatrix comprised of specified column vectors)
template<typename eT>
inline
const subview<eT>
subview<eT>::cols(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "subview::cols(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  const uword base_col1 = aux_col1 + in_col1;
  
  return subview<eT>(m, aux_row1, base_col1, n_rows, subview_n_cols);
  }



//! creation of subview (submatrix)
template<typename eT>
inline
subview<eT>
subview<eT>::submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 >  in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "subview::submat(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  const uword base_row1 = aux_row1 + in_row1;
  const uword base_col1 = aux_col1 + in_col1;
  
  return subview<eT>(*m_ptr, base_row1, base_col1, subview_n_rows, subview_n_cols);
  }



//! creation of subview (generic submatrix)
template<typename eT>
inline
const subview<eT>
subview<eT>::submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 >  in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "subview::submat(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  const uword base_row1 = aux_row1 + in_row1;
  const uword base_col1 = aux_col1 + in_col1;
  
  return subview<eT>(m, base_row1, base_col1, subview_n_rows, subview_n_cols);
  }



//! creation of subview (submatrix)
template<typename eT>
inline
subview<eT>
subview<eT>::submat(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  const bool col_all = col_span.whole;
  
  const uword local_n_rows = n_rows;
  const uword local_n_cols = n_cols;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword submat_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword submat_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  arma_debug_check
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "subview::submat(): indices out of bounds or incorrectly used"
    );
  
  const uword base_row1 = aux_row1 + in_row1;
  const uword base_col1 = aux_col1 + in_col1;
  
  return subview<eT>(*m_ptr, base_row1, base_col1, submat_n_rows, submat_n_cols);
  }



//! creation of subview (generic submatrix)
template<typename eT>
inline
const subview<eT>
subview<eT>::submat(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  const bool col_all = col_span.whole;
  
  const uword local_n_rows = n_rows;
  const uword local_n_cols = n_cols;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword submat_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword submat_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  arma_debug_check
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "subview::submat(): indices out of bounds or incorrectly used"
    );
  
  const uword base_row1 = aux_row1 + in_row1;
  const uword base_col1 = aux_col1 + in_col1;
  
  return subview<eT>(m, base_row1, base_col1, submat_n_rows, submat_n_cols);
  }



template<typename eT>
inline
subview<eT>
subview<eT>::operator()(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return (*this).submat(row_span, col_span);
  }



template<typename eT>
inline
const subview<eT>
subview<eT>::operator()(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).submat(row_span, col_span);
  }



//! creation of diagview (diagonal)
template<typename eT>
inline
diagview<eT>
subview<eT>::diag(const sword in_id)
  {
  arma_extra_debug_sigprint();
  
  const uword row_offset = (in_id < 0) ? uword(-in_id) : 0;
  const uword col_offset = (in_id > 0) ? uword( in_id) : 0;
  
  arma_debug_check
    (
    ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
    "subview::diag(): requested diagonal out of bounds"
    );
  
  const uword len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  const uword base_row_offset = aux_row1 + row_offset;
  const uword base_col_offset = aux_col1 + col_offset;
  
  return diagview<eT>(*m_ptr, base_row_offset, base_col_offset, len);
  }



//! creation of diagview (diagonal)
template<typename eT>
inline
const diagview<eT>
subview<eT>::diag(const sword in_id) const
  {
  arma_extra_debug_sigprint();
  
  const uword row_offset = (in_id < 0) ? -in_id : 0;
  const uword col_offset = (in_id > 0) ?  in_id : 0;
  
  arma_debug_check
    (
    ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
    "subview::diag(): requested diagonal out of bounds"
    );
  
  const uword len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  const uword base_row_offset = aux_row1 + row_offset;
  const uword base_col_offset = aux_col1 + col_offset;
  
  return diagview<eT>(m, base_row_offset, base_col_offset, len);
  }



template<typename eT>
inline
void
subview<eT>::swap_rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 >= n_rows) || (in_row2 >= n_rows),
    "subview::swap_rows(): out of bounds"
    );
  
  eT* mem = (*m_ptr).memptr();
  
  for(uword col=0; col<n_cols; ++col)
    {
    const uword offset = (aux_col1 + col) * m.n_rows;
    const uword pos1   = aux_row1 + in_row1 + offset;
    const uword pos2   = aux_row1 + in_row2 + offset;
    
    const eT tmp          = mem[pos1];
    access::rw(mem[pos1]) = mem[pos2];
    access::rw(mem[pos2]) = tmp;
    }
  }



template<typename eT>
inline
void
subview<eT>::swap_cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 >= n_cols) || (in_col2 >= n_cols),
    "subview::swap_cols(): out of bounds"
    );
  
  if(n_elem > 0)
    {
    eT* ptr1 = colptr(in_col1);
    eT* ptr2 = colptr(in_col2);
    
    for(uword row=0; row<n_rows; ++row)
      {
      const eT tmp = ptr1[row];
      ptr1[row]    = ptr2[row];
      ptr2[row]    = tmp;
      }
    }
  }



// template<typename eT>
// inline
// subview<eT>::iter::iter(const subview<eT>& S)
//   : mem       (S.m.mem)
//   , n_rows    (S.m.n_rows)
//   , row_start (S.aux_row1)
//   , row_end_p1(row_start + S.n_rows)
//   , row       (row_start)
//   , col       (S.aux_col1)
//   , i         (row + col*n_rows)
//   {
//   arma_extra_debug_sigprint();
//   }
// 
// 
// 
// template<typename eT>
// arma_inline
// eT
// subview<eT>::iter::operator*() const
//   {
//   return mem[i];
//   }
// 
// 
// 
// template<typename eT>
// inline
// void
// subview<eT>::iter::operator++()
//   {
//   ++row;
//   
//   if(row < row_end_p1)
//     {
//     ++i;
//     }
//   else
//     {
//     row = row_start;
//     ++col;
//     
//     i = row + col*n_rows;
//     }
//   }
// 
// 
// 
// template<typename eT>
// inline
// void
// subview<eT>::iter::operator++(int)
//   {
//   operator++();
//   }



//
//
//



template<typename eT>
inline
subview_col<eT>::subview_col(const Mat<eT>& in_m, const uword in_col)
  : subview<eT>(in_m, 0, in_col, in_m.n_rows, 1)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_col<eT>::subview_col(Mat<eT>& in_m, const uword in_col)
  : subview<eT>(in_m, 0, in_col, in_m.n_rows, 1)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_col<eT>::subview_col(const Mat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows)
  : subview<eT>(in_m, in_row1, in_col, in_n_rows, 1)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_col<eT>::subview_col(Mat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows)
  : subview<eT>(in_m, in_row1, in_col, in_n_rows, 1)
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



template<typename eT>
inline
subview_col<eT>
subview_col<eT>::rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= subview<eT>::n_rows) ), "subview_col::rows(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  const uword base_row1 = this->aux_row1 + in_row1;
  
  return subview_col<eT>(*(this->m_ptr), this->aux_col1, base_row1, subview_n_rows);
  }



template<typename eT>
inline
const subview_col<eT>
subview_col<eT>::rows(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= subview<eT>::n_rows) ), "subview_col::rows(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  const uword base_row1 = this->aux_row1 + in_row1;
  
  return subview_col<eT>(this->m, this->aux_col1, base_row1, subview_n_rows);
  }



template<typename eT>
inline
subview_col<eT>
subview_col<eT>::subvec(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= subview<eT>::n_rows) ), "subview_col::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  const uword base_row1 = this->aux_row1 + in_row1;
  
  return subview_col<eT>(*(this->m_ptr), this->aux_col1, base_row1, subview_n_rows);
  }



template<typename eT>
inline
const subview_col<eT>
subview_col<eT>::subvec(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_row1 > in_row2) || (in_row2 >= subview<eT>::n_rows) ), "subview_col::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  const uword base_row1 = this->aux_row1 + in_row1;
  
  return subview_col<eT>(this->m, this->aux_col1, base_row1, subview_n_rows);
  }



//
//
//



template<typename eT>
inline
subview_row<eT>::subview_row(const Mat<eT>& in_m, const uword in_row)
  : subview<eT>(in_m, in_row, 0, 1, in_m.n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_row<eT>::subview_row(Mat<eT>& in_m, const uword in_row)
  : subview<eT>(in_m, in_row, 0, 1, in_m.n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_row<eT>::subview_row(const Mat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols)
  : subview<eT>(in_m, in_row, in_col1, 1, in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_row<eT>::subview_row(Mat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols)
  : subview<eT>(in_m, in_row, in_col1, 1, in_n_cols)
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



template<typename eT>
inline
subview_row<eT>
subview_row<eT>::cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= subview<eT>::n_cols) ), "subview_row::cols(): indices out of bounds or incorrectly used" );
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  const uword base_col1 = this->aux_col1 + in_col1;
  
  return subview_row<eT>(*(this->m_ptr), this->aux_row1, base_col1, subview_n_cols);
  }



template<typename eT>
inline
const subview_row<eT>
subview_row<eT>::cols(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= subview<eT>::n_cols) ), "subview_row::cols(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  const uword base_col1 = this->aux_col1 + in_col1;
  
  return subview_row<eT>(this->m, this->aux_row1, base_col1, subview_n_cols);
  }



template<typename eT>
inline
subview_row<eT>
subview_row<eT>::subvec(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= subview<eT>::n_cols) ), "subview_row::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  const uword base_col1 = this->aux_col1 + in_col1;
  
  return subview_row<eT>(*(this->m_ptr), this->aux_row1, base_col1, subview_n_cols);
  }



template<typename eT>
inline
const subview_row<eT>
subview_row<eT>::subvec(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in_col1 > in_col2) || (in_col2 >= subview<eT>::n_cols) ), "subview_row::subvec(): indices out of bounds or incorrectly used");
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  const uword base_col1 = this->aux_col1 + in_col1;
  
  return subview_row<eT>(this->m, this->aux_row1, base_col1, subview_n_cols);
  }



//! @}
