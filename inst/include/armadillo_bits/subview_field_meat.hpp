// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup subview_field
//! @{


template<typename oT>
inline
subview_field<oT>::~subview_field()
  {
  arma_extra_debug_sigprint();
  }



template<typename oT>
arma_inline
subview_field<oT>::subview_field
  (
  const field<oT>& in_f,
  const uword      in_row1,
  const uword      in_col1,
  const uword      in_n_rows,
  const uword      in_n_cols
  )
  : f(in_f)
  , aux_row1(in_row1)
  , aux_col1(in_col1)
  , n_rows(in_n_rows)
  , n_cols(in_n_cols)
  , n_elem(in_n_rows*in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename oT>
inline
void
subview_field<oT>::operator= (const field<oT>& x)
  {
  arma_extra_debug_sigprint();
  
  subview_field<oT>& t = *this;
  
  arma_debug_check( (t.n_rows != x.n_rows) || (t.n_cols != x.n_cols), "incompatible field dimensions");
  
  for(uword col=0; col<t.n_cols; ++col)
    {
    for(uword row=0; row<t.n_rows; ++row)
      {
      t.at(row,col) = x.at(row,col);
      }
    }
  }



//! x.subfield(...) = y.subfield(...)
template<typename oT>
inline
void
subview_field<oT>::operator= (const subview_field<oT>& x_in)
  {
  arma_extra_debug_sigprint();
  
  const bool overlap = check_overlap(x_in);
        
        field<oT>*         tmp_field   = overlap ? new field<oT>(x_in.f) : 0;
  const subview_field<oT>* tmp_subview = overlap ? new subview_field<oT>(*tmp_field, x_in.aux_row1, x_in.aux_col1, x_in.n_rows, x_in.n_cols) : 0;
  const subview_field<oT>& x           = overlap ? (*tmp_subview) : x_in;
  
  subview_field<oT>& t = *this;
  
  arma_debug_check( (t.n_rows != x.n_rows) || (t.n_cols != x.n_cols), "incompatible field dimensions");
  
  for(uword col=0; col<t.n_cols; ++col)
    {
    for(uword row=0; row<t.n_rows; ++row)
      {
      t.at(row,col) = x.at(row,col);
      }
    }
    
  if(overlap)
    {
    delete tmp_subview;
    delete tmp_field;
    }
  }



template<typename oT>
arma_inline
oT&
subview_field<oT>::operator[](const uword i)
  {
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
    
  const uword index = (in_col + aux_col1)*f.n_rows + aux_row1 + in_row;
  
  return *((const_cast< field<oT>& >(f)).mem[index]);
  }



template<typename oT>
arma_inline
const oT&
subview_field<oT>::operator[](const uword i) const
  {
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  const uword index = (in_col + aux_col1)*f.n_rows + aux_row1 + in_row;
  
  return *(f.mem[index]);
  }



template<typename oT>
arma_inline
oT&
subview_field<oT>::operator()(const uword i)
  {
  arma_debug_check( (i >= n_elem), "subview_field::operator(): index out of bounds");
  
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  const uword index = (in_col + aux_col1)*f.n_rows + aux_row1 + in_row;
  
  return *((const_cast< field<oT>& >(f)).mem[index]);
  }



template<typename oT>
arma_inline
const oT&
subview_field<oT>::operator()(const uword i) const
  {
  arma_debug_check( (i >= n_elem), "subview_field::operator(): index out of bounds");
  
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  const uword index = (in_col + aux_col1)*f.n_rows + aux_row1 + in_row;
  
  return *(f.mem[index]);
  }



template<typename oT>
arma_inline
oT&
subview_field<oT>::operator()(const uword in_row, const uword in_col)
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "subview_field::operator(): index out of bounds");
  
  const uword index = (in_col + aux_col1)*f.n_rows + aux_row1 + in_row;
  
  return *((const_cast< field<oT>& >(f)).mem[index]);
  }



template<typename oT>
arma_inline
const oT&
subview_field<oT>::operator()(const uword in_row, const uword in_col) const
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "subview_field::operator(): index out of bounds");
  
  const uword index = (in_col + aux_col1)*f.n_rows + aux_row1 + in_row;
  
  return *(f.mem[index]);
  }



template<typename oT>
arma_inline
oT&
subview_field<oT>::at(const uword in_row, const uword in_col)
  {
  const uword index = (in_col + aux_col1)*f.n_rows + aux_row1 + in_row;
  
  return *((const_cast< field<oT>& >(f)).mem[index]);
  }



template<typename oT>
arma_inline
const oT&
subview_field<oT>::at(const uword in_row, const uword in_col) const
  {
  //arma_extra_debug_sigprint();
  
  const uword index = (in_col + aux_col1)*f.n_rows + aux_row1 + in_row;
  
  return *(f.mem[index]);
  }



template<typename oT>
inline
bool
subview_field<oT>::check_overlap(const subview_field<oT>& x) const
  {
  const subview_field<oT>& t = *this;
  
  if(&t.f != &x.f)
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



//! X = Y.subfield(...)
template<typename oT>
inline
void
subview_field<oT>::extract(field<oT>& actual_out, const subview_field<oT>& in)
  {
  arma_extra_debug_sigprint();
  
  //
  const bool alias = (&actual_out == &in.f);
  
  field<oT>* tmp = (alias) ? new field<oT> : 0;
  field<oT>& out = (alias) ? (*tmp)        : actual_out;
  
  //
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  out.set_size(n_rows, n_cols);
  
  arma_extra_debug_print(arma_boost::format("out.n_rows = %d   out.n_cols = %d    in.m.n_rows = %d  in.m.n_cols = %d") % out.n_rows % out.n_cols % in.f.n_rows % in.f.n_cols );
  
  for(uword col = 0; col<n_cols; ++col)
    {
    for(uword row = 0; row<n_rows; ++row)
      {
      out.at(row,col) = in.at(row,col);
      }
    }
  
  
  if(alias)
    {
    actual_out = out;
    delete tmp;
    }
  
  }



//! @}
