// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup field
//! @{


template<typename oT>
inline
field<oT>::~field()
  {
  arma_debug_sigprint_this(this);
  
  delete_objects();
  
  if(n_elem > 0)  { delete [] mem; }
  
  // try to expose buggy user code that accesses deleted objects
  mem = nullptr;
  }



template<typename oT>
inline
field<oT>::field()
  : n_rows(0)
  , n_cols(0)
  , n_slices(0)
  , n_elem(0)
  , mem(nullptr)
  {
  arma_debug_sigprint_this(this);
  }



//! construct a field from a given field
template<typename oT>
inline
field<oT>::field(const field& x)
  : n_rows(0)
  , n_cols(0)
  , n_slices(0)
  , n_elem(0)
  , mem(nullptr)
  {
  arma_debug_sigprint(arma_str::format("this: %x; x: %x") % this % &x);
  
  init(x);
  }



//! construct a field from a given field
template<typename oT>
inline
field<oT>&
field<oT>::operator=(const field& x)
  {
  arma_debug_sigprint();
  
  init(x);
  
  return *this;
  }



//! construct a field from subview_field (eg. construct a field from a delayed subfield operation)
template<typename oT>
inline
field<oT>::field(const subview_field<oT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_slices(0)
  , n_elem(0)
  , mem(nullptr)
  {
  arma_debug_sigprint_this(this);
  
  this->operator=(X);
  }



//! construct a field from subview_field (eg. construct a field from a delayed subfield operation)
template<typename oT>
inline
field<oT>&
field<oT>::operator=(const subview_field<oT>& X)
  {
  arma_debug_sigprint();
  
  subview_field<oT>::extract(*this, X);
  
  return *this;
  }



//! construct the field with the specified number of elements,
//! assuming a column-major layout
template<typename oT>
inline
field<oT>::field(const uword n_elem_in)
  : n_rows(0)
  , n_cols(0)
  , n_slices(0)
  , n_elem(0)
  , mem(nullptr)
  {
  arma_debug_sigprint_this(this);
  
  init(n_elem_in, 1);
  }



//! construct the field with the specified dimensions
template<typename oT>
inline
field<oT>::field(const uword n_rows_in, const uword n_cols_in)
  : n_rows(0)
  , n_cols(0)
  , n_slices(0)
  , n_elem(0)
  , mem(nullptr)
  {
  arma_debug_sigprint_this(this);
  
  init(n_rows_in, n_cols_in);
  }



//! construct the field with the specified dimensions
template<typename oT>
inline
field<oT>::field(const uword n_rows_in, const uword n_cols_in, const uword n_slices_in)
  : n_rows(0)
  , n_cols(0)
  , n_slices(0)
  , n_elem(0)
  , mem(nullptr)
  {
  arma_debug_sigprint_this(this);
  
  init(n_rows_in, n_cols_in, n_slices_in);
  }



template<typename oT>
inline
field<oT>::field(const SizeMat& s)
  : n_rows(0)
  , n_cols(0)
  , n_slices(0)
  , n_elem(0)
  , mem(nullptr)
  {
  arma_debug_sigprint_this(this);
  
  init(s.n_rows, s.n_cols);
  }



template<typename oT>
inline
field<oT>::field(const SizeCube& s)
  : n_rows(0)
  , n_cols(0)
  , n_slices(0)
  , n_elem(0)
  , mem(nullptr)
  {
  arma_debug_sigprint_this(this);
  
  init(s.n_rows, s.n_cols, s.n_slices);
  }



//! change the field to have the specified number of elements,
//! assuming a column-major layout (data is not preserved)
template<typename oT>
inline
field<oT>&
field<oT>::set_size(const uword n_elem_in)
  {
  arma_debug_sigprint(arma_str::format("n_elem_in: %u") % n_elem_in);
  
  init(n_elem_in, 1);
  
  return *this;
  }



//! change the field to have the specified dimensions (data is not preserved)
template<typename oT>
inline
field<oT>&
field<oT>::set_size(const uword n_rows_in, const uword n_cols_in)
  {
  arma_debug_sigprint(arma_str::format("n_rows_in: %u; n_cols_in: %u") % n_rows_in % n_cols_in);
  
  init(n_rows_in, n_cols_in);
  
  return *this;
  }



//! change the field to have the specified dimensions (data is not preserved)
template<typename oT>
inline
field<oT>&
field<oT>::set_size(const uword n_rows_in, const uword n_cols_in, const uword n_slices_in)
  {
  arma_debug_sigprint(arma_str::format("n_rows_in: %u; n_cols_in: %u; n_slices_in: %u") % n_rows_in % n_cols_in % n_slices_in);
  
  init(n_rows_in, n_cols_in, n_slices_in);
  
  return *this;
  }



template<typename oT>
inline
field<oT>&
field<oT>::set_size(const SizeMat& s)
  {
  arma_debug_sigprint();
  
  init(s.n_rows, s.n_cols);
  
  return *this;
  }



template<typename oT>
inline
field<oT>&
field<oT>::set_size(const SizeCube& s)
  {
  arma_debug_sigprint();
  
  init(s.n_rows, s.n_cols, s.n_slices);
  
  return *this;
  }



template<typename oT>
inline
field<oT>&
field<oT>::reshape(const uword n_elem_in)
  {
  arma_debug_sigprint();
  
  return (*this).reshape(n_elem_in, 1, 1);
  }



template<typename oT>
inline
field<oT>&
field<oT>::reshape(const uword n_rows_in, const uword n_cols_in)
  {
  arma_debug_sigprint();
  
  return (*this).reshape(n_rows_in, n_cols_in, 1);
  }



template<typename oT>
inline
field<oT>&
field<oT>::reshape(const uword n_rows_in, const uword n_cols_in, const uword n_slices_in)
  {
  arma_debug_sigprint(arma_str::format("n_rows_in: %u; n_cols_in: %u; n_slices_in: %u") % n_rows_in % n_cols_in % n_slices_in);
  
  if((n_rows == n_rows_in) && (n_cols == n_cols_in) && (n_slices == n_slices_in))
    {
    // do nothing
    }
  else
  if((n_elem == 0) || ((n_rows == n_cols_in) && (n_cols == n_rows_in) && (n_slices == n_slices_in)))
    {
    init(n_rows_in, n_cols_in, n_slices_in);
    }
  else
    {
    field<oT> tmp(n_rows_in, n_cols_in, n_slices_in);
    
    const uword n_elem_to_copy = (std::min)((*this).n_elem, tmp.n_elem);
    
    for(uword i=0; i < n_elem_to_copy; ++i)  { tmp.at(i) = std::move((*this).at(i)); }
    
    (*this) = std::move(tmp);
    }
  
  return *this;
  }



template<typename oT>
inline
field<oT>&
field<oT>::reshape(const SizeMat& s)
  {
  arma_debug_sigprint();
  
  return (*this).reshape(s.n_rows, s.n_cols, 1);
  }



template<typename oT>
inline
field<oT>&
field<oT>::reshape(const SizeCube& s)
  {
  arma_debug_sigprint();
  
  return (*this).reshape(s.n_rows, s.n_cols, s.n_slices);
  }



template<typename oT>
inline
field<oT>&
field<oT>::resize(const uword n_elem_in)
  {
  arma_debug_sigprint();
  
  return (*this).resize(n_elem_in, 1, 1);
  }



template<typename oT>
inline
field<oT>&
field<oT>::resize(const uword n_rows_in, const uword n_cols_in)
  {
  arma_debug_sigprint();
  
  return (*this).resize(n_rows_in, n_cols_in, 1);
  }



template<typename oT>
inline
field<oT>&
field<oT>::resize(const uword n_rows_in, const uword n_cols_in, const uword n_slices_in)
  {
  arma_debug_sigprint(arma_str::format("n_rows_in: %u; n_cols_in: %u; n_slices_in: %u") % n_rows_in % n_cols_in % n_slices_in);
  
  if((n_rows == n_rows_in) && (n_cols == n_cols_in) && (n_slices == n_slices_in))
    {
    // do nothing
    }
  else
  if(n_elem == 0)
    {
    (*this).set_size(n_rows_in, n_cols_in, n_slices_in);
    }
  else
    {
    // better-than-nothing implementation
    
    field<oT> tmp(n_rows_in, n_cols_in, n_slices_in);
    
    if(tmp.n_elem > 0)
      {
      const uword end_row   = (std::min)(n_rows_in,   n_rows  ) - 1;
      const uword end_col   = (std::min)(n_cols_in,   n_cols  ) - 1;
      const uword end_slice = (std::min)(n_slices_in, n_slices) - 1;
      
      tmp.subfield(0, 0, 0, end_row, end_col, end_slice) = (*this).subfield(0, 0, 0, end_row, end_col, end_slice);
      }
    
    (*this) = std::move(tmp);
    }
  
  return *this;
  }



template<typename oT>
inline
field<oT>&
field<oT>::resize(const SizeMat& s)
  {
  arma_debug_sigprint();
  
  return (*this).resize(s.n_rows, s.n_cols, 1);
  }



template<typename oT>
inline
field<oT>&
field<oT>::resize(const SizeCube& s)
  {
  arma_debug_sigprint();
  
  return (*this).resize(s.n_rows, s.n_cols, s.n_slices);
  }



template<typename oT>
inline
field<oT>::field(const std::vector<oT>& x)
  : n_rows  (0)
  , n_cols  (0)
  , n_slices(0)
  , n_elem  (0)
  , mem     (nullptr)
  {
  arma_debug_sigprint_this(this);
  
  (*this).operator=(x);
  }



template<typename oT>
inline
field<oT>&
field<oT>::operator=(const std::vector<oT>& x)
  {
  arma_debug_sigprint();
  
  const uword N = uword(x.size());
  
  set_size(N, 1);
  
  for(uword i=0; i<N; ++i)  { operator[](i) = x[i]; }
  
  return *this;
  }



template<typename oT>
inline
field<oT>::field(const std::initializer_list<oT>& list)
  : n_rows  (0)
  , n_cols  (0)
  , n_slices(0)
  , n_elem  (0)
  , mem     (nullptr)
  {
  arma_debug_sigprint_this(this);
  
  (*this).operator=(list);
  }



template<typename oT>
inline
field<oT>&
field<oT>::operator=(const std::initializer_list<oT>& list)
  {
  arma_debug_sigprint();
  
  const uword N = uword(list.size());
  
  set_size(1, N);
  
  const oT* item_ptr = list.begin();
  
  for(uword i=0; i<N; ++i)  { operator[](i) = item_ptr[i]; }
  
  return *this;
  }



template<typename oT>
inline
field<oT>::field(const std::initializer_list< std::initializer_list<oT> >& list)
  : n_rows  (0)
  , n_cols  (0)
  , n_slices(0)
  , n_elem  (0)
  , mem     (nullptr)
  {
  arma_debug_sigprint_this(this);
  
  (*this).operator=(list);
  }



template<typename oT>
inline
field<oT>&
field<oT>::operator=(const std::initializer_list< std::initializer_list<oT> >& list)
  {
  arma_debug_sigprint();
  
  uword x_n_rows = uword(list.size());
  uword x_n_cols = 0;
  
  auto it     = list.begin();
  auto it_end = list.end();
  
  for(; it != it_end; ++it)  { x_n_cols = (std::max)(x_n_cols, uword((*it).size())); }
  
  field<oT>& t = (*this);
  
  t.set_size(x_n_rows, x_n_cols);
  
  uword row_num = 0;
  
  auto row_it     = list.begin();
  auto row_it_end = list.end();
  
  for(; row_it != row_it_end; ++row_it)
    {
    uword col_num = 0;
    
    auto col_it     = (*row_it).begin();
    auto col_it_end = (*row_it).end();
    
    for(; col_it != col_it_end; ++col_it)
      {
      t.at(row_num, col_num) = (*col_it);
      ++col_num;
      }
    
    for(uword c=col_num; c < x_n_cols; ++c)
      {
      t.at(row_num, c) = oT();
      }
    
    ++row_num;
    }
  
  return *this;
  }



template<typename oT>
inline
field<oT>::field(field<oT>&& X)
  : n_rows  (X.n_rows  )
  , n_cols  (X.n_cols  )
  , n_slices(X.n_slices)
  , n_elem  (X.n_elem  )
  , mem     (X.mem     )
  {
  arma_debug_sigprint(arma_str::format("this: %x; X: %x") % this % &X);
  
  access::rw(X.n_rows  ) = 0;
  access::rw(X.n_cols  ) = 0;
  access::rw(X.n_slices) = 0;
  access::rw(X.n_elem  ) = 0;
  access::rw(X.mem     ) = nullptr;
  }



template<typename oT>
inline
field<oT>&
field<oT>::operator=(field<oT>&& X)
  {
  arma_debug_sigprint(arma_str::format("this: %x; X: %x") % this % &X);
  
  if(this == &X)  { return *this; }
  
  reset();
  
  access::rw(n_rows  ) = X.n_rows;
  access::rw(n_cols  ) = X.n_cols;
  access::rw(n_slices) = X.n_slices;
  access::rw(n_elem  ) = X.n_elem;
  
  mem = X.mem;
  
  access::rw(X.n_rows  ) = 0;
  access::rw(X.n_cols  ) = 0;
  access::rw(X.n_elem  ) = 0;
  access::rw(X.n_slices) = 0;
  access::rw(X.mem     ) = nullptr;
  
  return *this;
  }



//! change the field to have the specified dimensions (data is not preserved)
template<typename oT>
template<typename oT2>
inline
field<oT>&
field<oT>::copy_size(const field<oT2>& x)
  {
  arma_debug_sigprint();
  
  init(x.n_rows, x.n_cols, x.n_slices);
  
  return *this;
  }



//! linear element accessor (treats the field as a vector); no bounds check
template<typename oT>
arma_inline
oT&
field<oT>::operator[] (const uword i)
  {
  return (*mem[i]);
  }



//! linear element accessor (treats the field as a vector); no bounds check
template<typename oT>
arma_inline
const oT&
field<oT>::operator[] (const uword i) const
  {
  return (*mem[i]);
  }



//! linear element accessor (treats the field as a vector); no bounds check
template<typename oT>
arma_inline
oT&
field<oT>::at(const uword i)
  {
  return (*mem[i]);
  }



//! linear element accessor (treats the field as a vector); no bounds check
template<typename oT>
arma_inline
const oT&
field<oT>::at(const uword i) const
  {
  return (*mem[i]);
  }



//! linear element accessor (treats the field as a vector); bounds checking not done when ARMA_DONT_CHECK_CONFORMANCE is defined
template<typename oT>
arma_inline
oT&
field<oT>::operator() (const uword i)
  {
  arma_conform_check_bounds( (i >= n_elem), "field::operator(): index out of bounds" );
  
  return (*mem[i]);
  }



//! linear element accessor (treats the field as a vector); bounds checking not done when ARMA_DONT_CHECK_CONFORMANCE is defined
template<typename oT>
arma_inline
const oT&
field<oT>::operator() (const uword i) const
  {
  arma_conform_check_bounds( (i >= n_elem), "field::operator(): index out of bounds" );
  
  return (*mem[i]);
  }



//! element accessor; bounds checking not done when ARMA_DONT_CHECK_CONFORMANCE is defined
template<typename oT>
arma_inline
oT&
field<oT>::operator() (const uword in_row, const uword in_col)
  {
  arma_conform_check_bounds( ((in_row >= n_rows) || (in_col >= n_cols) || (0 >= n_slices) ), "field::operator(): index out of bounds" );
  
  return (*mem[in_row + in_col*n_rows]);
  }



//! element accessor; bounds checking not done when ARMA_DONT_CHECK_CONFORMANCE is defined
template<typename oT>
arma_inline
const oT&
field<oT>::operator() (const uword in_row, const uword in_col) const
  {
  arma_conform_check_bounds( ((in_row >= n_rows) || (in_col >= n_cols) || (0 >= n_slices) ), "field::operator(): index out of bounds" );
  
  return (*mem[in_row + in_col*n_rows]);
  }



//! element accessor; bounds checking not done when ARMA_DONT_CHECK_CONFORMANCE is defined
template<typename oT>
arma_inline
oT&
field<oT>::operator() (const uword in_row, const uword in_col, const uword in_slice)
  {
  arma_conform_check_bounds( ((in_row >= n_rows) || (in_col >= n_cols) || (in_slice >= n_slices)), "field::operator(): index out of bounds" );
  
  return (*mem[in_row + in_col*n_rows + in_slice*(n_rows*n_cols)]);
  }



//! element accessor; bounds checking not done when ARMA_DONT_CHECK_CONFORMANCE is defined
template<typename oT>
arma_inline
const oT&
field<oT>::operator() (const uword in_row, const uword in_col, const uword in_slice) const
  {
  arma_conform_check_bounds( ((in_row >= n_rows) || (in_col >= n_cols) || (in_slice >= n_slices)), "field::operator(): index out of bounds" );
  
  return (*mem[in_row + in_col*n_rows + in_slice*(n_rows*n_cols)]);
  }



#if defined(__cpp_multidimensional_subscript)
  
  //! element accessor; no bounds check
  template<typename oT>
  arma_inline
  oT&
  field<oT>::operator[] (const uword in_row, const uword in_col)
    {
    return (*mem[in_row + in_col*n_rows]);
    }
  
  
  
  //! element accessor; no bounds check
  template<typename oT>
  arma_inline
  const oT&
  field<oT>::operator[] (const uword in_row, const uword in_col) const
    {
    return (*mem[in_row + in_col*n_rows]);
    }
  
#endif



//! element accessor; no bounds check
template<typename oT>
arma_inline
oT&
field<oT>::at(const uword in_row, const uword in_col)
  {
  return (*mem[in_row + in_col*n_rows]);
  }



//! element accessor; no bounds check
template<typename oT>
arma_inline
const oT&
field<oT>::at(const uword in_row, const uword in_col) const
  {
  return (*mem[in_row + in_col*n_rows]);
  }



#if defined(__cpp_multidimensional_subscript)
  
  //! element accessor; no bounds check
  template<typename oT>
  arma_inline
    oT&
  field<oT>::operator[] (const uword in_row, const uword in_col, const uword in_slice)
    {
    return (*mem[in_row + in_col*n_rows + in_slice*(n_rows*n_cols)]);
    }
  
  
  
  //! element accessor; no bounds check
  template<typename oT>
  arma_inline
    const oT&
  field<oT>::operator[] (const uword in_row, const uword in_col, const uword in_slice) const
    {
    return (*mem[in_row + in_col*n_rows + in_slice*(n_rows*n_cols)]);
    }
  
#endif



//! element accessor; no bounds check
template<typename oT>
arma_inline
oT&
field<oT>::at(const uword in_row, const uword in_col, const uword in_slice)
  {
  return (*mem[in_row + in_col*n_rows + in_slice*(n_rows*n_cols)]);
  }



//! element accessor; no bounds check
template<typename oT>
arma_inline
const oT&
field<oT>::at(const uword in_row, const uword in_col, const uword in_slice) const
  {
  return (*mem[in_row + in_col*n_rows + in_slice*(n_rows*n_cols)]);
  }



template<typename oT>
arma_inline
oT&
field<oT>::front()
  {
  arma_conform_check( (n_elem == 0), "field::front(): field is empty" );
  
  return (*mem[0]);
  }



template<typename oT>
arma_inline
const oT&
field<oT>::front() const
  {
  arma_conform_check( (n_elem == 0), "field::front(): field is empty" );
  
  return (*mem[0]);
  }



template<typename oT>
arma_inline
oT&
field<oT>::back()
  {
  arma_conform_check( (n_elem == 0), "field::back(): field is empty" );
  
  return (*mem[n_elem-1]);
  }



template<typename oT>
arma_inline
const oT&
field<oT>::back() const
  {
  arma_conform_check( (n_elem == 0), "field::back(): field is empty" );
  
  return (*mem[n_elem-1]);
  }



template<typename oT>
inline
field_injector< field<oT> >
field<oT>::operator<<(const oT& val)
  {
  return field_injector< field<oT> >(*this, val);
  }



template<typename oT>
inline
field_injector< field<oT> >
field<oT>::operator<<(const injector_end_of_row<>& x)
  {
  return field_injector< field<oT> >(*this, x);
  }



//! creation of subview_field (row of a field)
template<typename oT>
inline
subview_field<oT>
field<oT>::row(const uword row_num)
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::row(): field must be 2D" );

  arma_conform_check_bounds( (row_num >= n_rows), "field::row(): row out of bounds" );
  
  return subview_field<oT>(*this, row_num, 0, 1, n_cols);
  }



//! creation of subview_field (row of a field)
template<typename oT>
inline
const subview_field<oT>
field<oT>::row(const uword row_num) const
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::row(): field must be 2D" );

  arma_conform_check_bounds( (row_num >= n_rows), "field::row(): row out of bounds" );
  
  return subview_field<oT>(*this, row_num, 0, 1, n_cols);
  }



//! creation of subview_field (column of a field)
template<typename oT>
inline
subview_field<oT>
field<oT>::col(const uword col_num)
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::col(): field must be 2D" );

  arma_conform_check_bounds( (col_num >= n_cols), "field::col(): out of bounds" );
  
  return subview_field<oT>(*this, 0, col_num, n_rows, 1);
  }



//! creation of subview_field (column of a field)
template<typename oT>
inline
const subview_field<oT>
field<oT>::col(const uword col_num) const
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::col(): field must be 2D" );

  arma_conform_check_bounds( (col_num >= n_cols), "field::col(): out of bounds" );
  
  return subview_field<oT>(*this, 0, col_num, n_rows, 1);
  }



//! creation of subview_field (slice of a field)
template<typename oT>
inline
subview_field<oT>
field<oT>::slice(const uword slice_num)
  {
  arma_debug_sigprint();
  
  arma_conform_check_bounds( (slice_num >= n_slices), "field::slice(): out of bounds" );
  
  return subview_field<oT>(*this, 0, 0, slice_num, n_rows, n_cols, 1);
  }



//! creation of subview_field (slice of a field)
template<typename oT>
inline
const subview_field<oT>
field<oT>::slice(const uword slice_num) const
  {
  arma_debug_sigprint();
  
  arma_conform_check_bounds( (slice_num >= n_slices), "field::slice(): out of bounds" );
  
  return subview_field<oT>(*this, 0, 0, slice_num, n_rows, n_cols, 1);
  }



//! creation of subview_field (subfield comprised of specified rows)
template<typename oT>
inline
subview_field<oT>
field<oT>::rows(const uword in_row1, const uword in_row2)
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::rows(): field must be 2D" );
  
  arma_conform_check_bounds
    (
    ( (in_row1 > in_row2) || (in_row2 >= n_rows) ),
    "field::rows(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_rows = in_row2 - in_row1 + 1;
  
  return subview_field<oT>(*this, in_row1, 0, sub_n_rows, n_cols);
  }



//! creation of subview_field (subfield comprised of specified rows)
template<typename oT>
inline
const subview_field<oT>
field<oT>::rows(const uword in_row1, const uword in_row2) const
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::rows(): field must be 2D" );
  
  arma_conform_check_bounds
    (
    ( (in_row1 > in_row2) || (in_row2 >= n_rows) ),
    "field::rows(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_rows = in_row2 - in_row1 + 1;
  
  return subview_field<oT>(*this, in_row1, 0, sub_n_rows, n_cols);
  }



//! creation of subview_field (subfield comprised of specified columns)
template<typename oT>
inline
subview_field<oT>
field<oT>::cols(const uword in_col1, const uword in_col2)
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::cols(): field must be 2D" );
  
  arma_conform_check_bounds
    (
    ( (in_col1 > in_col2) || (in_col2 >= n_cols) ),
    "field::cols(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_cols = in_col2 - in_col1 + 1;
  
  return subview_field<oT>(*this, 0, in_col1, n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield comprised of specified columns)
template<typename oT>
inline
const subview_field<oT>
field<oT>::cols(const uword in_col1, const uword in_col2) const
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::cols(): field must be 2D" );
  
  arma_conform_check_bounds
    (
    ( (in_col1 > in_col2) || (in_col2 >= n_cols) ),
    "field::cols(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_cols = in_col2 - in_col1 + 1;
  
  return subview_field<oT>(*this, 0, in_col1, n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield comprised of specified slices)
template<typename oT>
inline
subview_field<oT>
field<oT>::slices(const uword in_slice1, const uword in_slice2)
  {
  arma_debug_sigprint();
  
  arma_conform_check_bounds
    (
    ( (in_slice1 > in_slice2) || (in_slice2 >= n_slices) ),
    "field::slices(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_slices = in_slice2 - in_slice1 + 1;
  
  return subview_field<oT>(*this, 0, 0, in_slice1, n_rows, n_cols, sub_n_slices);
  }



//! creation of subview_field (subfield comprised of specified slices)
template<typename oT>
inline
const subview_field<oT>
field<oT>::slices(const uword in_slice1, const uword in_slice2) const
  {
  arma_debug_sigprint();
  
  arma_conform_check_bounds
    (
    ( (in_slice1 > in_slice2) || (in_slice2 >= n_slices) ),
    "field::slices(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_slices = in_slice2 - in_slice1 + 1;
  
  return subview_field<oT>(*this, 0, 0, in_slice1, n_rows, n_cols, sub_n_slices);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
subview_field<oT>
field<oT>::subfield(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2)
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::subfield(): field must be 2D" );
  
  arma_conform_check_bounds
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_rows = in_row2 - in_row1 + 1;
  const uword sub_n_cols = in_col2 - in_col1 + 1;
  
  return subview_field<oT>(*this, in_row1, in_col1, sub_n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
const subview_field<oT>
field<oT>::subfield(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::subfield(): field must be 2D" );
  
  arma_conform_check_bounds
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_rows = in_row2 - in_row1 + 1;
  const uword sub_n_cols = in_col2 - in_col1 + 1;
  
  return subview_field<oT>(*this, in_row1, in_col1, sub_n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
subview_field<oT>
field<oT>::subfield(const uword in_row1, const uword in_col1, const uword in_slice1, const uword in_row2, const uword in_col2, const uword in_slice2)
  {
  arma_debug_sigprint();
  
  arma_conform_check_bounds
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_slice1 > in_slice2) || (in_row2 >= n_rows) || (in_col2 >= n_cols) || (in_slice2 >= n_slices),
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_rows   = in_row2   - in_row1   + 1;
  const uword sub_n_cols   = in_col2   - in_col1   + 1;
  const uword sub_n_slices = in_slice2 - in_slice1 + 1;
  
  return subview_field<oT>(*this, in_row1, in_col1, in_slice1, sub_n_rows, sub_n_cols, sub_n_slices);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
const subview_field<oT>
field<oT>::subfield(const uword in_row1, const uword in_col1, const uword in_slice1, const uword in_row2, const uword in_col2, const uword in_slice2) const
  {
  arma_debug_sigprint();
  
  arma_conform_check_bounds
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_slice1 > in_slice2) || (in_row2 >= n_rows) || (in_col2 >= n_cols) || (in_slice2 >= n_slices),
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  const uword sub_n_rows   = in_row2   - in_row1   + 1;
  const uword sub_n_cols   = in_col2   - in_col1   + 1;
  const uword sub_n_slices = in_slice2 - in_slice1 + 1;
  
  return subview_field<oT>(*this, in_row1, in_col1, in_slice1, sub_n_rows, sub_n_cols, sub_n_slices);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
subview_field<oT>
field<oT>::subfield(const uword in_row1, const uword in_col1, const SizeMat& s)
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::subfield(): field must be 2D" );
  
  const uword l_n_rows = n_rows;
  const uword l_n_cols = n_cols;
  
  const uword s_n_rows = s.n_rows;
  const uword s_n_cols = s.n_cols;
  
  arma_conform_check_bounds
    (
    ((in_row1 >= l_n_rows) || (in_col1 >= l_n_cols) || ((in_row1 + s_n_rows) > l_n_rows) || ((in_col1 + s_n_cols) > l_n_cols)),
    "field::subfield(): indices or size out of bounds"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, s_n_rows, s_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
const subview_field<oT>
field<oT>::subfield(const uword in_row1, const uword in_col1, const SizeMat& s) const
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::subfield(): field must be 2D" );
  
  const uword l_n_rows = n_rows;
  const uword l_n_cols = n_cols;
  
  const uword s_n_rows = s.n_rows;
  const uword s_n_cols = s.n_cols;
  
  arma_conform_check_bounds
    (
    ((in_row1 >= l_n_rows) || (in_col1 >= l_n_cols) || ((in_row1 + s_n_rows) > l_n_rows) || ((in_col1 + s_n_cols) > l_n_cols)),
    "field::subfield(): indices or size out of bounds"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, s_n_rows, s_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
subview_field<oT>
field<oT>::subfield(const uword in_row1, const uword in_col1, const uword in_slice1, const SizeCube& s)
  {
  arma_debug_sigprint();
  
  const uword l_n_rows   = n_rows;
  const uword l_n_cols   = n_cols;
  const uword l_n_slices = n_slices;
  
  const uword s_n_rows     = s.n_rows;
  const uword s_n_cols     = s.n_cols;
  const uword sub_n_slices = s.n_slices;
  
  arma_conform_check_bounds
    (
    ((in_row1 >= l_n_rows) || (in_col1 >= l_n_cols) || (in_slice1 >= l_n_slices) || ((in_row1 + s_n_rows) > l_n_rows) || ((in_col1 + s_n_cols) > l_n_cols) || ((in_slice1 + sub_n_slices) > l_n_slices)),
    "field::subfield(): indices or size out of bounds"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, in_slice1, s_n_rows, s_n_cols, sub_n_slices);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
const subview_field<oT>
field<oT>::subfield(const uword in_row1, const uword in_col1, const uword in_slice1, const SizeCube& s) const
  {
  arma_debug_sigprint();
  
  const uword l_n_rows = n_rows;
  const uword l_n_cols = n_cols;
  const uword l_n_slices = n_slices;
  
  const uword s_n_rows = s.n_rows;
  const uword s_n_cols = s.n_cols;
  const uword sub_n_slices = s.n_slices;
  
  arma_conform_check_bounds
    (
    ((in_row1 >= l_n_rows) || (in_col1 >= l_n_cols) || (in_slice1 >= l_n_slices) || ((in_row1 + s_n_rows) > l_n_rows) || ((in_col1 + s_n_cols) > l_n_cols) || ((in_slice1 + sub_n_slices) > l_n_slices)),
    "field::subfield(): indices or size out of bounds"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, in_slice1, s_n_rows, s_n_cols, sub_n_slices);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
subview_field<oT>
field<oT>::subfield(const span& row_span, const span& col_span)
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::subfield(): field must be 2D" );
  
  const bool row_all = row_span.whole;
  const bool col_all = col_span.whole;
  
  const uword local_n_rows = n_rows;
  const uword local_n_cols = n_cols;
  
  const uword in_row1    = row_all ? 0            : row_span.a;
  const uword in_row2    =                          row_span.b;
  const uword sub_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const uword in_col1    = col_all ? 0            : col_span.a;
  const uword in_col2    =                          col_span.b;
  const uword sub_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  arma_conform_check_bounds
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, sub_n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
const subview_field<oT>
field<oT>::subfield(const span& row_span, const span& col_span) const
  {
  arma_debug_sigprint();
  
  arma_conform_check( (n_slices >= 2), "field::subfield(): field must be 2D" );
  
  const bool row_all = row_span.whole;
  const bool col_all = col_span.whole;
  
  const uword local_n_rows = n_rows;
  const uword local_n_cols = n_cols;
  
  const uword in_row1    = row_all ? 0            : row_span.a;
  const uword in_row2    =                          row_span.b;
  const uword sub_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const uword in_col1    = col_all ? 0            : col_span.a;
  const uword in_col2    =                          col_span.b;
  const uword sub_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  arma_conform_check_bounds
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, sub_n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
subview_field<oT>
field<oT>::subfield(const span& row_span, const span& col_span, const span& slice_span)
  {
  arma_debug_sigprint();
  
  const bool row_all   = row_span.whole;
  const bool col_all   = col_span.whole;
  const bool slice_all = slice_span.whole;
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  const uword in_row1    = row_all ? 0            : row_span.a;
  const uword in_row2    =                          row_span.b;
  const uword sub_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const uword in_col1    = col_all ? 0            : col_span.a;
  const uword in_col2    =                          col_span.b;
  const uword sub_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  const uword in_slice1    = slice_all ? 0              : slice_span.a;
  const uword in_slice2    =                              slice_span.b;
  const uword sub_n_slices = slice_all ? local_n_slices : in_slice2 - in_slice1 + 1;

  arma_conform_check_bounds
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ||
    ( slice_all ? false : ((in_slice1 > in_slice2) || (in_slice2 >= local_n_slices)) )
    ,
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, in_slice1, sub_n_rows, sub_n_cols, sub_n_slices);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
const subview_field<oT>
field<oT>::subfield(const span& row_span, const span& col_span, const span& slice_span) const
  {
  arma_debug_sigprint();
  
  const bool row_all   = row_span.whole;
  const bool col_all   = col_span.whole;
  const bool slice_all = slice_span.whole;
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  const uword in_row1    = row_all ? 0            : row_span.a;
  const uword in_row2    =                          row_span.b;
  const uword sub_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const uword in_col1    = col_all ? 0            : col_span.a;
  const uword in_col2    =                          col_span.b;
  const uword sub_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  const uword in_slice1    = slice_all ? 0              : slice_span.a;
  const uword in_slice2    =                              slice_span.b;
  const uword sub_n_slices = slice_all ? local_n_slices : in_slice2 - in_slice1 + 1;

  arma_conform_check_bounds
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ||
    ( slice_all ? false : ((in_slice1 > in_slice2) || (in_slice2 >= local_n_slices)) )
    ,
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, in_slice1, sub_n_rows, sub_n_cols, sub_n_slices);
  }



template<typename oT>
inline
subview_field<oT>
field<oT>::operator()(const span& row_span, const span& col_span)
  {
  arma_debug_sigprint();
  
  return (*this).subfield(row_span, col_span);
  }



template<typename oT>
inline
const subview_field<oT>
field<oT>::operator()(const span& row_span, const span& col_span) const
  {
  arma_debug_sigprint();
  
  return (*this).subfield(row_span, col_span);
  }



template<typename oT>
inline
subview_field<oT>
field<oT>::operator()(const span& row_span, const span& col_span, const span& slice_span)
  {
  arma_debug_sigprint();
  
  return (*this).subfield(row_span, col_span, slice_span);
  }



template<typename oT>
inline
const subview_field<oT>
field<oT>::operator()(const span& row_span, const span& col_span, const span& slice_span) const
  {
  arma_debug_sigprint();
  
  return (*this).subfield(row_span, col_span, slice_span);
  }



template<typename oT>
inline
subview_field<oT>
field<oT>::operator()(const uword in_row1, const uword in_col1, const SizeMat& s)
  {
  arma_debug_sigprint();
  
  return (*this).subfield(in_row1, in_col1, s);
  }



template<typename oT>
inline
const subview_field<oT>
field<oT>::operator()(const uword in_row1, const uword in_col1, const SizeMat& s) const
  {
  arma_debug_sigprint();
  
  return (*this).subfield(in_row1, in_col1, s);
  }



template<typename oT>
inline
subview_field<oT>
field<oT>::operator()(const uword in_row1, const uword in_col1, const uword in_slice1, const SizeCube& s)
  {
  arma_debug_sigprint();
  
  return (*this).subfield(in_row1, in_col1, in_slice1, s);
  }



template<typename oT>
inline
const subview_field<oT>
field<oT>::operator()(const uword in_row1, const uword in_col1, const uword in_slice1, const SizeCube& s) const
  {
  arma_debug_sigprint();
  
  return (*this).subfield(in_row1, in_col1, in_slice1, s);
  }



//! print contents of the field (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the field class preserves the stream's flags
//! but the associated operator<< function for type oT 
//! may still modify the stream's parameters.
//! NOTE: this function assumes that type oT can be printed,
//! ie. the function "std::ostream& operator<< (std::ostream&, const oT&)"
//! has been defined.

template<typename oT>
inline
void
field<oT>::print(const std::string extra_text) const
  {
  arma_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = get_cout_stream().width();
    
    get_cout_stream() << extra_text << '\n';
    
    get_cout_stream().width(orig_width);
    }
  
  arma_ostream::print(get_cout_stream(), *this);
  }



//! print contents of the field to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the field class preserves the stream's flags
//! but the associated operator<< function for type oT 
//! may still modify the stream's parameters.
//! NOTE: this function assumes that type oT can be printed,
//! ie. the function "std::ostream& operator<< (std::ostream&, const oT&)"
//! has been defined.

template<typename oT>
inline
void
field<oT>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
  
    user_stream.width(orig_width);
    }
  
  arma_ostream::print(user_stream, *this);
  }



//! apply a lambda function to each object
template<typename oT>
inline
field<oT>&
field<oT>::for_each(const std::function< void(oT&) >& F)
  {
  arma_debug_sigprint();
  
  for(uword i=0; i < n_elem; ++i)  { F(operator[](i)); }
  
  return *this;
  }



template<typename oT>
inline
const field<oT>&
field<oT>::for_each(const std::function< void(const oT&) >& F) const
  {
  arma_debug_sigprint();
  
  for(uword i=0; i < n_elem; ++i)  { F(operator[](i)); }
  
  return *this;
  }



//! fill the field with an object
template<typename oT>
inline
field<oT>&
field<oT>::fill(const oT& x)
  {
  arma_debug_sigprint();
  
  field<oT>& t = *this;
  
  for(uword i=0; i<n_elem; ++i)  { t[i] = x; }
  
  return *this;
  }



//! reset the field to an empty state (ie. the field will have no objects)
template<typename oT>
inline
void
field<oT>::reset()
  {
  arma_debug_sigprint();
  
  init(0,0,0);
  }



//! reset each object
template<typename oT>
inline
void
field<oT>::reset_objects()
  {
  arma_debug_sigprint();
  
  field_aux::reset_objects(*this);
  }



//! returns true if the field has no objects
template<typename oT>
arma_inline
bool
field<oT>::is_empty() const
  {
  return (n_elem == 0);
  }



//! returns true if the given index is currently in range
template<typename oT>
arma_inline
bool
field<oT>::in_range(const uword i) const
  {
  return (i < n_elem);
  }



//! returns true if the given start and end indices are currently in range
template<typename oT>
arma_inline
bool
field<oT>::in_range(const span& x) const
  {
  arma_debug_sigprint();
  
  if(x.whole)
    {
    return true;
    }
  else
    {
    const uword a = x.a;
    const uword b = x.b;
    
    return ( (a <= b) && (b < n_elem) );
    }
  }



//! returns true if the given location is currently in range
template<typename oT>
arma_inline
bool
field<oT>::in_range(const uword in_row, const uword in_col) const
  {
  return ( (in_row < n_rows) && (in_col < n_cols) );
  }



template<typename oT>
arma_inline
bool
field<oT>::in_range(const span& row_span, const uword in_col) const
  {
  arma_debug_sigprint();
  
  if(row_span.whole)
    {
    return (in_col < n_cols);
    }
  else
    {
    const uword in_row1 = row_span.a;
    const uword in_row2 = row_span.b;
    
    return ( (in_row1 <= in_row2) && (in_row2 < n_rows) && (in_col < n_cols) );
    }
  }



template<typename oT>
arma_inline
bool
field<oT>::in_range(const uword in_row, const span& col_span) const
  {
  arma_debug_sigprint();
  
  if(col_span.whole)
    {
    return (in_row < n_rows);
    }
  else
    {
    const uword in_col1 = col_span.a;
    const uword in_col2 = col_span.b;
  
    return ( (in_row < n_rows) && (in_col1 <= in_col2) && (in_col2 < n_cols) );
    }
  }



template<typename oT>
arma_inline
bool
field<oT>::in_range(const span& row_span, const span& col_span) const
  {
  arma_debug_sigprint();
  
  const uword in_row1 = row_span.a;
  const uword in_row2 = row_span.b;
  
  const uword in_col1 = col_span.a;
  const uword in_col2 = col_span.b;
  
  const bool rows_ok = row_span.whole ? true : ( (in_row1 <= in_row2) && (in_row2 < n_rows) );
  const bool cols_ok = col_span.whole ? true : ( (in_col1 <= in_col2) && (in_col2 < n_cols) );
  
  return ( rows_ok && cols_ok );
  }



template<typename oT>
arma_inline
bool
field<oT>::in_range(const uword in_row, const uword in_col, const SizeMat& s) const
  {
  const uword l_n_rows = n_rows;
  const uword l_n_cols = n_cols;
  
  if( (in_row >= l_n_rows) || (in_col >= l_n_cols) || ((in_row + s.n_rows) > l_n_rows) || ((in_col + s.n_cols) > l_n_cols) )
    {
    return false;
    }
  else
    {
    return true;
    }
  }



template<typename oT>
arma_inline
bool
field<oT>::in_range(const uword in_row, const uword in_col, const uword in_slice) const
  {
  return ( (in_row < n_rows) && (in_col < n_cols) && (in_slice < n_slices) );
  }



template<typename oT>
arma_inline
bool
field<oT>::in_range(const span& row_span, const span& col_span, const span& slice_span) const
  {
  arma_debug_sigprint();
  
  const uword in_row1 = row_span.a;
  const uword in_row2 = row_span.b;

  const uword in_col1 = col_span.a;
  const uword in_col2 = col_span.b;
  
  const uword in_slice1 = slice_span.a;
  const uword in_slice2 = slice_span.b;
  
  const bool   rows_ok =   row_span.whole ? true : ( (in_row1   <= in_row2  ) && (in_row2   < n_rows  ) );
  const bool   cols_ok =   col_span.whole ? true : ( (in_col1   <= in_col2  ) && (in_col2   < n_cols  ) );
  const bool slices_ok = slice_span.whole ? true : ( (in_slice1 <= in_slice2) && (in_slice2 < n_slices) );
  
  return ( rows_ok && cols_ok && slices_ok );
  }



template<typename oT>
arma_inline
bool
field<oT>::in_range(const uword in_row, const uword in_col, const uword in_slice, const SizeCube& s) const
  {
  const uword l_n_rows = n_rows;
  const uword l_n_cols = n_cols;
  const uword l_n_slices = n_slices;
  
  if( (in_row >= l_n_rows) || (in_col >= l_n_cols) || (in_slice >= l_n_slices) || ((in_row + s.n_rows) > l_n_rows) || ((in_col + s.n_cols) > l_n_cols) || ((in_slice + s.n_slices) > l_n_slices) )
    {
    return false;
    }
  else
    {
    return true;
    }
  }



template<typename oT>
inline
bool
field<oT>::save(const std::string name, const file_type type) const
  {
  arma_debug_sigprint();
  
  std::string err_msg;
  
  const bool save_okay = field_aux::save(*this, name, type, err_msg);
  
  if(save_okay == false)
    {
    if(err_msg.length() > 0)
      {
      arma_warn(3, "field::save(): ", err_msg, "; file: ", name);
      }
    else
      {
      arma_warn(3, "field::save(): write failed; file: ", name);
      }
    }
  
  return save_okay;
  }



template<typename oT>
inline
bool
field<oT>::save(std::ostream& os, const file_type type) const
  {
  arma_debug_sigprint();
  
  std::string err_msg;
  
  const bool save_okay = field_aux::save(*this, os, type, err_msg);
  
  if(save_okay == false)
    {
    if(err_msg.length() > 0)
      {
      arma_warn(3, "field::save(): ", err_msg);
      }
    else
      {
      arma_warn(3, "field::save(): stream write failed");
      }
    }
  
  return save_okay;
  }



template<typename oT>
inline
bool
field<oT>::load(const std::string name, const file_type type)
  {
  arma_debug_sigprint();
  
  std::string err_msg;
  
  const bool load_okay = field_aux::load(*this, name, type, err_msg);
  
  if(load_okay == false)
    {
    if(err_msg.length() > 0)
      {
      arma_warn(3, "field::load(): ", err_msg, "; file: ", name);
      }
    else
      {
      arma_warn(3, "field::load(): read failed; file: ", name);
      }
    }
  
  if(load_okay == false)  { (*this).reset(); }
  
  return load_okay;
  }



template<typename oT>
inline
bool
field<oT>::load(std::istream& is, const file_type type)
  {
  arma_debug_sigprint();
  
  std::string err_msg;
  const bool load_okay = field_aux::load(*this, is, type, err_msg);
  
  if(load_okay == false)
    {
    if(err_msg.length() > 0)
      {
      arma_warn(3, "field::load(): ", err_msg);
      }
    else
      {
      arma_warn(3, "field::load(): stream read failed");
      }
    }
  
  if(load_okay == false)  { (*this).reset(); }
  
  return load_okay;
  }



template<typename oT>
inline
bool
field<oT>::quiet_save(const std::string name, const file_type type) const
  {
  arma_debug_sigprint();
  
  return (*this).save(name, type);
  }



template<typename oT>
inline
bool
field<oT>::quiet_save(std::ostream& os, const file_type type) const
  {
  arma_debug_sigprint();
  
  return (*this).save(os, type);
  }



template<typename oT>
inline
bool
field<oT>::quiet_load(const std::string name, const file_type type)
  {
  arma_debug_sigprint();
  
  return (*this).load(name, type);
  }



template<typename oT>
inline
bool
field<oT>::quiet_load(std::istream& is, const file_type type)
  {
  arma_debug_sigprint();
  
  return (*this).load(is, type);
  }



//! construct a field from a given field
template<typename oT>
inline
void
field<oT>::init(const field<oT>& x)
  {
  arma_debug_sigprint();
  
  if(this == &x)
    {
    arma_debug_print("field::init(): copy omitted");
    
    return;
    }
  
  field& t = (*this);
  
  t.init(x.n_rows, x.n_cols, x.n_slices);
  
  const uword t_n_elem = t.n_elem;
  
  for(uword i=0; i < t_n_elem; ++i)  { t.at(i) = x.at(i); }
  }



template<typename oT>
inline
void
field<oT>::init(const uword n_rows_in, const uword n_cols_in)
  {
  (*this).init(n_rows_in, n_cols_in, 1);
  }



template<typename oT>
inline
void
field<oT>::init(const uword n_rows_in, const uword n_cols_in, const uword n_slices_in)
  {
  arma_debug_sigprint( arma_str::format("n_rows_in: %u; n_cols_in: %u; n_slices_in: %u") % n_rows_in % n_cols_in % n_slices_in );
  
  #if defined(ARMA_64BIT_WORD)
    const char* error_message = "field::init(): requested size is too large";
  #else
    const char* error_message = "field::init(): requested size is too large; suggest to enable ARMA_64BIT_WORD";
  #endif
  
  arma_conform_check
    (
      (
      ( (n_rows_in > 0x0FFF) || (n_cols_in > 0x0FFF) || (n_slices_in > 0xFF) )
        ? ( (double(n_rows_in) * double(n_cols_in) * double(n_slices_in)) > double(ARMA_MAX_UWORD) )
        : false
      ),
    error_message
    );
  
  const uword n_elem_new = n_rows_in * n_cols_in * n_slices_in;
  
  if(n_elem == n_elem_new)
    {
    // delete_objects();
    // create_objects();
    access::rw(n_rows)   = n_rows_in;
    access::rw(n_cols)   = n_cols_in;
    access::rw(n_slices) = n_slices_in;
    }
  else
    {
    delete_objects();
    
    if(n_elem > 0)  { delete [] mem; }
    
    mem = nullptr;
    
    if(n_elem_new > 0)
      {
      mem = new(std::nothrow) oT* [n_elem_new];
      
      arma_check_bad_alloc( (mem == nullptr), "field::init(): out of memory" );
      }
    
    access::rw(n_rows)   = n_rows_in;
    access::rw(n_cols)   = n_cols_in;
    access::rw(n_slices) = n_slices_in;
    access::rw(n_elem)   = n_elem_new;
    
    create_objects();
    }
  }



template<typename oT>
inline
void
field<oT>::delete_objects()
  {
  arma_debug_sigprint( arma_str::format("n_elem: %u") % n_elem );
  
  for(uword i=0; i<n_elem; ++i)
    {
    if(mem[i] != nullptr)  { delete mem[i]; mem[i] = nullptr; }
    }
  }



template<typename oT>
inline
void
field<oT>::create_objects()
  {
  arma_debug_sigprint( arma_str::format("n_elem: %u") % n_elem );
  
  for(uword i=0; i<n_elem; ++i)  { mem[i] = new oT; }
  }



template<typename oT>
inline
field<oT>::iterator::iterator(field<oT>& in_M, const bool at_end)
  : M(in_M)
  , i( (at_end == false) ? 0 : in_M.n_elem )
  {
  arma_debug_sigprint();
  }



template<typename oT>
inline
oT&
field<oT>::iterator::operator*()
  {
  return M[i];
  }



template<typename oT>
inline
typename field<oT>::iterator&
field<oT>::iterator::operator++()
  {
  ++i;
  
  return *this;
  }



template<typename oT>
inline
void
field<oT>::iterator::operator++(int)
  {
  operator++();
  }



template<typename oT>
inline
typename field<oT>::iterator&
field<oT>::iterator::operator--()
  {
  if(i > 0)  { --i; }
  
  return *this;
  }



template<typename oT>
inline
void
field<oT>::iterator::operator--(int)
  {
  operator--();
  }



template<typename oT>
inline
bool
field<oT>::iterator::operator!=(const typename field<oT>::iterator& X) const
  {
  return (i != X.i);
  }



template<typename oT>
inline
bool
field<oT>::iterator::operator==(const typename field<oT>::iterator& X) const
  {
  return (i == X.i);
  }



template<typename oT>
inline
field<oT>::const_iterator::const_iterator(const field<oT>& in_M, const bool at_end)
  : M(in_M)
  , i( (at_end == false) ? 0 : in_M.n_elem )
  {
  arma_debug_sigprint();
  }



template<typename oT>
inline
field<oT>::const_iterator::const_iterator(const typename field<oT>::iterator& X)
  : M(X.M)
  , i(X.i)
  {
  arma_debug_sigprint();
  }



template<typename oT>
inline
const oT&
field<oT>::const_iterator::operator*() const
  {
  return M[i];
  }



template<typename oT>
inline
typename field<oT>::const_iterator&
field<oT>::const_iterator::operator++()
  {
  ++i;
  
  return *this;
  }



template<typename oT>
inline
void
field<oT>::const_iterator::operator++(int)
  {
  operator++();
  }



template<typename oT>
inline
typename field<oT>::const_iterator&
field<oT>::const_iterator::operator--()
  {
  if(i > 0)  { --i; }
  
  return *this;
  }



template<typename oT>
inline
void
field<oT>::const_iterator::operator--(int)
  {
  operator--();
  }



template<typename oT>
inline
bool
field<oT>::const_iterator::operator!=(const typename field<oT>::const_iterator& X) const
  {
  return (i != X.i);
  }



template<typename oT>
inline
bool
field<oT>::const_iterator::operator==(const typename field<oT>::const_iterator& X) const
  {
  return (i == X.i);
  }



template<typename oT>
inline
typename field<oT>::iterator
field<oT>::begin()
  {
  arma_debug_sigprint();
  
  return field<oT>::iterator(*this);
  }



template<typename oT>
inline
typename field<oT>::const_iterator
field<oT>::begin() const
  {
  arma_debug_sigprint();
  
  return field<oT>::const_iterator(*this);
  }



template<typename oT>
inline
typename field<oT>::const_iterator
field<oT>::cbegin() const
  {
  arma_debug_sigprint();
  
  return field<oT>::const_iterator(*this);
  }



template<typename oT>
inline
typename field<oT>::iterator
field<oT>::end()
  {
  arma_debug_sigprint();
  
  return field<oT>::iterator(*this, true);
  }



template<typename oT>
inline
typename field<oT>::const_iterator
field<oT>::end() const
  {
  arma_debug_sigprint();
  
  return field<oT>::const_iterator(*this, true);
  }
  


template<typename oT>
inline
typename field<oT>::const_iterator
field<oT>::cend() const
  {
  arma_debug_sigprint();
  
  return field<oT>::const_iterator(*this, true);
  }



template<typename oT>
inline
void
field<oT>::clear()
  {
  reset();
  }



template<typename oT>
inline
bool
field<oT>::empty() const
  {
  return (n_elem == 0);
  }



template<typename oT>
inline
uword
field<oT>::size() const
  {
  return n_elem;
  }



//
//
//



template<typename oT>
inline
void
field_aux::reset_objects(field<oT>& x)
  {
  arma_debug_sigprint();
  
  x.delete_objects();
  x.create_objects();
  }



template<typename eT>
inline
void
field_aux::reset_objects(field< Mat<eT> >& x)
  {
  arma_debug_sigprint();
  
  for(uword i=0; i < x.n_elem; ++i)  { (*(x.mem[i])).reset(); }
  }



template<typename eT>
inline
void
field_aux::reset_objects(field< Col<eT> >& x)
  {
  arma_debug_sigprint();
  
  for(uword i=0; i < x.n_elem; ++i)  { (*(x.mem[i])).reset(); }
  }



template<typename eT>
inline
void
field_aux::reset_objects(field< Row<eT> >& x)
  {
  arma_debug_sigprint();
  
  for(uword i=0; i < x.n_elem; ++i)  { (*(x.mem[i])).reset(); }
  }



template<typename eT>
inline
void
field_aux::reset_objects(field< Cube<eT> >& x)
  {
  arma_debug_sigprint();
  
  for(uword i=0; i < x.n_elem; ++i)  { (*(x.mem[i])).reset(); }
  }



inline
void
field_aux::reset_objects(field< std::string >& x)
  {
  arma_debug_sigprint();
  
  for(uword i=0; i < x.n_elem; ++i)  { (*(x.mem[i])).clear(); }
  }



//
//
//



template<typename oT>
inline
bool
field_aux::save(const field<oT>&, const std::string&, const file_type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  err_msg = "saving/loading this type of field is currently not supported";
  
  return false;
  }



template<typename oT>
inline
bool
field_aux::save(const field<oT>&, std::ostream&, const file_type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  err_msg = "saving/loading this type of field is currently not supported";
  
  return false;
  }



template<typename oT>
inline
bool
field_aux::load(field<oT>&, const std::string&, const file_type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  err_msg = "saving/loading this type of field is currently not supported";
  
  return false;
  }



template<typename oT>
inline
bool
field_aux::load(field<oT>&, std::istream&, const file_type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  err_msg = "saving/loading this type of field is currently not supported";
  
  return false;
  }



template<typename eT>
inline
bool
field_aux::save(const field< Mat<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, name);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, name);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Mat<eT> >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, os);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, os);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Mat<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      return diskio::load_auto_detect(x, name, err_msg);
      break;
    
    case arma_binary:
      return diskio::load_arma_binary(x, name, err_msg);
      break;
      
    case ppm_binary:
      return diskio::load_ppm_binary(x, name, err_msg);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Mat<eT> >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      return diskio::load_auto_detect(x, is, err_msg);
      break;
    
    case arma_binary:
      return diskio::load_arma_binary(x, is, err_msg);
      break;
      
    case ppm_binary:
      return diskio::load_ppm_binary(x, is, err_msg);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Col<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, name);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, name);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Col<eT> >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, os);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, os);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Col<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      return diskio::load_auto_detect(x, name, err_msg);
      break;
    
    case arma_binary:
      return diskio::load_arma_binary(x, name, err_msg);
      break;
      
    case ppm_binary:
      return diskio::load_ppm_binary(x, name, err_msg);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Col<eT> >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      return diskio::load_auto_detect(x, is, err_msg);
      break;
    
    case arma_binary:
      return diskio::load_arma_binary(x, is, err_msg);
      break;
      
    case ppm_binary:
      return diskio::load_ppm_binary(x, is, err_msg);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Row<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, name);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, name);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Row<eT> >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, os);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, os);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Row<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      return diskio::load_auto_detect(x, name, err_msg);
      break;
    
    case arma_binary:
      return diskio::load_arma_binary(x, name, err_msg);
      break;
      
    case ppm_binary:
      return diskio::load_ppm_binary(x, name, err_msg);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Row<eT> >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      return diskio::load_auto_detect(x, is, err_msg);
      break;
    
    case arma_binary:
      return diskio::load_arma_binary(x, is, err_msg);
      break;
      
    case ppm_binary:
      return diskio::load_ppm_binary(x, is, err_msg);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Cube<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, name);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Cube<eT> >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, os);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Cube<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
    case arma_binary:
      return diskio::load_arma_binary(x, name, err_msg);
      break;
    
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Cube<eT> >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
    case arma_binary:
      return diskio::load_arma_binary(x, is, err_msg);
      break;
      
    default:
      err_msg = "unsupported type";
      return false;
    }
  }



inline
bool
field_aux::save(const field< std::string >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  arma_ignore(type);
  
  err_msg.clear();
  
  return diskio::save_std_string(x, name);
  }



inline
bool
field_aux::save(const field< std::string >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  arma_ignore(type);
  
  err_msg.clear();
  
  return diskio::save_std_string(x, os);
  }



inline
bool
field_aux::load(field< std::string >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  arma_ignore(type);
  
  return diskio::load_std_string(x, name, err_msg);
  }



inline
bool
field_aux::load(field< std::string >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_debug_sigprint();
  
  arma_ignore(type);
  
  return diskio::load_std_string(x, is, err_msg);
  }



#if defined(ARMA_EXTRA_FIELD_MEAT)
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_FIELD_MEAT)
#endif



//! @}
