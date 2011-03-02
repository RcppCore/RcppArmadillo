// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// Copyright (C) 2009-2010 Ian Cullinan
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup field
//! @{


template<typename oT>
inline
field<oT>::~field()
  {
  arma_extra_debug_sigprint_this(this);
  
  delete_objects();
  
  if(n_elem > sizeof(mem_local)/sizeof(oT*) )
    {
    delete [] mem;
    }
  
  if(arma_config::debug == true)
    {
    // try to expose buggy user code that accesses deleted objects
    access::rw(n_rows) = 0;
    access::rw(n_cols) = 0;
    access::rw(n_elem) = 0;
    mem = 0;
    }
  }



template<typename oT>
inline
field<oT>::field()
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , mem(0)
  {
  arma_extra_debug_sigprint_this(this);
  }



//! construct a field from a given field
template<typename oT>
inline
field<oT>::field(const field& x)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , mem(0)
  {
  arma_extra_debug_sigprint(arma_boost::format("this = %x   x = %x") % this % &x);
  
  init(x);
  }



//! construct a field from a given field
template<typename oT>
inline
const field<oT>&
field<oT>::operator=(const field& x)
  {
  arma_extra_debug_sigprint();
  
  init(x);
  return *this;
  }



//! construct a field from subview_field (e.g. construct a field from a delayed subfield operation)
template<typename oT>
inline
field<oT>::field(const subview_field<oT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , mem(0)
  {
  arma_extra_debug_sigprint_this(this);
  
  this->operator=(X);
  }



//! construct a field from subview_field (e.g. construct a field from a delayed subfield operation)
template<typename oT>
inline
const field<oT>&
field<oT>::operator=(const subview_field<oT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_field<oT>::extract(*this, X);
  return *this;
  }



//! construct the field with the specified number of elements,
//! assuming a column-major layout
template<typename oT>
inline
field<oT>::field(const u32 n_elem_in)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , mem(0)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(n_elem_in, 1);
  }



//! construct the field with the specified dimensions
template<typename oT>
inline
field<oT>::field(const u32 n_rows_in, const u32 n_cols_in)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , mem(0)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(n_rows_in, n_cols_in);
  }



//! change the field to have the specified number of elements,
//! assuming a column-major layout (data is not preserved)
template<typename oT>
inline
void
field<oT>::set_size(const u32 n_elem_in)
  {
  arma_extra_debug_sigprint(arma_boost::format("n_elem_in = %d") % n_elem_in);
  
  init(n_elem_in, 1);
  }



//! change the field to have the specified dimensions (data is not preserved)
template<typename oT>
inline
void
field<oT>::set_size(const u32 n_rows_in, const u32 n_cols_in)
  {
  arma_extra_debug_sigprint(arma_boost::format("n_rows_in = %d, n_cols_in = %d") % n_rows_in % n_cols_in);
  
  init(n_rows_in, n_cols_in);
  }



//! change the field to have the specified dimensions (data is not preserved)
template<typename oT>
template<typename oT2>
inline
void
field<oT>::copy_size(const field<oT2>& x)
  {
  arma_extra_debug_sigprint();
  
  init(x.n_rows, x.n_cols);
  }



//! linear element accessor (treats the field as a vector); no bounds check
template<typename oT>
arma_inline
oT&
field<oT>::operator[] (const u32 i)
  {
  return (*mem[i]);
  }
  
  
  
//! linear element accessor (treats the field as a vector); no bounds check
template<typename oT>
arma_inline
const oT&
field<oT>::operator[] (const u32 i) const
  {
  return (*mem[i]);
  }



//! linear element accessor (treats the field as a vector); no bounds check
template<typename oT>
arma_inline
oT&
field<oT>::at(const u32 i)
  {
  return (*mem[i]);
  }
  
  
  
//! linear element accessor (treats the field as a vector); no bounds check
template<typename oT>
arma_inline
const oT&
field<oT>::at(const u32 i) const
  {
  return (*mem[i]);
  }



//! linear element accessor (treats the field as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename oT>
arma_inline
oT&
field<oT>::operator() (const u32 i)
  {
  arma_debug_check( (i >= n_elem), "field::operator(): index out of bounds");
  return (*mem[i]);
  }
  
  
  
//! linear element accessor (treats the field as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename oT>
arma_inline
const oT&
field<oT>::operator() (const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "field::operator(): index out of bounds");
  return (*mem[i]);
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename oT>
arma_inline
oT&
field<oT>::operator() (const u32 in_row, const u32 in_col)
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "field::operator(): index out of bounds");
  return (*mem[in_row + in_col*n_rows]);
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename oT>
arma_inline
const oT&
field<oT>::operator() (const u32 in_row, const u32 in_col) const
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "field::operator(): index out of bounds");
  return (*mem[in_row + in_col*n_rows]);
  }



//! element accessor; no bounds check
template<typename oT>
arma_inline
oT&
field<oT>::at(const u32 in_row, const u32 in_col)
  {
  return (*mem[in_row + in_col*n_rows]);
  }



//! element accessor; no bounds check
template<typename oT>
arma_inline
const oT&
field<oT>::at(const u32 in_row, const u32 in_col) const
  {
  return (*mem[in_row + in_col*n_rows]);
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
field<oT>::operator<<(const injector_end_of_row& x)
  {
  return field_injector< field<oT> >(*this, x);
  }



//! creation of subview_field (row of a field)
template<typename oT>
inline
subview_field<oT>
field<oT>::row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= n_rows), "field::row(): row out of bounds" );
  
  return subview_field<oT>(*this, row_num, 0, 1, n_cols);
  }



//! creation of subview_field (row of a field)
template<typename oT>
inline
const subview_field<oT>
field<oT>::row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= n_rows), "field::row(): row out of bounds" );
  
  return subview_field<oT>(*this, row_num, 0, 1, n_cols);
  }



//! creation of subview_field (column of a field)
template<typename oT>
inline
subview_field<oT>
field<oT>::col(const u32 col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (col_num >= n_cols), "field::col(): out of bounds");
  
  return subview_field<oT>(*this, 0, col_num, n_rows, 1);
  }



//! creation of subview_field (column of a field)
template<typename oT>
inline
const subview_field<oT>
field<oT>::col(const u32 col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (col_num >= n_cols), "field::col(): out of bounds");
  
  return subview_field<oT>(*this, 0, col_num, n_rows, 1);
  }



//! creation of subview_field (subfield comprised of specified rows)
template<typename oT>
inline
subview_field<oT>
field<oT>::rows(const u32 in_row1, const u32 in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    ( (in_row1 > in_row2) || (in_row2 >= n_rows) ),
    "field::rows(): indicies out of bounds or incorrectly used"
    );
  
  const u32 sub_n_rows = in_row2 - in_row1 + 1;
  
  return subview_field<oT>(*this, in_row1, 0, sub_n_rows, n_cols);
  }



//! creation of subview_field (subfield comprised of specified rows)
template<typename oT>
inline
const subview_field<oT>
field<oT>::rows(const u32 in_row1, const u32 in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    ( (in_row1 > in_row2) || (in_row2 >= n_rows) ),
    "field::rows(): indicies out of bounds or incorrectly used"
    );
  
  const u32 sub_n_rows = in_row2 - in_row1 + 1;
  
  return subview_field<oT>(*this, in_row1, 0, sub_n_rows, n_cols);
  }



//! creation of subview_field (subfield comprised of specified columns)
template<typename oT>
inline
subview_field<oT>
field<oT>::cols(const u32 in_col1, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    ( (in_col1 > in_col2) || (in_col2 >= n_cols) ),
    "field::cols(): indicies out of bounds or incorrectly used"
    );
  
  const u32 sub_n_cols = in_col2 - in_col1 + 1;
  
  return subview_field<oT>(*this, 0, in_col1, n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield comprised of specified columns)
template<typename oT>
inline
const subview_field<oT>
field<oT>::cols(const u32 in_col1, const u32 in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    ( (in_col1 > in_col2) || (in_col2 >= n_cols) ),
    "field::cols(): indicies out of bounds or incorrectly used"
    );
  
  const u32 sub_n_cols = in_col2 - in_col1 + 1;
  
  return subview_field<oT>(*this, 0, in_col1, n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
subview_field<oT>
field<oT>::subfield(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  const u32 sub_n_rows = in_row2 - in_row1 + 1;
  const u32 sub_n_cols = in_col2 - in_col1 + 1;
  
  return subview_field<oT>(*this, in_row1, in_col1, sub_n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
const subview_field<oT>
field<oT>::subfield(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  const u32 sub_n_rows = in_row2 - in_row1 + 1;
  const u32 sub_n_cols = in_col2 - in_col1 + 1;
  
  return subview_field<oT>(*this, in_row1, in_col1, sub_n_rows, sub_n_cols);
  }



//! creation of subview_field (subfield with arbitrary dimensions)
template<typename oT>
inline
subview_field<oT>
field<oT>::subfield(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  const bool col_all = col_span.whole;
  
  const u32 local_n_rows = n_rows;
  const u32 local_n_cols = n_cols;
  
  const u32 in_row1    = row_all ? 0            : row_span.a;
  const u32 in_row2    =                          row_span.b;
  const u32 sub_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const u32 in_col1    = col_all ? 0            : col_span.a;
  const u32 in_col2    =                          col_span.b;
  const u32 sub_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  arma_debug_check
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
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  const bool col_all = col_span.whole;
  
  const u32 local_n_rows = n_rows;
  const u32 local_n_cols = n_cols;
  
  const u32 in_row1    = row_all ? 0            : row_span.a;
  const u32 in_row2    =                          row_span.b;
  const u32 sub_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  const u32 in_col1    = col_all ? 0            : col_span.a;
  const u32 in_col2    =                          col_span.b;
  const u32 sub_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  arma_debug_check
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "field::subfield(): indices out of bounds or incorrectly used"
    );
  
  return subview_field<oT>(*this, in_row1, in_col1, sub_n_rows, sub_n_cols);
  }



template<typename oT>
inline
subview_field<oT>
field<oT>::operator()(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return (*this).subfield(row_span, col_span);
  }



template<typename oT>
inline
const subview_field<oT>
field<oT>::operator()(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).subfield(row_span, col_span);
  }



//! print contents of the field (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the field class preserves the stream's flags
//! but the associated operator<< function for type oT 
//! may still modify the stream's parameters.
//! NOTE: this function assumes that type oT can be printed,
//! i.e. the function "std::ostream& operator<< (std::ostream&, const oT&)"
//! has been defined.

template<typename oT>
inline
void
field<oT>::print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = cout.width();
    
    cout << extra_text << '\n';
  
    cout.width(orig_width);
    }
  
  arma_ostream::print(cout, *this);
  }



//! print contents of the field to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the field class preserves the stream's flags
//! but the associated operator<< function for type oT 
//! may still modify the stream's parameters.
//! NOTE: this function assumes that type oT can be printed,
//! i.e. the function "std::ostream& operator<< (std::ostream&, const oT&)"
//! has been defined.

template<typename oT>
inline
void
field<oT>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
  
    user_stream.width(orig_width);
    }
  
  arma_ostream::print(user_stream, *this);
  }



//! fill the field with an object
template<typename oT>
inline
void
field<oT>::fill(const oT& x)
  {
  arma_extra_debug_sigprint();
  
  field<oT>& t = *this;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    t[i] = x;
    }
  }



//! reset the field to an empty state (i.e. the field will have no objects)
template<typename oT>
inline
void
field<oT>::reset()
  {
  arma_extra_debug_sigprint();
  
  init(0,0);
  }



//! reset each object
template<typename oT>
inline
void
field<oT>::reset_objects()
  {
  arma_extra_debug_sigprint();
  
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
arma_warn_unused
bool
field<oT>::in_range(const u32 i) const
  {
  return (i < n_elem);
  }



//! returns true if the given start and end indices are currently in range
template<typename oT>
arma_inline
arma_warn_unused
bool
field<oT>::in_range(const span& x) const
  {
  arma_extra_debug_sigprint();
  
  if(x.whole == true)
    {
    return true;
    }
  else
    {
    const u32 a = x.a;
    const u32 b = x.b;
    
    return ( (a <= b) && (b < n_elem) );
    }
  }



//! returns true if the given location is currently in range
template<typename oT>
arma_inline
arma_warn_unused
bool
field<oT>::in_range(const u32 in_row, const u32 in_col) const
  {
  return ( (in_row < n_rows) && (in_col < n_cols) );
  }



template<typename oT>
arma_inline
arma_warn_unused
bool
field<oT>::in_range(const span& row_span, const u32 in_col) const
  {
  arma_extra_debug_sigprint();
  
  if(row_span.whole == true)
    {
    return (in_col < n_cols);
    }
  else
    {
    const u32 in_row1 = row_span.a;
    const u32 in_row2 = row_span.b;
    
    return ( (in_row1 <= in_row2) && (in_row2 < n_rows) && (in_col < n_cols) );
    }
  }



template<typename oT>
arma_inline
arma_warn_unused
bool
field<oT>::in_range(const u32 in_row, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  if(col_span.whole == true)
    {
    return (in_row < n_rows);
    }
  else
    {
    const u32 in_col1 = col_span.a;
    const u32 in_col2 = col_span.b;
  
    return ( (in_row < n_rows) && (in_col1 <= in_col2) && (in_col2 < n_cols) );
    }
  }



template<typename oT>
arma_inline
arma_warn_unused
bool
field<oT>::in_range(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const u32 in_row1 = row_span.a;
  const u32 in_row2 = row_span.b;
  
  const u32 in_col1 = col_span.a;
  const u32 in_col2 = col_span.b;
  
  const bool rows_ok = row_span.whole ? true : ( (in_row1 <= in_row2) && (in_row2 < n_rows) );
  const bool cols_ok = col_span.whole ? true : ( (in_col1 <= in_col2) && (in_col2 < n_cols) );
  
  return ( (rows_ok == true) && (cols_ok == true) );
  }



template<typename oT>
inline
bool
field<oT>::save(const std::string name, const file_type type, const bool print_status) const
  {
  arma_extra_debug_sigprint();
  
  std::string err_msg;
  const bool save_okay = field_aux::save(*this, name, type, err_msg);
  
  if( (print_status == true) && (save_okay == false) )
    {
    if(err_msg.length() > 0)
      {
      arma_print("field::save(): ", err_msg, name);
      }
    else
      {
      arma_print("field::save(): couldn't write to ", name);
      }
    }
  
  return save_okay;
  }



template<typename oT>
inline
bool
field<oT>::save(std::ostream& os, const file_type type, const bool print_status) const
  {
  arma_extra_debug_sigprint();
  
  std::string err_msg;
  const bool save_okay = field_aux::save(*this, os, type, err_msg);
  
  if( (print_status == true) && (save_okay == false) )
    {
    if(err_msg.length() > 0)
      {
      arma_print("field::save(): ", err_msg, "[ostream]");
      }
    else
      {
      arma_print("field::save(): couldn't write to [ostream]");
      }
    }
  
  return save_okay;
  }



template<typename oT>
inline
bool
field<oT>::load(const std::string name, const file_type type, const bool print_status)
  {
  arma_extra_debug_sigprint();
  
  std::string err_msg;
  const bool load_okay = field_aux::load(*this, name, type, err_msg);
  
  if( (print_status == true) && (load_okay == false) )
    {
    if(err_msg.length() > 0)
      {
      arma_print("field::load(): ", err_msg, name);
      }
    else
      {
      arma_print("field::load(): couldn't read from ", name);
      }
    }
  
  if(load_okay == false)
    {
    (*this).reset();
    }
  
  return load_okay;
  }



template<typename oT>
inline
bool
field<oT>::load(std::istream& is, const file_type type, const bool print_status)
  {
  arma_extra_debug_sigprint();
  
  std::string err_msg;
  const bool load_okay = field_aux::load(*this, is, type, err_msg);
  
  if( (print_status == true) && (load_okay == false) )
    {
    if(err_msg.length() > 0)
      {
      arma_print("field::load(): ", err_msg, "[istream]");
      }
    else
      {
      arma_print("field::load(): couldn't read from [istream]");
      }
    }
  
  if(load_okay == false)
    {
    (*this).reset();
    }
  
  return load_okay;
  }



template<typename oT>
inline
bool
field<oT>::quiet_save(const std::string name, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).save(name, type, false);
  }



template<typename oT>
inline
bool
field<oT>::quiet_save(std::ostream& os, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).save(os, type, false);
  }



template<typename oT>
inline
bool
field<oT>::quiet_load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  return (*this).load(name, type, false);
  }



template<typename oT>
inline
bool
field<oT>::quiet_load(std::istream& is, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  return (*this).load(is, type, false);
  }



//! construct a field from a given field
template<typename oT>
inline
void
field<oT>::init(const field<oT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    init(x.n_rows, x.n_cols);
    
    field& t = *this;
    
    for(u32 col=0; col<x.n_cols; ++col)
    for(u32 row=0; row<x.n_rows; ++row)
      {
      t.at(row,col) = x.at(row,col);
      }
    }
  
  }



//! internal field construction; if the requested size is small enough, memory from the stack is used. otherwise memory is allocated via 'new'
template<typename oT>
inline
void
field<oT>::init(const u32 n_rows_in, const u32 n_cols_in)
  {
  arma_extra_debug_sigprint( arma_boost::format("n_rows_in = %d, n_cols_in = %d") % n_rows_in % n_cols_in );
  
  const u32 n_elem_new = n_rows_in * n_cols_in;

  if(n_elem == n_elem_new)
    {
    // delete_objects();
    // create_objects();
    access::rw(n_rows) = n_rows_in;
    access::rw(n_cols) = n_cols_in;
    }
  else
    {
    delete_objects();
    
    if(n_elem > sizeof(mem_local)/sizeof(oT*) )
      {
      delete [] mem;
      }
    
    if(n_elem_new <= sizeof(mem_local)/sizeof(oT*) )
      {
      mem = mem_local;
      }
    else
      {
      mem = new(std::nothrow) oT* [n_elem_new];
      arma_check( (mem == 0), "field::init(): out of memory" );
      }
    
    access::rw(n_elem) = n_elem_new;
    
    if(n_elem_new == 0)
      {
      access::rw(n_rows) = 0;
      access::rw(n_cols) = 0;
      }
    else
      {
      access::rw(n_rows) = n_rows_in;
      access::rw(n_cols) = n_cols_in;
      }
    
    create_objects();
    
    }
  
  }



template<typename oT>
inline
void
field<oT>::delete_objects()
  {
  arma_extra_debug_sigprint( arma_boost::format("n_elem = %d") % n_elem );
  
  for(u32 i=0; i<n_elem; ++i)
    {
    if(mem[i] != 0)
      {
      delete mem[i];
      mem[i] = 0;
      }
    }
  
  }



template<typename oT>
inline
void
field<oT>::create_objects()
  {
  arma_extra_debug_sigprint( arma_boost::format("n_elem = %d") % n_elem );
  
  for(u32 i=0; i<n_elem; ++i)
    {
    mem[i] = new oT;
    }
  
  }



template<typename oT>
inline
field<oT>::iterator::iterator(field<oT>& in_M, const bool at_end)
  : M(in_M)
  , i( (at_end == false) ? 0 : in_M.n_elem )
  {
  arma_extra_debug_sigprint();
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
  if(i > 0)
    {
    --i;
    }
  
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
  arma_extra_debug_sigprint();
  }



template<typename oT>
inline
field<oT>::const_iterator::const_iterator(const typename field<oT>::iterator& X)
  : M(X.M)
  , i(X.i)
  {
  arma_extra_debug_sigprint();
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
  if(i > 0)
    {
    --i;
    }
  
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
  arma_extra_debug_sigprint();
  
  return field<oT>::iterator(*this);
  }



template<typename oT>
inline
typename field<oT>::const_iterator
field<oT>::begin() const
  {
  arma_extra_debug_sigprint();
  
  return field<oT>::const_iterator(*this);
  }



template<typename oT>
inline
typename field<oT>::iterator
field<oT>::end()
  {
  arma_extra_debug_sigprint();
  
  return field<oT>::iterator(*this, true);
  }



template<typename oT>
inline
typename field<oT>::const_iterator
field<oT>::end() const
  {
  arma_extra_debug_sigprint();
  
  return field<oT>::const_iterator(*this, true);
  }
  


//
//
//



template<typename oT>
inline
void
field_aux::reset_objects(field<oT>& x)
  {
  arma_extra_debug_sigprint();
  
  x.delete_objects();
  x.create_objects();
  }



template<typename eT>
inline
void
field_aux::reset_objects(field< Mat<eT> >& x)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<x.n_elem; ++i)
    {
    (*(x.mem[i])).reset();
    }
  }



template<typename eT>
inline
void
field_aux::reset_objects(field< Col<eT> >& x)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<x.n_elem; ++i)
    {
    (*(x.mem[i])).reset();
    }
  }
  
  
  
template<typename eT>
inline
void
field_aux::reset_objects(field< Row<eT> >& x)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<x.n_elem; ++i)
    {
    (*(x.mem[i])).reset();
    }
  }



template<typename eT>
inline
void
field_aux::reset_objects(field< Cube<eT> >& x)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<x.n_elem; ++i)
    {
    (*(x.mem[i])).reset();
    }
  }



inline
void
field_aux::reset_objects(field< std::string >& x)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<x.n_elem; ++i)
    {
    (*(x.mem[i])).clear();
    }
  }



//
//
//



template<typename oT>
inline
bool
field_aux::save(const field<oT>& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  err_msg = " [sorry, saving/loading this type of field is currently not supported] filename = ";
  
  return false;
  }



template<typename oT>
inline
bool
field_aux::save(const field<oT>& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  err_msg = " [sorry, saving/loading this type of field is currently not supported] filename = ";
  
  return false;
  }



template<typename oT>
inline
bool
field_aux::load(field<oT>& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  err_msg = " [sorry, saving/loading this type of field is currently not supported] filename = ";
  
  return false;
  }



template<typename oT>
inline
bool
field_aux::load(field<oT>& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  err_msg = " [sorry, saving/loading this type of field is currently not supported] filename = ";
  
  return false;
  }



template<typename eT>
inline
bool
field_aux::save(const field< Mat<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, name);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, name);
      break;
    
    default:
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Mat<eT> >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, os);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, os);
      break;
    
    default:
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Mat<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
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
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Mat<eT> >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
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
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Col<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, name);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, name);
      break;
    
    default:
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Col<eT> >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, os);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, os);
      break;
    
    default:
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Col<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
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
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Col<eT> >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
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
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Row<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, name, err_msg);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, name, err_msg);
      break;
    
    default:
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Row<eT> >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, os, err_msg);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, os, err_msg);
      break;
    
    default:
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Row<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
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
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Row<eT> >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
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
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Cube<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, name, err_msg);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, name, err_msg);
      break;
    
    default:
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::save(const field< Cube<eT> >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case arma_binary:
      return diskio::save_arma_binary(x, os, err_msg);
      break;
      
    case ppm_binary:
      return diskio::save_ppm_binary(x, os, err_msg);
      break;
    
    default:
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Cube<eT> >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
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
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



template<typename eT>
inline
bool
field_aux::load(field< Cube<eT> >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
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
      err_msg = " [unsupported type] filename = ";
      return false;
    }
  }



inline
bool
field_aux::save(const field< std::string >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(type);
  
  err_msg.clear();
  
  return diskio::save_std_string(x, name);
  }



inline
bool
field_aux::save(const field< std::string >& x, std::ostream& os, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(type);
  
  err_msg.clear();
  
  return diskio::save_std_string(x, os);
  }



inline
bool
field_aux::load(field< std::string >& x, const std::string& name, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(type);
  
  return diskio::load_std_string(x, name, err_msg);
  }



inline
bool
field_aux::load(field< std::string >& x, std::istream& is, const file_type type, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(type);
  
  return diskio::load_std_string(x, is, err_msg);
  }



#ifdef ARMA_EXTRA_FIELD_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_FIELD_MEAT)
#endif



//! @}
