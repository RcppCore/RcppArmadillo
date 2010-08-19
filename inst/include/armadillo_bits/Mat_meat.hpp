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


//! \addtogroup Mat
//! @{


template<typename eT>
inline
Mat<eT>::~Mat()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(use_aux_mem == false)
    {
    if(n_elem > sizeof(mem_local)/sizeof(eT) )
      {
      delete [] mem;
      }
    }
    
  if(arma_config::debug == true)
    {
    // try to expose buggy user code that accesses deleted objects
    access::rw(n_rows) = 0;
    access::rw(n_cols) = 0;
    access::rw(n_elem) = 0;
    access::rw(mem)    = 0;
    }
  
  isnt_supported_elem_type<eT>::check();
  }



template<typename eT>
inline
Mat<eT>::Mat()
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  }


//! construct the matrix to have user specified dimensions
template<typename eT>
inline
Mat<eT>::Mat(const u32 in_n_rows, const u32 in_n_cols)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(in_n_rows, in_n_cols);
  }



//! internal matrix construction; if the requested size is small enough, memory from the stack is used. otherwise memory is allocated via 'new'
template<typename eT>
inline
void
Mat<eT>::init(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_n_rows = %d, in_n_cols = %d") % in_n_rows % in_n_cols );
  
  const u32 new_n_elem = in_n_rows * in_n_cols;

  if(n_elem == new_n_elem)
    {
    access::rw(n_rows) = in_n_rows;
    access::rw(n_cols) = in_n_cols;
    }
  else
    {
    arma_debug_check
      (
      (use_aux_mem == true),
      "Mat::init(): can't change the amount of memory as auxiliary memory is in use"
      );

    if(n_elem > sizeof(mem_local)/sizeof(eT) )
      {
      delete [] mem;
      }
    
    if(new_n_elem <= sizeof(mem_local)/sizeof(eT) )
      {
      access::rw(mem) = mem_local;
      }
    else
      {
      access::rw(mem) = new(std::nothrow) eT[new_n_elem];
      arma_check( (mem == 0), "Mat::init(): out of memory" );
      }
    
    access::rw(n_elem) = new_n_elem;

    if(new_n_elem == 0)
      {
      access::rw(n_rows) = 0;
      access::rw(n_cols) = 0;
      }
    else
      {
      access::rw(n_rows) = in_n_rows;
      access::rw(n_cols) = in_n_cols;
      }
    
    }
  }


//! create the matrix from a textual description
template<typename eT>
inline
Mat<eT>::Mat(const char* text)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init( std::string(text) );
  }
  
  
  
//! create the matrix from a textual description
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  init( std::string(text) );
  return *this;
  }
  
  

//! create the matrix from a textual description
template<typename eT>
inline
Mat<eT>::Mat(const std::string& text)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(text);
  }
  
  
  
//! create the matrix from a textual description
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  init(text);
  return *this;
  }



//! internal function to create the matrix from a textual description
template<typename eT>
inline 
void
Mat<eT>::init(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  //
  // work out the size
  
  u32 t_n_rows = 0;
  u32 t_n_cols = 0;
  
  bool t_n_cols_found = false;
  
  std::string token;
  
  std::string::size_type line_start = 0;
  std::string::size_type   line_end = 0;
  
  while( line_start < text.length() )
    {
    
    line_end = text.find(';', line_start);
    
    if(line_end == std::string::npos)
      line_end = text.length()-1;
    
    std::string::size_type line_len = line_end - line_start + 1;
    std::stringstream line_stream( text.substr(line_start,line_len) );
    
    
    u32 line_n_cols = 0;
    while(line_stream >> token)
      {
      ++line_n_cols;
      }
    
    
    if(line_n_cols > 0)
      {
      if(t_n_cols_found == false)
        {
        t_n_cols = line_n_cols;
        t_n_cols_found = true;
        }
      else
        arma_check( (line_n_cols != t_n_cols), "Mat::init(): inconsistent number of columns in given string");
      
      ++t_n_rows;
      }
    line_start = line_end+1;
    
    }
    
  Mat<eT>& x = *this;
  x.set_size(t_n_rows, t_n_cols);
  
  line_start = 0;
  line_end = 0;
  
  u32 row = 0;
  
  while( line_start < text.length() )
    {
    
    line_end = text.find(';', line_start);
    
    if(line_end == std::string::npos)
      line_end = text.length()-1;
    
    std::string::size_type line_len = line_end - line_start + 1;
    std::stringstream line_stream( text.substr(line_start,line_len) );
    
//     u32 col = 0;
//     while(line_stream >> token)
//       {
//       x.at(row,col) = strtod(token.c_str(), 0);
//       ++col;
//       }
    
    u32 col = 0;
    eT val;
    while(line_stream >> val)
      {
      x.at(row,col) = val;
      ++col;
      }
    
    ++row;
    line_start = line_end+1;
    }
  
  }



//! Set the matrix to be equal to the specified scalar.
//! NOTE: the size of the matrix will be 1x1
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  init(1,1);
  access::rw(mem[0]) = val;
  return *this;
  }



//! In-place addition of a scalar to all elements of the matrix
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator+=(const eT val)
  {
  arma_extra_debug_sigprint();
  
        eT* local_ptr    = memptr();
  const u32 local_n_elem = n_elem;
    
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    local_ptr[i] += val;
    local_ptr[j] += val;
    }
  
  if(i < local_n_elem)
    {
    local_ptr[i] += val;
    }
  
  return *this;
  }



//! In-place subtraction of a scalar from all elements of the matrix
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();
  
        eT* local_ptr    = memptr();
  const u32 local_n_elem = n_elem;
    
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    local_ptr[i] -= val;
    local_ptr[j] -= val;
    }
  
  if(i < local_n_elem)
    {
    local_ptr[i] -= val;
    }
  
  return *this;
  }



//! In-place multiplication of all elements of the matrix with a scalar
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
        eT* local_ptr    = memptr();
  const u32 local_n_elem = n_elem;
    
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    local_ptr[i] *= val;
    local_ptr[j] *= val;
    }
  
  if(i < local_n_elem)
    {
    local_ptr[i] *= val;
    }
  
  return *this;
  }



//! In-place division of all elements of the matrix with a scalar
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
        eT* local_ptr    = memptr();
  const u32 local_n_elem = n_elem;
    
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    local_ptr[i] /= val;
    local_ptr[j] /= val;
    }
  
  if(i < local_n_elem)
    {
    local_ptr[i] /= val;
    }
  
  return *this;
  }



//! construct a matrix from a given matrix
template<typename eT>
inline
Mat<eT>::Mat(const Mat<eT>& in_mat)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint(arma_boost::format("this = %x   in_mat = %x") % this % &in_mat);
  
  init(in_mat);
  }



//! construct a matrix from a given matrix
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const Mat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  init(x);
  return *this;
  }



//! construct a matrix from a given matrix
template<typename eT>
inline
void
Mat<eT>::init(const Mat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    init(x.n_rows, x.n_cols);
    syslib::copy_elem( memptr(), x.mem, n_elem );
    }
  }



//! construct a matrix from a given auxiliary array of eTs.
//! if copy_aux_mem is true, new memory is allocated and the array is copied.
//! if copy_aux_mem is false, the auxiliary array is used directly (without allocating memory and copying).
//! the default is to copy the array.

template<typename eT>
inline
Mat<eT>::Mat(eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const bool copy_aux_mem)
  : n_rows     (copy_aux_mem ? 0     : aux_n_rows           )
  , n_cols     (copy_aux_mem ? 0     : aux_n_cols           )
  , n_elem     (copy_aux_mem ? 0     : aux_n_rows*aux_n_cols)
  , use_aux_mem(copy_aux_mem ? false : true                 )
  , mem        (copy_aux_mem ? mem   : aux_mem              )
  {
  arma_extra_debug_sigprint_this(this);
  
  if(copy_aux_mem == true)
    {
    init(aux_n_rows, aux_n_cols);
    syslib::copy_elem( memptr(), aux_mem, n_elem );
    }
  }



//! construct a matrix from a given auxiliary read-only array of eTs.
//! the array is copied.
template<typename eT>
inline
Mat<eT>::Mat(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(aux_n_rows, aux_n_cols);
  syslib::copy_elem( memptr(), aux_mem, n_elem );
  }



//! DANGEROUS! Construct a temporary matrix, using auxiliary memory.
//! This constructor is NOT intended for usage by user code.
//! Its sole purpose is to be used by the Cube class.

template<typename eT>
inline
Mat<eT>::Mat(const char junk, const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols)
  : n_rows     (aux_n_rows           )
  , n_cols     (aux_n_cols           )
  , n_elem     (aux_n_rows*aux_n_cols)
  , use_aux_mem(true                 )
  , mem        (aux_mem              )
  {
  arma_extra_debug_sigprint_this(this);
  }



//! in-place matrix addition
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator+=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "matrix addition");
  
  const u32 local_n_elem = m.n_elem;
  
        eT* out_mem = (*this).memptr();
  const eT* m_mem   = m.memptr();
  
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    out_mem[i] += m_mem[i];
    out_mem[j] += m_mem[j];
    }
  
  if(i < local_n_elem)
    {
    out_mem[i] += m_mem[i];
    }
  
  return *this;
  }



//! in-place matrix subtraction
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator-=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "matrix subtraction");
  
  const u32 local_n_elem = m.n_elem;
  
        eT* out_mem = (*this).memptr();
  const eT* m_mem   = m.memptr();
  
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    out_mem[i] -= m_mem[i];
    out_mem[j] -= m_mem[j];
    }
  
  if(i < local_n_elem)
    {
    out_mem[i] -= m_mem[i];
    }
  
  return *this;
  }



//! in-place matrix multiplication
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator*=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  glue_times::apply_inplace(*this, m);
  return *this;
  }



//! in-place element-wise matrix multiplication
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator%=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "element-wise matrix multplication");
  
  const u32 local_n_elem = m.n_elem;
  
        eT* out_mem = (*this).memptr();
  const eT* m_mem   = m.memptr();
  
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    out_mem[i] *= m_mem[i];
    out_mem[j] *= m_mem[j];
    }
  
  if(i < local_n_elem)
    {
    out_mem[i] *= m_mem[i];
    }
  
  return *this;
  }



//! in-place element-wise matrix division
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator/=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "element-wise matrix division");
  
  const u32 local_n_elem = m.n_elem;
  
        eT* out_mem = (*this).memptr();
  const eT* m_mem   = m.memptr();
  
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    out_mem[i] /= m_mem[i];
    out_mem[j] /= m_mem[j];
    }
  
  if(i < local_n_elem)
    {
    out_mem[i] /= m_mem[i];
    }
  
  return *this;
  }



//! for constructing a complex matrix out of two non-complex matrices
template<typename eT>
template<typename T1, typename T2>
inline
Mat<eT>::Mat
  (
  const Base<typename Mat<eT>::pod_type,T1>& A,
  const Base<typename Mat<eT>::pod_type,T2>& B
  )
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  arma_type_check< is_complex<eT>::value == false >::apply();   //!< compile-time abort if eT isn't std::complex
  
  typedef typename T1::elem_type T;
  arma_type_check< is_complex<T>::value == true >::apply();   //!< compile-time abort if T is std::complex
  
  isnt_same_type<std::complex<T>, eT>::check();   //!< compile-time abort if types are not compatible
  
  const unwrap<T1> tmp_A(A.get_ref());
  const unwrap<T2> tmp_B(B.get_ref());
  
  const Mat<T>& X = tmp_A.M;
  const Mat<T>& Y = tmp_B.M;
  
  arma_assert_same_size(X, Y, "Mat()");
  
  init(X.n_rows, Y.n_cols);
  
  const T* X_mem = X.mem;
  const T* Y_mem = Y.mem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) = std::complex<T>(X_mem[i], Y_mem[i]);
    }
  }



//! construct a matrix from subview (e.g. construct a matrix from a delayed submatrix operation)
template<typename eT>
inline
Mat<eT>::Mat(const subview<eT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  this->operator=(X);
  }



//! construct a matrix from subview (e.g. construct a matrix from a delayed submatrix operation)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::extract(*this, X);
  return *this;
  }


//! in-place matrix addition (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator+=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::plus_inplace(*this, X);
  return *this;
  }


//! in-place matrix subtraction (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator-=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::minus_inplace(*this, X);
  return *this;
  }



//! in-place matrix mutiplication (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator*=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  glue_times::apply_inplace(*this, X);
  return *this;
  }



//! in-place element-wise matrix mutiplication (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator%=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::schur_inplace(*this, X);
  return *this;
  }



//! in-place element-wise matrix division (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator/=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::div_inplace(*this, X);
  return *this;
  }



//! construct a matrix from a subview_cube instance
template<typename eT>
inline
Mat<eT>::Mat(const subview_cube<eT>& x)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  this->operator=(x);
  }



//! construct a matrix from a subview_cube instance
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::extract(*this, X);
  return *this;
  }



//! in-place matrix addition (using a single-slice subcube on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator+=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();

  subview_cube<eT>::plus_inplace(*this, X);
  return *this;
  }



//! in-place matrix subtraction (using a single-slice subcube on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator-=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::minus_inplace(*this, X);
  return *this;
  }



//! in-place matrix mutiplication (using a single-slice subcube on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator*=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();

  const Mat<eT> tmp(X);
  glue_times::apply_inplace(*this, tmp);
  return *this;
  }



//! in-place element-wise matrix mutiplication (using a single-slice subcube on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator%=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::schur_inplace(*this, X);
  return *this;
  }



//! in-place element-wise matrix division (using a single-slice subcube on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator/=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::div_inplace(*this, X);
  return *this;
  }



//! construct a matrix from diagview (e.g. construct a matrix from a delayed diag operation)
template<typename eT>
inline
Mat<eT>::Mat(const diagview<eT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  this->operator=(X);
  }



//! construct a matrix from diagview (e.g. construct a matrix from a delayed diag operation)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>::extract(*this, X);
  return *this;
  }



//! in-place matrix addition (using a diagview on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator+=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>::plus_inplace(*this, X);
  return *this;
  }


//! in-place matrix subtraction (using a diagview on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator-=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>::minus_inplace(*this, X);
  return *this;
  }



//! in-place matrix mutiplication (using a diagview on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator*=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  glue_times::apply_inplace(*this, X);
  return *this;
  }



//! in-place element-wise matrix mutiplication (using a diagview on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator%=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>::schur_inplace(*this, X);
  return *this;
  }



//! in-place element-wise matrix division (using a diagview on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator/=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>::div_inplace(*this, X);
  return *this;
  }



template<typename eT>
inline
mat_injector< Mat<eT> >
Mat<eT>::operator<<(const eT val)
  {
  return mat_injector< Mat<eT> >(*this, val);
  }



template<typename eT>
inline
mat_injector< Mat<eT> >
Mat<eT>::operator<<(const injector_helper x)
  {
  return mat_injector< Mat<eT> >(*this, x);
  }



//! creation of subview (row vector)
template<typename eT>
arma_inline
subview_row<eT>
Mat<eT>::row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= n_rows, "Mat::row(): row out of bounds" );
  
  return subview_row<eT>(*this, row_num);
  }



//! creation of subview (row vector)
template<typename eT>
arma_inline
const subview_row<eT>
Mat<eT>::row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= n_rows, "Mat::row(): row out of bounds" );
  
  return subview_row<eT>(*this, row_num);
  }



//! creation of subview (column vector)
template<typename eT>
arma_inline
subview_col<eT>
Mat<eT>::col(const u32 col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "Mat::col(): out of bounds");
  
  return subview_col<eT>(*this, col_num);
  }



//! creation of subview (column vector)
template<typename eT>
arma_inline
const subview_col<eT>
Mat<eT>::col(const u32 col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "Mat::col(): out of bounds");
  
  return subview_col<eT>(*this, col_num);
  }



//! creation of subview (submatrix comprised of specified row vectors)
template<typename eT>
arma_inline
subview<eT>
Mat<eT>::rows(const u32 in_row1, const u32 in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "Mat::rows(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, in_row1, 0, in_row2, ((n_cols>0) ? n_cols-1 : 0) );
  }



//! creation of subview (submatrix comprised of specified row vectors)
template<typename eT>
arma_inline
const subview<eT>
Mat<eT>::rows(const u32 in_row1, const u32 in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "Mat::rows(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, in_row1, 0, in_row2, ((n_cols>0) ? n_cols-1 : 0) );
  }



//! creation of subview (submatrix comprised of specified column vectors)
template<typename eT>
arma_inline
subview<eT>
Mat<eT>::cols(const u32 in_col1, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "Mat::cols(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, 0, in_col1, ((n_rows>0) ? n_rows-1 : 0), in_col2);
  }



//! creation of subview (submatrix comprised of specified column vectors)
template<typename eT>
arma_inline
const subview<eT>
Mat<eT>::cols(const u32 in_col1, const u32 in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "Mat::cols(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, 0, in_col1, ((n_rows>0) ? n_rows-1 : 0), in_col2);
  }



//! creation of subview (submatrix)
template<typename eT>
arma_inline
subview<eT>
Mat<eT>::submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 >  in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "Mat::submat(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, in_row1, in_col1, in_row2, in_col2);
  }



//! creation of subview (generic submatrix)
template<typename eT>
arma_inline
const subview<eT>
Mat<eT>::submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 >  in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "Mat::submat(): indices out of bounds or incorrectly used"
    );
    
  return subview<eT>(*this, in_row1, in_col1, in_row2, in_col2);
  }



//! creation of diagview (diagonal)
template<typename eT>
arma_inline
diagview<eT>
Mat<eT>::diag(const s32 in_id)
  {
  arma_extra_debug_sigprint();
  
  const u32 row_offset = (in_id < 0) ? -in_id : 0;
  const u32 col_offset = (in_id > 0) ?  in_id : 0;
  
  arma_debug_check
    (
    (row_offset >= n_rows) || (col_offset >= n_cols),
    "Mat::diag(): requested diagonal out of bounds"
    );
  
  const u32 len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  return diagview<eT>(*this, row_offset, col_offset, len);
  }



//! creation of diagview (diagonal)
template<typename eT>
arma_inline
const diagview<eT>
Mat<eT>::diag(const s32 in_id) const
  {
  arma_extra_debug_sigprint();
  
  const u32 row_offset = (in_id < 0) ? -in_id : 0;
  const u32 col_offset = (in_id > 0) ?  in_id : 0;
  
  arma_debug_check
    (
    (row_offset >= n_rows) || (col_offset >= n_cols),
    "Mat::diag(): requested diagonal out of bounds"
    );
  
  
  const u32 len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  return diagview<eT>(*this, row_offset, col_offset, len);
  }



template<typename eT>
inline
void
Mat<eT>::swap_rows(const u32 in_row1, const u32 in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 >= n_rows) || (in_row2 >= n_rows),
    "Mat::swap_rows(): out of bounds"
    );
  
  for(u32 col=0; col<n_cols; ++col)
    {
    const u32 offset = col*n_rows;
    const u32 pos1   = in_row1 + offset;
    const u32 pos2   = in_row2 + offset;
    
    const eT tmp          = mem[pos1];
    access::rw(mem[pos1]) = mem[pos2];
    access::rw(mem[pos2]) = tmp;
    }
  
  }



template<typename eT>
inline
void
Mat<eT>::swap_cols(const u32 in_col1, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 >= n_cols) || (in_col2 >= n_cols),
    "Mat::swap_cols(): out of bounds"
    );
  
  eT* ptr1 = colptr(in_col1);
  eT* ptr2 = colptr(in_col2);
  
  for(u32 row=0; row<n_rows; ++row)
    {
    const eT tmp = ptr1[row];
    ptr1[row]    = ptr2[row];
    ptr2[row]    = tmp;
    }
  
  }



//! create a matrix from Op, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename op_type>
inline
Mat<eT>::Mat(const Op<T1, op_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  op_type::apply(*this, X);
  }



//! create a matrix from Op, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  op_type::apply(*this, X);
  
  return *this;
  }



//! in-place matrix addition, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator+=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const Mat<eT> m(X);
  
  return (*this).operator+=(m);
  }



//! in-place matrix subtraction, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator-=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const Mat<eT> m(X);
  
  return (*this).operator-=(m);
  }



//! in-place matrix multiplication, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator*=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  glue_times::apply_inplace(*this, X);
  
  return *this;
  }



//! in-place matrix element-wise multiplication, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator%=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const Mat<eT> m(X);
  
  return (*this).operator%=(m);
  }



//! in-place matrix element-wise division, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator/=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const Mat<eT> m(X);
  
  return (*this).operator/=(m);
  }



//! create a matrix from eOp, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename eop_type>
inline
Mat<eT>::Mat(const eOp<T1, eop_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply(*this, X);
  }



//! create a matrix from eOp, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename eop_type>
inline
const Mat<eT>&
Mat<eT>::operator=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Mat<eT>&
Mat<eT>::operator+=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply_inplace_plus(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Mat<eT>&
Mat<eT>::operator-=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply_inplace_minus(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Mat<eT>&
Mat<eT>::operator*=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  glue_times::apply_inplace(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Mat<eT>&
Mat<eT>::operator%=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply_inplace_schur(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1, typename eop_type>
inline
const Mat<eT>&
Mat<eT>::operator/=(const eOp<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply_inplace_div(*this, X);
  
  return *this;
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
Mat<eT>::Mat(const mtOp<eT, T1, op_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  op_type::apply(*this, X);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  op_type::apply(*this, X);
  
  return *this;
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator+=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator+=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator-=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator-=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator*=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator*=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator%=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator%=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator/=(const mtOp<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator/=(m);
  }



//! create a matrix from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Mat<eT>::Mat(const Glue<T1, T2, glue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_type::apply(*this, X);
  }



//! create a matrix from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_type::apply(*this, X);
  
  return *this;
  }



//! in-place matrix addition, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator+=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  const Mat<eT> m(X);
  
  return (*this).operator+=(m);
  }



//! in-place matrix subtraction, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator-=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  const Mat<eT> m(X);
  
  return (*this).operator-=(m);
  }



//! in-place matrix multiplications, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator*=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_times::apply_inplace(*this, X);
  
  return *this;
  }



//! in-place matrix element-wise multiplication, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator%=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  const Mat<eT> m(X);
  
  return (*this).operator%=(m);
  }



//! in-place matrix element-wise division, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator/=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  const Mat<eT> m(X);
  
  return (*this).operator/=(m);
  }



template<typename eT>
template<typename T1, typename T2>
inline
const Mat<eT>&
Mat<eT>::operator+=(const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  glue_times::apply_inplace_plus(*this, X, s32(+1));
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2>
inline
const Mat<eT>&
Mat<eT>::operator-=(const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  glue_times::apply_inplace_plus(*this, X, s32(-1));
  
  return *this;
  }



//! create a matrix from eGlue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
Mat<eT>::Mat(const eGlue<T1, T2, eglue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply(*this, X);
  }



//! create a matrix from eGlue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Mat<eT>&
Mat<eT>::operator=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply(*this, X);
  
  return *this;
  }



//! in-place matrix addition, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Mat<eT>&
Mat<eT>::operator+=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply_inplace_plus(*this, X);
  
  return *this;
  }



//! in-place matrix subtraction, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Mat<eT>&
Mat<eT>::operator-=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply_inplace_minus(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Mat<eT>&
Mat<eT>::operator*=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_times::apply_inplace(*this, X);
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Mat<eT>&
Mat<eT>::operator%=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply_inplace_schur(*this, X);
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Mat<eT>&
Mat<eT>::operator/=(const eGlue<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply_inplace_div(*this, X);
  return *this;
  }



//! EXPERIMENTAL: create a matrix from mtGlue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Mat<eT>::Mat(const mtGlue<eT, T1, T2, glue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , use_aux_mem(false)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  glue_type::apply(*this, X);
  }



//! EXPERIMENTAL: create a matrix from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  glue_type::apply(*this, X);
  
  return *this;
  }



//! EXPERIMENTAL: in-place matrix addition, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator+=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator+=(m);
  }



//! EXPERIMENTAL: in-place matrix subtraction, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator-=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator-=(m);
  }



//! EXPERIMENTAL: in-place matrix multiplications, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator*=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  glue_times::apply_inplace(*this, m);
  
  return *this;
  }



//! EXPERIMENTAL: in-place matrix element-wise multiplication, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator%=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator%=(m);
  }



//! EXPERIMENTAL: in-place matrix element-wise division, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator/=(const mtGlue<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> m(X);
  
  return (*this).operator/=(m);
  }



//! linear element accessor (treats the matrix as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT&
Mat<eT>::operator() (const u32 i)
  {
  arma_debug_check( (i >= n_elem), "Mat::operator(): out of bounds");
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the matrix as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT
Mat<eT>::operator() (const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "Mat::operator(): out of bounds");
  return mem[i];
  }


//! linear element accessor (treats the matrix as a vector); no bounds check.  
template<typename eT>
arma_inline
eT&
Mat<eT>::operator[] (const u32 i)
  {
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the matrix as a vector); no bounds check
template<typename eT>
arma_inline
eT
Mat<eT>::operator[] (const u32 i) const
  {
  return mem[i];
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT&
Mat<eT>::operator() (const u32 in_row, const u32 in_col)
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "Mat::operator(): out of bounds");
  return access::rw(mem[in_row + in_col*n_rows]);
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT
Mat<eT>::operator() (const u32 in_row, const u32 in_col) const
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "Mat::operator(): out of bounds");
  return mem[in_row + in_col*n_rows];
  }



//! element accessor; no bounds check
template<typename eT>
arma_inline
eT&
Mat<eT>::at(const u32 in_row, const u32 in_col)
  {
  return access::rw( mem[in_row + in_col*n_rows] );
  }



//! element accessor; no bounds check
template<typename eT>
arma_inline
eT
Mat<eT>::at(const u32 in_row, const u32 in_col) const
  {
  return mem[in_row + in_col*n_rows];
  }



//! prefix ++
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator++()
  {
  Mat_aux::prefix_pp(*this);
  return *this;
  }



//! postfix ++  (must not return the object by reference)
template<typename eT>
arma_inline
void
Mat<eT>::operator++(int)
  {
  Mat_aux::postfix_pp(*this);
  }



//! prefix --
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator--()
  {
  Mat_aux::prefix_mm(*this);
  return *this;
  }



//! postfix --  (must not return the object by reference)
template<typename eT>
arma_inline
void
Mat<eT>::operator--(int)
  {
  Mat_aux::postfix_mm(*this);
  }



//! returns true if the object can be interpreted as a column or row vector
template<typename eT>
arma_inline
bool
Mat<eT>::is_vec() const
  {
  return ( (n_rows == 1) || (n_cols == 1) );
  }



//! returns true if the object has the same number of non-zero rows and columnns
template<typename eT>
arma_inline
bool
Mat<eT>::is_square() const
  {
  return ( (n_rows == n_cols) && (n_elem > 0) );
  }



//! returns true if all of the elements are finite
template<typename eT>
arma_inline
bool
Mat<eT>::is_finite() const
  {
  for(u32 i=0; i<n_elem; ++i)
    {
    if(arma_isfinite(mem[i]) == false)
      {
      return false;
      }
    }

  return true;
  }



//! returns a pointer to array of eTs for a specified column; no bounds check
template<typename eT>
arma_inline
eT*
Mat<eT>::colptr(const u32 in_col)
  {
  return & access::rw(mem[in_col*n_rows]);
  }



//! returns a pointer to array of eTs for a specified column; no bounds check
template<typename eT>
arma_inline
const eT*
Mat<eT>::colptr(const u32 in_col) const
  {
  return & mem[in_col*n_rows];
  }



//! returns a pointer to array of eTs used by the matrix
template<typename eT>
arma_inline
eT*
Mat<eT>::memptr()
  {
  return const_cast<eT*>(mem);
  }



//! returns a pointer to array of eTs used by the matrix
template<typename eT>
arma_inline
const eT*
Mat<eT>::memptr() const
  {
  return mem;
  }



//! print contents of the matrix (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's flags are restored to their original values.
template<typename eT>
inline
void
Mat<eT>::print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = cout.width();
    
    cout << extra_text << '\n';
  
    cout.width(orig_width);
    }
  
  arma_ostream::print(cout, *this, true);
  }



//! print contents of the matrix to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's flags are restored to their original values.
template<typename eT>
inline
void
Mat<eT>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
    
    user_stream.width(orig_width);
    }
  
  arma_ostream::print(user_stream, *this, true);
  }



//! print contents of the transposed version of the matrix (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's flags are restored to their original values.
template<typename eT>
inline
void
Mat<eT>::print_trans(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  op_trans::apply_noalias(tmp, *this);
  
  tmp.print(extra_text);
  }



//! print contents of the transposed version of matrix to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's flags are restored to their original values.
template<typename eT>
inline
void
Mat<eT>::print_trans(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  op_trans::apply_noalias(tmp, *this);
  
  tmp.print(user_stream, extra_text);
  }



//! print contents of the matrix (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the stream's flags are used as is and are not modified
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Mat<eT>::raw_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = cout.width();
    
    cout << extra_text << '\n';
  
    cout.width(orig_width);
    }
  
  arma_ostream::print(cout, *this, false);
  }



//! print contents of the matrix to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the stream's flags are used as is and are not modified.
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Mat<eT>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
  
    user_stream << extra_text << '\n';
  
    user_stream.width(orig_width);
    }
  
  arma_ostream::print(user_stream, *this, false);
  }



//! print contents of the transposed version of the matrix (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the stream's flags are used as is and are not modified
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Mat<eT>::raw_print_trans(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  op_trans::apply_noalias(tmp, *this);
  
  tmp.raw_print(extra_text);
  }



//! print contents of the transposed version of the matrix to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the stream's flags are used as is and are not modified.
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Mat<eT>::raw_print_trans(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  op_trans::apply_noalias(tmp, *this);
  
  tmp.raw_print(user_stream, extra_text);
  }



//! change the matrix to have user specified dimensions (data is not preserved)
template<typename eT>
inline
void
Mat<eT>::set_size(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint();
  
  init(in_rows, in_cols);
  }



//! change the matrix to have user specified dimensions (data is preserved)
template<typename eT>
inline
void
Mat<eT>::reshape(const u32 in_rows, const u32 in_cols, const u32 dim)
  {
  arma_extra_debug_sigprint();
  
  *this = arma::reshape(*this, in_rows, in_cols, dim);
  }



//! change the matrix (without preserving data) to have the same dimensions as the given matrix 
template<typename eT>
template<typename eT2>
inline
void
Mat<eT>::copy_size(const Mat<eT2>& m)
  {
  arma_extra_debug_sigprint();
  
  init(m.n_rows, m.n_cols);
  }



//! fill the matrix with the specified value
template<typename eT>
arma_hot
inline
void
Mat<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
        eT* local_ptr    = memptr();
  const u32 local_n_elem = n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<local_n_elem; i+=2, j+=2)
    {
    local_ptr[i] = val;
    local_ptr[j] = val;
    }
  
  if(i < local_n_elem)
    {
    local_ptr[i] = val;
    }
  }



template<typename eT>
inline
void
Mat<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(0));
  }



template<typename eT>
inline
void
Mat<eT>::zeros(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_rows = %d, in_cols = %d") % in_rows % in_cols );

  set_size(in_rows, in_cols);
  fill(eT(0));
  }



template<typename eT>
inline
void
Mat<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(1));
  }



template<typename eT>
inline
void
Mat<eT>::ones(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_rows = %d, in_cols = %d") % in_rows % in_cols );

  set_size(in_rows, in_cols);
  fill(eT(1));
  }



template<typename eT>
inline
void
Mat<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  init(0,0);
  }



//! save the matrix to a file
template<typename eT>
inline
bool
Mat<eT>::save(const std::string name, const file_type type, const bool print_status) const
  {
  arma_extra_debug_sigprint();
  
  bool save_okay;
  
  switch(type)
    {
    case raw_ascii:
      save_okay = diskio::save_raw_ascii(*this, name);
      break;
    
    case arma_ascii:
      save_okay = diskio::save_arma_ascii(*this, name);
      break;
    
    case arma_binary:
      save_okay = diskio::save_arma_binary(*this, name);
      break;
      
    case pgm_binary:
      save_okay = diskio::save_pgm_binary(*this, name);
      break;
    
    default:
      arma_warn(print_status, "Mat::save(): unsupported file type");
      save_okay = false;
    }
  
  arma_warn( (print_status && (save_okay == false)), "Mat::save(): couldn't write to ", name);
  
  return save_okay;
  }



//! save the matrix to a stream
template<typename eT>
inline
bool
Mat<eT>::save(std::ostream& os, const file_type type, const bool print_status) const
  {
  arma_extra_debug_sigprint();
  
  bool save_okay;
  
  switch(type)
    {
    case raw_ascii:
      save_okay = diskio::save_raw_ascii(*this, os);
      break;
    
    case arma_ascii:
      save_okay = diskio::save_arma_ascii(*this, os);
      break;
    
    case arma_binary:
      save_okay = diskio::save_arma_binary(*this, os);
      break;
      
    case pgm_binary:
      save_okay = diskio::save_pgm_binary(*this, os);
      break;
    
    default:
      arma_warn(print_status, "Mat::save(): unsupported file type");
      save_okay = false;
    }
  
  arma_warn( (print_status && (save_okay == false)), "Mat::save(): couldn't write to the given stream");
  
  return save_okay;
  }



//! load a matrix from a file
template<typename eT>
inline
bool
Mat<eT>::load(const std::string name, const file_type type, const bool print_status)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay;
  std::string err_msg;
  
  switch(type)
    {
    case auto_detect:
      load_okay = diskio::load_auto_detect(*this, name, err_msg);
      break;
    
    case raw_ascii:
      load_okay = diskio::load_raw_ascii(*this, name, err_msg);
      break;
    
    case arma_ascii:
      load_okay = diskio::load_arma_ascii(*this, name, err_msg);
      break;
    
    case arma_binary:
      load_okay = diskio::load_arma_binary(*this, name, err_msg);
      break;
      
    case pgm_binary:
      load_okay = diskio::load_pgm_binary(*this, name, err_msg);
      break;
    
    default:
      arma_warn(print_status, "Mat::load(): unsupported file type");
      load_okay = false;
    }
  
  if( (print_status == true) && (load_okay == false) )
    {
    if(err_msg.length() > 0)
      {
      arma_print("Mat::load(): ", err_msg, name);
      }
    else
      {
      arma_print("Mat::load(): couldn't read ", name);
      }
    }
  
  if(load_okay == false)
    {
    (*this).reset();
    }
    
  return load_okay;
  }



//! load a matrix from a stream
template<typename eT>
inline
bool
Mat<eT>::load(std::istream& is, const file_type type, const bool print_status)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay;
  std::string err_msg;
  
  switch(type)
    {
    case auto_detect:
      load_okay = diskio::load_auto_detect(*this, is, err_msg);
      break;
    
    case raw_ascii:
      load_okay = diskio::load_raw_ascii(*this, is, err_msg);
      break;
    
    case arma_ascii:
      load_okay = diskio::load_arma_ascii(*this, is, err_msg);
      break;
    
    case arma_binary:
      load_okay = diskio::load_arma_binary(*this, is, err_msg);
      break;
      
    case pgm_binary:
      load_okay = diskio::load_pgm_binary(*this, is, err_msg);
      break;
    
    default:
      arma_warn(print_status, "Mat::load(): unsupported file type");
      load_okay = false;
    }
  
  
  if( (print_status == true) && (load_okay == false) )
    {
    if(err_msg.length() > 0)
      {
      arma_print("Mat::load(): ", err_msg, "the given stream");
      }
    else
      {
      arma_print("Mat::load(): couldn't load from the given stream");
      }
    }
  
  if(load_okay == false)
    {
    (*this).reset();
    }
    
  return load_okay;
  }



//! save the matrix to a file, without printing any error messages
template<typename eT>
inline
bool
Mat<eT>::quiet_save(const std::string name, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).save(name, type, false);
  }



//! save the matrix to a stream, without printing any error messages
template<typename eT>
inline
bool
Mat<eT>::quiet_save(std::ostream& os, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).save(os, type, false);
  }



//! load a matrix from a file, without printing any error messages
template<typename eT>
inline
bool
Mat<eT>::quiet_load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  return (*this).load(name, type, false);
  }



//! load a matrix from a stream, without printing any error messages
template<typename eT>
inline
bool
Mat<eT>::quiet_load(std::istream& is, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  return (*this).load(is, type, false);
  }



template<typename eT>
inline
Mat<eT>::row_iterator::row_iterator(Mat<eT>& in_M, const u32 in_row)
  : M  (in_M  )
  , row(in_row)
  , col(0     )
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
eT&
Mat<eT>::row_iterator::operator*()
  {
  return M.at(row,col);
  }



template<typename eT>
inline
typename Mat<eT>::row_iterator&
Mat<eT>::row_iterator::operator++()
  {
  ++col;
  
  if(col >= M.n_cols)
    {
    col = 0;
    ++row;
    }
  
  return *this;
  }



template<typename eT>
inline
void
Mat<eT>::row_iterator::operator++(int)
  {
  operator++();
  }



template<typename eT>
inline
typename Mat<eT>::row_iterator&
Mat<eT>::row_iterator::operator--()
  {
  if(col > 0)
    {
    --col;
    }
  else
    {
    if(row > 0)
      {
      col = M.n_cols - 1;
      --row;
      }
    }
  
  return *this;
  }



template<typename eT>
inline
void
Mat<eT>::row_iterator::operator--(int)
  {
  operator--();
  }



template<typename eT>
inline
bool
Mat<eT>::row_iterator::operator!=(const typename Mat<eT>::row_iterator& X) const
  {
  return ( (row != X.row) || (col != X.col) ) ? true : false;
  }



template<typename eT>
inline
bool
Mat<eT>::row_iterator::operator==(const typename Mat<eT>::row_iterator& X) const
  {
  return ( (row == X.row) && (col == X.col) ) ? true : false;
  }



template<typename eT>
inline
Mat<eT>::const_row_iterator::const_row_iterator(const Mat<eT>& in_M, const u32 in_row)
  : M  (in_M  )
  , row(in_row)
  , col(0     )
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
Mat<eT>::const_row_iterator::const_row_iterator(const typename Mat<eT>::row_iterator& X)
  : M  (X.M)
  , row(X.row)
  , col(X.col)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
eT
Mat<eT>::const_row_iterator::operator*() const
  {
  return M.at(row,col);
  }



template<typename eT>
inline
typename Mat<eT>::const_row_iterator&
Mat<eT>::const_row_iterator::operator++()
  {
  ++col;
  
  if(col >= M.n_cols)
    {
    col = 0;
    ++row;
    }
  
  return *this;
  }



template<typename eT>
inline
void
Mat<eT>::const_row_iterator::operator++(int)
  {
  operator++();
  }



template<typename eT>
inline
typename Mat<eT>::const_row_iterator&
Mat<eT>::const_row_iterator::operator--()
  {
  if(col > 0)
    {
    --col;
    }
  else
    {
    if(row > 0)
      {
      col = M.n_cols - 1;
      --row;
      }
    }
  
  return *this;
  }



template<typename eT>
inline
void
Mat<eT>::const_row_iterator::operator--(int)
  {
  operator--();
  }



template<typename eT>
inline
bool
Mat<eT>::const_row_iterator::operator!=(const typename Mat<eT>::const_row_iterator& X) const
  {
  return ( (row != X.row) || (col != X.col) ) ? true : false;
  }



template<typename eT>
inline
bool
Mat<eT>::const_row_iterator::operator==(const typename Mat<eT>::const_row_iterator& X) const
  {
  return ( (row == X.row) && (col == X.col) ) ? true : false;
  }



template<typename eT>
inline
typename Mat<eT>::iterator
Mat<eT>::begin()
  {
  arma_extra_debug_sigprint();
  
  return memptr();
  }



template<typename eT>
inline
typename Mat<eT>::const_iterator
Mat<eT>::begin() const
  {
  arma_extra_debug_sigprint();
  
  return memptr();
  }



template<typename eT>
inline
typename Mat<eT>::iterator
Mat<eT>::end()
  {
  arma_extra_debug_sigprint();
  
  return memptr() + n_elem;
  }



template<typename eT>
inline
typename Mat<eT>::const_iterator
Mat<eT>::end() const
  {
  arma_extra_debug_sigprint();
  
  return memptr() + n_elem;
  }
  


template<typename eT>
inline
typename Mat<eT>::col_iterator
Mat<eT>::begin_col(const u32 col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (col_num >= n_cols), "begin_col(): index out of bounds");
  
  return colptr(col_num);
  }



template<typename eT>
inline
typename Mat<eT>::const_col_iterator
Mat<eT>::begin_col(const u32 col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (col_num >= n_cols), "begin_col(): index out of bounds");
  
  return colptr(col_num);
  }



template<typename eT>
inline
typename Mat<eT>::col_iterator
Mat<eT>::end_col(const u32 col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (col_num >= n_cols), "end_col(): index out of bounds");
  
  return colptr(col_num) + n_rows;
  }



template<typename eT>
inline
typename Mat<eT>::const_col_iterator
Mat<eT>::end_col(const u32 col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (col_num >= n_cols), "end_col(): index out of bounds");
  
  return colptr(col_num) + n_rows;
  }
  


template<typename eT>
inline
typename Mat<eT>::row_iterator
Mat<eT>::begin_row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= n_rows), "Mat::begin_row(): index out of bounds" );
  
  return typename Mat<eT>::row_iterator(*this, row_num);
  }



template<typename eT>
inline
typename Mat<eT>::const_row_iterator
Mat<eT>::begin_row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= n_rows), "Mat::begin_row(): index out of bounds" );
  
  return typename Mat<eT>::const_row_iterator(*this, row_num);
  }



template<typename eT>
inline
typename Mat<eT>::row_iterator
Mat<eT>::end_row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= n_rows), "Mat::end_row(): index out of bounds" );
  
  return typename Mat<eT>::row_iterator(*this, row_num + 1);
  }



template<typename eT>
inline
typename Mat<eT>::const_row_iterator
Mat<eT>::end_row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (row_num >= n_rows), "Mat::end_row(): index out of bounds" );
  
  return typename Mat<eT>::const_row_iterator(*this, row_num + 1);
  }



//! prefix ++
template<typename eT>
arma_inline
void
Mat_aux::prefix_pp(Mat<eT>& x)
  {
        eT* memptr = x.memptr();
  const u32 n_elem = x.n_elem;
  
  u32 i,j;

  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    ++(memptr[i]);
    ++(memptr[j]);
    }
  
  if(i < n_elem)
    {
    ++(memptr[i]);
    }
  }



//! prefix ++ for complex numbers (work around for limitations of the std::complex class)
template<typename T>
arma_inline
void
Mat_aux::prefix_pp(Mat< std::complex<T> >& x)
  {
  x += T(1);
  }



//! postfix ++
template<typename eT>
arma_inline
void
Mat_aux::postfix_pp(Mat<eT>& x)
  {
        eT* memptr = x.memptr();
  const u32 n_elem = x.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    (memptr[i])++;
    (memptr[j])++;
    }
  
  if(i < n_elem)
    {
    (memptr[i])++;
    }
  }



//! postfix ++ for complex numbers (work around for limitations of the std::complex class)
template<typename T>
arma_inline
void
Mat_aux::postfix_pp(Mat< std::complex<T> >& x)
  {
  x += T(1);
  }



//! prefix --
template<typename eT>
arma_inline
void
Mat_aux::prefix_mm(Mat<eT>& x)
  {
        eT* memptr = x.memptr();
  const u32 n_elem = x.n_elem;

  u32 i,j;

  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    --(memptr[i]);
    --(memptr[j]);
    }
  
  if(i < n_elem)
    {
    --(memptr[i]);
    }
  }



//! prefix -- for complex numbers (work around for limitations of the std::complex class)
template<typename T>
arma_inline
void
Mat_aux::prefix_mm(Mat< std::complex<T> >& x)
  {
  x -= T(1);
  }



//! postfix --
template<typename eT>
arma_inline
void
Mat_aux::postfix_mm(Mat<eT>& x)
  {
        eT* memptr = x.memptr();
  const u32 n_elem = x.n_elem;

  u32 i,j;

  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    (memptr[i])--;
    (memptr[j])--;
    }
  
  if(i < n_elem)
    {
    (memptr[i])--;
    }
  }



//! postfix ++ for complex numbers (work around for limitations of the std::complex class)
template<typename T>
arma_inline
void
Mat_aux::postfix_mm(Mat< std::complex<T> >& x)
  {
  x -= T(1);
  }


#ifdef ARMA_EXTRA_MAT_MEAT
#include ARMA_EXTRA_MAT_MEAT
#endif


//! @}
