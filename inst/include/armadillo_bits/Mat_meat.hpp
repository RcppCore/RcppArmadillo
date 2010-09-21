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
  
  if(mem_state == 0)
    {
    if(n_elem > Mat_prealloc::mem_n_elem)
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
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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
  
  const bool same_size = ( (n_rows == in_n_rows) && (n_cols == in_n_cols) );
  
  if(same_size == false)
    {
    const u32 t_vec_state = vec_state;
    const u32 t_mem_state = mem_state;
    
    arma_debug_check
      (
      (t_mem_state == 3),
      "Mat::init(): size is fixed and hence cannot be changed"
      );
    
    arma_debug_check
      (
        (
        (t_vec_state > 0) &&
          (
          ((t_vec_state == 1) && (in_n_cols > 1)) ||
          ((t_vec_state == 2) && (in_n_rows > 1))
          )
        ),
      "Mat::init(): object is a row or column vector; requested size is not compatible"
      );
    
    const u32 old_n_elem = n_elem;
    const u32 new_n_elem = in_n_rows * in_n_cols;
    
    if(old_n_elem == new_n_elem)
      {
      arma_extra_debug_print("Mat::init(): reusing memory");
      
      access::rw(n_rows) = in_n_rows;
      access::rw(n_cols) = in_n_cols;
      }
    else
      {
      arma_debug_check
        (
        (t_mem_state == 2),
        "Mat::init(): requested size is not compatible with the size of auxiliary memory"
        );
      
      if(t_mem_state == 0)
        {
        if(old_n_elem > Mat_prealloc::mem_n_elem )
          {
          arma_extra_debug_print("Mat::init(): freeing memory");
          
          delete [] mem;
          }
        }
      
      
      if(new_n_elem <= Mat_prealloc::mem_n_elem )
        {
        access::rw(mem) = mem_local;
        }
      else
        {
        arma_extra_debug_print("Mat::init(): allocating memory");
        
        access::rw(mem) = new(std::nothrow) eT[new_n_elem];
        
        arma_check( (mem == 0), "Mat::init(): out of memory" );
        }
      
      access::rw(n_elem)    = new_n_elem;
      access::rw(n_rows)    = in_n_rows;
      access::rw(n_cols)    = in_n_cols;
      access::rw(mem_state) = 0;
      }
    
    
    if(new_n_elem == 0)
      {
      access::rw(n_rows) = 0;
      access::rw(n_cols) = 0;
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
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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
  
  arrayops::inplace_plus( memptr(), val, n_elem );
  
  return *this;
  }



//! In-place subtraction of a scalar from all elements of the matrix
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_minus( memptr(), val, n_elem );
  
  return *this;
  }



//! In-place multiplication of all elements of the matrix with a scalar
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_mul( memptr(), val, n_elem );
  
  return *this;
  }



//! In-place division of all elements of the matrix with a scalar
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_div( memptr(), val, n_elem );
  
  return *this;
  }



//! construct a matrix from a given matrix
template<typename eT>
inline
Mat<eT>::Mat(const Mat<eT>& in_mat)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , vec_state(0)
  , mem_state(0)
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
    syslib::copy_elem( memptr(), x.mem, x.n_elem );
    }
  }



//! for constructing a complex matrix out of two non-complex matrices
template<typename eT>
template<typename T1, typename T2>
inline
void
Mat<eT>::init
  (
  const Base<typename Mat<eT>::pod_type, T1>& A,
  const Base<typename Mat<eT>::pod_type, T2>& B
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type      T;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  arma_type_check< is_complex<eT>::value == false >::apply();   //!< compile-time abort if eT isn't std::complex
  arma_type_check< is_complex< T>::value == true  >::apply();   //!< compile-time abort if T is std::complex
  
  isnt_same_type<std::complex<T>, eT>::check();   //!< compile-time abort if types are not compatible
  
  const Proxy<T1> X(A.get_ref());
  const Proxy<T2> Y(B.get_ref());
  
  arma_assert_same_size(X, Y, "Mat()");
  
  init(X.get_n_rows(), X.get_n_cols());
  
  const u32      N       = n_elem;
        eT*      out_mem = memptr();
        ea_type1 PX      = X.get_ea();
        ea_type2 PY      = Y.get_ea();
  
  for(u32 i=0; i<N; ++i)
    {
    out_mem[i] = std::complex<T>(PX[i], PY[i]);
    }
  }



//! try to steal the memory from a given matrix; 
//! if memory can't be stolen, copy the given matrix
template<typename eT>
inline
void
Mat<eT>::steal_mem(Mat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    const u32 x_n_rows    = x.n_rows;
    const u32 x_n_cols    = x.n_cols;
    const u32 x_n_elem    = x.n_elem;
    const u32 x_vec_state = x.vec_state;
    const u32 x_mem_state = x.mem_state;
    
    const u32 t_vec_state = vec_state;
    
    bool layout_ok = false;
    
    if(t_vec_state == x_vec_state)
      {
      layout_ok = true;
      }
    else
      {
      if( (t_vec_state == 1) && ( x_n_cols <= 1) )
        {
        layout_ok = true;
        }
      
      if( (t_vec_state == 2) && ( x_n_rows <= 1) )
        {
        layout_ok = true;
        }
      }
    
    
    if( (x_mem_state == 0) && (x_n_elem > Mat_prealloc::mem_n_elem) && (layout_ok == true) )
      {
      reset();
      
      access::rw(n_rows) = x_n_rows;
      access::rw(n_cols) = x_n_cols;
      access::rw(n_elem) = x_n_elem;
      access::rw(mem)    = x.mem;
      
      access::rw(x.n_rows) = 0;
      access::rw(x.n_cols) = 0;
      access::rw(x.n_elem) = 0;
      access::rw(x.mem)    = 0;
      }
    else
      {
      init(x);
      }
    }
  }



//! construct a matrix from a given auxiliary array of eTs.
//! if copy_aux_mem is true, new memory is allocated and the array is copied.
//! if copy_aux_mem is false, the auxiliary array is used directly (without allocating memory and copying).
//! the default is to copy the array.

template<typename eT>
inline
Mat<eT>::Mat(eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const bool copy_aux_mem, const bool strict)
  : n_rows   (copy_aux_mem ? 0   : aux_n_rows           )
  , n_cols   (copy_aux_mem ? 0   : aux_n_cols           )
  , n_elem   (copy_aux_mem ? 0   : aux_n_rows*aux_n_cols)
  , vec_state(               0                          )
  , mem_state(copy_aux_mem ? 0   : ( strict ? 2 : 1 )   )
  , mem      (copy_aux_mem ? mem : aux_mem              )
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
  , vec_state(0)
  , mem_state(0)
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
  : n_rows   (aux_n_rows           )
  , n_cols   (aux_n_cols           )
  , n_elem   (aux_n_rows*aux_n_cols)
  , vec_state(0                    )
  , mem_state(3                    )
  , mem      (aux_mem              )
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
  
  arrayops::inplace_plus( memptr(), m.memptr(), n_elem );
  
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
  
  arrayops::inplace_minus( memptr(), m.memptr(), n_elem );
  
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
  
  arrayops::inplace_mul( memptr(), m.memptr(), n_elem );
  
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
  
  arrayops::inplace_div( memptr(), m.memptr(), n_elem );
  
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
  , vec_state(0)
  , mem_state(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(A,B);
  }



//! construct a matrix from subview (e.g. construct a matrix from a delayed submatrix operation)
template<typename eT>
inline
Mat<eT>::Mat(const subview<eT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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



//! creation of subview (submatrix)
template<typename eT>
arma_inline
subview<eT>
Mat<eT>::submat(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const u32 in_row1 = row_span.a;
  const u32 in_row2 = row_span.b;
  
  const u32 in_col1 = col_span.a;
  const u32 in_col2 = col_span.b;
  
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
Mat<eT>::submat(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const u32 in_row1 = row_span.a;
  const u32 in_row2 = row_span.b;
  
  const u32 in_col1 = col_span.a;
  const u32 in_col2 = col_span.b;
  
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



//! remove specified row
template<typename eT>
inline
void
Mat<eT>::shed_row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= n_rows, "Mat::shed_row(): out of bounds");
  
  shed_rows(row_num, row_num);
  }



//! remove specified column
template<typename eT>
inline
void
Mat<eT>::shed_col(const u32 col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "Mat::shed_col(): out of bounds");
  
  shed_cols(col_num, col_num);
  }



//! remove specified rows
template<typename eT>
inline
void
Mat<eT>::shed_rows(const u32 in_row1, const u32 in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "Mat::shed_rows(): indices out of bounds or incorrectly used"
    );
  
  const u32 n_keep_front = in_row1;
  const u32 n_keep_back  = n_rows - (in_row2 + 1);
  
  Mat<eT> X(n_keep_front + n_keep_back, n_cols);
  
  if(n_keep_front > 0)
    {
    X.rows( 0, (n_keep_front-1) ) = rows( 0, (in_row1-1) );
    }
  
  if(n_keep_back > 0)
    {
    X.rows( n_keep_front,  (n_keep_front+n_keep_back-1) ) = rows( (in_row2+1), (n_rows-1) );
    }
  
  steal_mem(X);
  }



//! remove specified columns
template<typename eT>
inline
void
Mat<eT>::shed_cols(const u32 in_col1, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "Mat::shed_cols(): indices out of bounds or incorrectly used"
    );
  
  const u32 n_keep_front = in_col1;
  const u32 n_keep_back  = n_cols - (in_col2 + 1);
  
  Mat<eT> X(n_rows, n_keep_front + n_keep_back);
  
  if(n_keep_front > 0)
    {
    X.cols( 0, (n_keep_front-1) ) = cols( 0, (in_col1-1) );
    }
  
  if(n_keep_back > 0)
    {
    X.cols( n_keep_front,  (n_keep_front+n_keep_back-1) ) = cols( (in_col2+1), (n_cols-1) );
    }
  
  steal_mem(X);
  }



//! insert N rows at the specified row position,
//! optionally setting the elements of the inserted rows to zero
template<typename eT>
inline
void
Mat<eT>::insert_rows(const u32 row_num, const u32 N, const bool set_to_zero)
  {
  arma_extra_debug_sigprint();
  
  const u32 t_n_rows = n_rows;
  const u32 t_n_cols = n_cols;
  
  const u32 A_n_rows = row_num;
  const u32 B_n_rows = t_n_rows - row_num;
  
  // insertion at row_num == n_rows is in effect an append operation
  arma_debug_check( (row_num > t_n_rows), "Mat::insert_rows(): out of bounds");
  
  if(N > 0)
    {
    Mat<eT> out(t_n_rows + N, t_n_cols);
    
    if(A_n_rows > 0)
      {
      out.rows(0, A_n_rows-1) = rows(0, A_n_rows-1);
      }
    
    if(B_n_rows > 0)
      {
      out.rows(row_num + N, t_n_rows + N - 1) = rows(row_num, t_n_rows-1);
      }
    
    if(set_to_zero == true)
      {
      out.rows(row_num, row_num + N - 1).zeros();
      }
    
    steal_mem(out);
    }
  }



//! insert N columns at the specified column position,
//! optionally setting the elements of the inserted columns to zero
template<typename eT>
inline
void
Mat<eT>::insert_cols(const u32 col_num, const u32 N, const bool set_to_zero)
  {
  arma_extra_debug_sigprint();
  
  const u32 t_n_rows = n_rows;
  const u32 t_n_cols = n_cols;
  
  const u32 A_n_cols = col_num;
  const u32 B_n_cols = t_n_cols - col_num;
  
  // insertion at col_num == n_cols is in effect an append operation
  arma_debug_check( (col_num > t_n_cols), "Mat::insert_cols(): out of bounds");
  
  if(N > 0)
    {
    Mat<eT> out(t_n_rows, t_n_cols + N);
    
    if(A_n_cols > 0)
      {
      out.cols(0, A_n_cols-1) = cols(0, A_n_cols-1);
      }
    
    if(B_n_cols > 0)
      {
      out.cols(col_num + N, t_n_cols + N - 1) = cols(col_num, t_n_cols-1);
      }
    
    if(set_to_zero == true)
      {
      out.cols(col_num, col_num + N - 1).zeros();
      }
    
    steal_mem(out);
    }
  }



//! insert the given object at the specified row position; 
//! the given object must have the same number of columns as the matrix
template<typename eT>
template<typename T1>
inline
void
Mat<eT>::insert_rows(const u32 row_num, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& C = tmp.M;
  
  const u32 N = C.n_rows;
  
  const u32 t_n_rows = n_rows;
  const u32 t_n_cols = n_cols;
  
  const u32 A_n_rows = row_num;
  const u32 B_n_rows = t_n_rows - row_num;

  // insertion at row_num == n_rows is in effect an append operation
  arma_debug_check( (row_num  >  t_n_rows), "Mat::insert_rows(): out of bounds");
  arma_debug_check( (C.n_cols != t_n_cols), "Mat::insert_rows(): given object has an incompatible number of columns");
  
  if(N > 0)
    {
    Mat<eT> out(t_n_rows + N, t_n_cols);
    
    if(A_n_rows > 0)
      {
      out.rows(0, A_n_rows-1) = rows(0, A_n_rows-1);
      }
    
    if(B_n_rows > 0)
      {
      out.rows(row_num + N, t_n_rows + N - 1) = rows(row_num, t_n_rows - 1);
      }
    
    out.rows(row_num, row_num + N - 1) = C;
    
    steal_mem(out);
    }
  }



//! insert the given object at the specified column position; 
//! the given object must have the same number of rows as the matrix
template<typename eT>
template<typename T1>
inline
void
Mat<eT>::insert_cols(const u32 col_num, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& C = tmp.M;
  
  const u32 N = C.n_cols;
  
  const u32 t_n_rows = n_rows;
  const u32 t_n_cols = n_cols;
  
  const u32 A_n_cols = col_num;
  const u32 B_n_cols = t_n_cols - col_num;

  // insertion at col_num == n_cols is in effect an append operation
  arma_debug_check( (col_num  >  t_n_cols), "Mat::insert_cols(): out of bounds");
  arma_debug_check( (C.n_rows != t_n_rows), "Mat::insert_cols(): given object has an incompatible number of rows");
  
  if(N > 0)
    {
    Mat<eT> out(t_n_rows, t_n_cols + N);
    
    if(A_n_cols > 0)
      {
      out.cols(0, A_n_cols-1) = cols(0, A_n_cols-1);
      }
    
    if(B_n_cols > 0)
      {
      out.cols(col_num + N, t_n_cols + N - 1) = cols(col_num, t_n_cols - 1);
      }
    
    out.cols(col_num, col_num + N - 1) = C;
    
    steal_mem(out);
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
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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
  , vec_state(0)
  , mem_state(0)
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



//! returns true if the matrix has no elements
template<typename eT>
arma_inline
bool
Mat<eT>::is_empty() const
  {
  return (n_elem == 0);
  }



//! returns true if the given index is currently in range
template<typename eT>
arma_inline
bool
Mat<eT>::in_range(const u32 i) const
  {
  return (i < n_elem);
  }



//! returns true if the given location is currently in range
template<typename eT>
arma_inline
bool
Mat<eT>::in_range(const u32 in_row, const u32 in_col) const
  {
  return ( (in_row < n_rows) && (in_col < n_cols) );
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
//! on return, the stream's state are restored to their original values.
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
//! on return, the stream's state are restored to their original values.
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
//! on return, the stream's state are restored to their original values.
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
//! on return, the stream's state are restored to their original values.
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
//! the stream's state are used as is and are not modified
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
//! the stream's state are used as is and are not modified.
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
//! the stream's state are used as is and are not modified
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
//! the stream's state are used as is and are not modified.
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
Mat<eT>::set_size(const u32 in_elem)
  {
  arma_extra_debug_sigprint();
  
  switch(vec_state)
    {
    case 0:
    case 1:
      init(in_elem, 1);
      break;
    
    case 2:
      init(1, in_elem);
      break;
      
    default:
      ;
    }
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
const Mat<eT>&
Mat<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set( memptr(), val, n_elem );
  
  return *this;
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  return fill(eT(0));
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::zeros(const u32 in_elem)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_elem);
  
  return fill(eT(0));
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::zeros(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint();

  set_size(in_rows, in_cols);
  
  return fill(eT(0));
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  return fill(eT(1));
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::ones(const u32 in_elem)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_elem);
  
  return fill(eT(1));
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::ones(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint();

  set_size(in_rows, in_cols);
  
  return fill(eT(1));
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::randu()
  {
  arma_extra_debug_sigprint();
  
  const u32 N   = n_elem;
        eT* ptr = memptr();
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    ptr[i] = eT(eop_aux_randu<eT>());
    ptr[j] = eT(eop_aux_randu<eT>());
    }
  
  if(i < N)
    {
    ptr[i] = eT(eop_aux_randu<eT>());
    }
  
  return *this;
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::randu(const u32 in_elem)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_elem);
  
  return (*this).randu();
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::randu(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_rows, in_cols);
  
  return (*this).randu();
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::randn()
  {
  arma_extra_debug_sigprint();
  
  const u32 N   = n_elem;
        eT* ptr = memptr();
  
  for(u32 i=0; i<N; ++i)
    {
    ptr[i] = eT(eop_aux_randn<eT>());
    }
  
  return *this;
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::randn(const u32 in_elem)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_elem);
  
  return (*this).randn();
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::randn(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_rows, in_cols);
  
  return (*this).randn();
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::eye()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(0));
  
  const u32 N = (std::min)(n_rows, n_cols);
  
  for(u32 i=0; i<N; ++i)
    {
    at(i,i) = eT(1);
    }
  
  return *this;
  }



template<typename eT>
inline
const Mat<eT>&
Mat<eT>::eye(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_rows, in_cols);
  
  return (*this).eye();
  }



template<typename eT>
inline
void
Mat<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  init(0,0);
  }



template<typename eT>
template<typename T1>
inline
void
Mat<eT>::set_real(const Base<typename Mat<eT>::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat_aux::set_real(*this, X);
  }



template<typename eT>
template<typename T1>
inline
void
Mat<eT>::set_imag(const Base<typename Mat<eT>::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat_aux::set_imag(*this, X);
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



template<typename eT>
template<u32 fixed_n_rows, u32 fixed_n_cols>
arma_inline
void
Mat<eT>::fixed<fixed_n_rows, fixed_n_cols>::mem_setup()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(fixed_n_elem > 0)
    {
    access::rw(Mat<eT>::n_rows)    = fixed_n_rows;
    access::rw(Mat<eT>::n_cols)    = fixed_n_cols;
    access::rw(Mat<eT>::n_elem)    = fixed_n_elem;
    access::rw(Mat<eT>::vec_state) = 0;
    access::rw(Mat<eT>::mem_state) = 3;
    access::rw(Mat<eT>::mem)       = (fixed_n_elem > Mat_prealloc::mem_n_elem) ? mem_local_extra : mem_local;
    }
  else
    {
    access::rw(Mat<eT>::n_rows)    = 0;
    access::rw(Mat<eT>::n_cols)    = 0;
    access::rw(Mat<eT>::n_elem)    = 0;
    access::rw(Mat<eT>::vec_state) = 0;
    access::rw(Mat<eT>::mem_state) = 3;
    access::rw(Mat<eT>::mem)       = 0;
    }
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



template<typename eT, typename T1>
inline
void
Mat_aux::set_real(Mat<eT>& out, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  arma_debug_assert_same_size( out, A, "Mat::set_real()" );
  
  out = A;
  }



template<typename eT, typename T1>
inline
void
Mat_aux::set_imag(Mat<eT>& out, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  }



template<typename T, typename T1>
inline
void
Mat_aux::set_real(Mat< std::complex<T> >& out, const Base<T,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T>    eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_assert_same_size( out, A, "Mat::set_real()" );
  
  const u32     n_elem  = out.n_elem;
        eT*     out_mem = out.memptr();
        ea_type PA      = A.get_ea();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    //out_mem[i].real() = PA[i];
    out_mem[i] = std::complex<T>( PA[i], out_mem[i].imag() );
    }
  }



template<typename T, typename T1>
inline
void
Mat_aux::set_imag(Mat< std::complex<T> >& out, const Base<T,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T>    eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> A(X.get_ref());
  
  arma_debug_assert_same_size( out, A, "Mat::set_imag()" );
  
  const u32     n_elem  = out.n_elem;
        eT*     out_mem = out.memptr();
        ea_type PA      = A.get_ea();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    //out_mem[i].imag() = PA[i];
    out_mem[i] = std::complex<T>( out_mem[i].real(), PA[i] );
    }
  }



#ifdef ARMA_EXTRA_MAT_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_MAT_MEAT)
#endif



//! @}
