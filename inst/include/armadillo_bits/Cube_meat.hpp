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


//! \addtogroup Cube
//! @{


template<typename eT>
inline
Cube<eT>::~Cube()
  {
  arma_extra_debug_sigprint_this(this);
  
  delete_mat();
  
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
    access::rw(n_rows)       = 0;
    access::rw(n_cols)       = 0;
    access::rw(n_elem_slice) = 0;
    access::rw(n_slices)     = 0;
    access::rw(n_elem)       = 0;
    access::rw(mat_ptrs)     = 0;
    access::rw(mem)          = 0;
    }
  
  isnt_supported_elem_type<eT>::check();
  }



template<typename eT>
inline
Cube<eT>::Cube()
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  }



//! construct the cube to have user specified dimensions
template<typename eT>
inline
Cube<eT>::Cube(const u32 in_n_rows, const u32 in_n_cols, const u32 in_n_slices)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(in_n_rows, in_n_cols, in_n_slices);
  }



//! internal cube construction; if the requested size is small enough, memory from the stack is used.
//! otherwise memory is allocated via 'new'
template<typename eT>
inline
void
Cube<eT>::init(const u32 in_n_rows, const u32 in_n_cols, const u32 in_n_slices)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_n_rows = %d, in_n_cols = %d, in_n_slices = %d") % in_n_rows % in_n_cols % in_n_slices );
  
  const u32 new_n_elem = in_n_rows * in_n_cols * in_n_slices;

  if(n_elem == new_n_elem)
    {
    if( (n_rows != in_n_rows) || (n_cols != in_n_cols) || (n_slices != in_n_slices) )
      {
      delete_mat();

      access::rw(n_rows)       = in_n_rows;
      access::rw(n_cols)       = in_n_cols;
      access::rw(n_elem_slice) = in_n_rows*in_n_cols;
      access::rw(n_slices)     = in_n_slices;
    
      create_mat();
      }
    }
  else
    {
    arma_debug_check
      (
      (use_aux_mem == true),
      "Cube::init(): can't change the amount of memory as auxiliary memory is in use"
      );
      
    delete_mat();

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
      arma_check( (mem == 0), "Cube::init(): out of memory" );
      }
    
    access::rw(n_elem) = new_n_elem;

    if(new_n_elem == 0)
      {
      access::rw(n_rows)       = 0;
      access::rw(n_cols)       = 0;
      access::rw(n_elem_slice) = 0;
      access::rw(n_slices)     = 0;
      }
    else
      {
      access::rw(n_rows)       = in_n_rows;
      access::rw(n_cols)       = in_n_cols;
      access::rw(n_elem_slice) = in_n_rows*in_n_cols;
      access::rw(n_slices)     = in_n_slices;
      }
      
    create_mat();
    }
  }


template<typename eT>
inline
void
Cube<eT>::delete_mat()
  {
  arma_extra_debug_sigprint();
  
  for(u32 slice = 0; slice < n_slices; ++slice)
    {
    delete access::rw(mat_ptrs[slice]);
    }

  if(n_slices > sizeof(mat_ptrs_local)/sizeof(Mat<eT>*) )
    {
    delete [] mat_ptrs;
    }
  }



template<typename eT>
inline
void
Cube<eT>::create_mat()
  {
  arma_extra_debug_sigprint();

  if( n_slices <= sizeof(mat_ptrs_local)/sizeof(Mat<eT>*) )
    {
    access::rw(mat_ptrs) = const_cast< const Mat<eT>** >(mat_ptrs_local);
    }
  else
    {
    access::rw(mat_ptrs) = new(std::nothrow) const Mat<eT>*[n_slices];
    arma_check( (mat_ptrs == 0), "Cube::create_mat(): out of memory" );
    }
    
  for(u32 slice = 0; slice < n_slices; ++slice)
    {
    mat_ptrs[slice] = new Mat<eT>('j', slice_memptr(slice), n_rows, n_cols);
    }
  }



//! Set the cube to be equal to the specified scalar.
//! NOTE: the size of the cube will be 1x1x1
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  init(1,1,1);
  access::rw(mem[0]) = val;
  return *this;
  }



//! In-place addition of a scalar to all elements of the cube
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator+=(const eT val)
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



//! In-place subtraction of a scalar from all elements of the cube
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator-=(const eT val)
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



//! In-place multiplication of all elements of the cube with a scalar
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator*=(const eT val)
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



//! In-place division of all elements of the cube with a scalar
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator/=(const eT val)
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



//! construct a cube from a given cube
template<typename eT>
inline
Cube<eT>::Cube(const Cube<eT>& in_cube)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint(arma_boost::format("this = %x   in_cube = %x") % this % &in_cube);
  
  init(in_cube);
  }



//! construct a cube from a given cube
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator=(const Cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  init(x);
  return *this;
  }



//! construct a cube from a given cube
template<typename eT>
inline
void
Cube<eT>::init(const Cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    init(x.n_rows, x.n_cols, x.n_slices);
    syslib::copy_elem( memptr(), x.mem, n_elem );
    }
  }



//! construct a cube from a given auxiliary array of eTs.
//! if copy_aux_mem is true, new memory is allocated and the array is copied.
//! if copy_aux_mem is false, the auxiliary array is used directly (without allocating memory and copying).
//! note that in the latter case 
//! the default is to copy the array.

template<typename eT>
inline
Cube<eT>::Cube(eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const u32 aux_n_slices, const bool copy_aux_mem)
  : n_rows      (copy_aux_mem ? 0     : aux_n_rows                        )
  , n_cols      (copy_aux_mem ? 0     : aux_n_cols                        )
  , n_elem_slice(copy_aux_mem ? 0     : aux_n_rows*aux_n_cols             )
  , n_slices    (copy_aux_mem ? 0     : aux_n_slices                      )
  , n_elem      (copy_aux_mem ? 0     : aux_n_rows*aux_n_cols*aux_n_slices)
  , use_aux_mem (copy_aux_mem ? false : true                              )
  , mat_ptrs    (mat_ptrs                                                 )
  , mem         (copy_aux_mem ? mem   : aux_mem                           )
  {
  arma_extra_debug_sigprint_this(this);
  
  if(copy_aux_mem == true)
    {
    init(aux_n_rows, aux_n_cols, aux_n_slices);
    syslib::copy_elem( memptr(), aux_mem, n_elem );
    }
  else
    {
    create_mat();
    }
  }



//! construct a cube from a given auxiliary read-only array of eTs.
//! the array is copied.
template<typename eT>
inline
Cube<eT>::Cube(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const u32 aux_n_slices)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(aux_n_rows, aux_n_cols, aux_n_slices);
  syslib::copy_elem( memptr(), aux_mem, n_elem );
  }



//! in-place cube addition
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator+=(const Cube<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "cube addition");
  
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



//! in-place cube subtraction
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator-=(const Cube<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "cube subtraction");
  
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



//! in-place element-wise cube multiplication
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator%=(const Cube<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "element-wise cube multiplication");
  
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



//! in-place element-wise cube division
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator/=(const Cube<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "element-wise cube division");
  
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



//! for constructing a complex cube out of two non-complex cubes
template<typename eT>
template<typename T1, typename T2>
inline
Cube<eT>::Cube
  (
  const BaseCube<typename Cube<eT>::pod_type,T1>& A,
  const BaseCube<typename Cube<eT>::pod_type,T2>& B
  )
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  arma_type_check< is_complex<eT>::value == false >::apply();   //!< compile-time abort if eT isn't std::complex
  
  typedef typename T1::elem_type T;
  arma_type_check< is_complex<T>::value == true >::apply();   //!< compile-time abort if T is std::complex
  
  isnt_same_type<std::complex<T>, eT>::check();   //!< compile-time abort if types are not compatible
  
  const unwrap_cube<T1> tmp_A(A.get_ref());
  const unwrap_cube<T2> tmp_B(B.get_ref());
  
  const Cube<T>& X = tmp_A.M;
  const Cube<T>& Y = tmp_B.M;
  
  arma_assert_same_size(X, Y, "Cube()");
  
  init(X.n_rows, X.n_cols, X.n_slices);
  
  const T* X_mem = X.mem;
  const T* Y_mem = Y.mem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) = std::complex<T>(X_mem[i], Y_mem[i]);
    }
  }



//! construct a cube from a subview_cube instance (e.g. construct a cube from a delayed subcube operation)
template<typename eT>
inline
Cube<eT>::Cube(const subview_cube<eT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  this->operator=(X);
  }



//! construct a cube from a subview_cube instance (e.g. construct a cube from a delayed subcube operation)
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::extract(*this, X);
  return *this;
  }



//! in-place cube addition (using a subcube on the right-hand-side)
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator+=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::plus_inplace(*this, X);
  return *this;
  }



//! in-place cube subtraction (using a subcube on the right-hand-side)
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator-=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::minus_inplace(*this, X);
  return *this;
  }



//! in-place element-wise cube mutiplication (using a subcube on the right-hand-side)
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator%=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::schur_inplace(*this, X);
  return *this;
  }



//! in-place element-wise cube division (using a subcube on the right-hand-side)
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator/=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>::div_inplace(*this, X);
  return *this;
  }



//! provide the reference to the matrix representing a single slice
template<typename eT>
arma_inline
Mat<eT>&
Cube<eT>::slice(const u32 in_slice)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_slice >= n_slices),
    "Cube::slice(): index out of bounds"
    );
  
  return const_cast< Mat<eT>& >( *(mat_ptrs[in_slice]) );
  }



//! provide the reference to the matrix representing a single slice
template<typename eT>
arma_inline
const Mat<eT>&
Cube<eT>::slice(const u32 in_slice) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_slice >= n_slices),
    "Cube::slice(): index out of bounds"
    );
   
  return *(mat_ptrs[in_slice]);
  }



//! creation of subview_cube (subcube comprised of specified slices)
template<typename eT>
arma_inline
subview_cube<eT>
Cube<eT>::slices(const u32 in_slice1, const u32 in_slice2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_slice1 > in_slice2) || (in_slice2 >= n_slices),
    "Cube::slices(): indices out of bounds or incorrectly used"
    );
  
  return subview_cube<eT>(*this, 0, 0, in_slice1, ( (n_rows>0) ? n_rows-1 : 0 ), ( (n_cols>0) ? n_cols-1 : 0 ), in_slice2);
  }



//! creation of subview_cube (subcube comprised of specified slices)
template<typename eT>
arma_inline
const subview_cube<eT>
Cube<eT>::slices(const u32 in_slice1, const u32 in_slice2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_slice1 > in_slice2) || (in_slice2 >= n_slices),
    "Cube::rows(): indices out of bounds or incorrectly used"
    );
  
  return subview_cube<eT>(*this, 0, 0, in_slice1, ( (n_rows>0) ? n_rows-1 : 0 ), ( (n_cols>0) ? n_cols-1 : 0 ), in_slice2);
  }



//! creation of subview_cube (generic subcube)
template<typename eT>
arma_inline
subview_cube<eT>
Cube<eT>::subcube(const u32 in_row1, const u32 in_col1, const u32 in_slice1, const u32 in_row2, const u32 in_col2, const u32 in_slice2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 >  in_row2) || (in_col1 >  in_col2) || (in_slice1 >  in_slice2) ||
    (in_row2 >= n_rows)  || (in_col2 >= n_cols)  || (in_slice2 >= n_slices),
    "Cube::subcube(): indices out of bounds or incorrectly used"
    );
  
  return subview_cube<eT>(*this, in_row1, in_col1, in_slice1, in_row2, in_col2, in_slice2);
  }



//! creation of subview_cube (generic subcube)
template<typename eT>
arma_inline
const subview_cube<eT>
Cube<eT>::subcube(const u32 in_row1, const u32 in_col1, const u32 in_slice1, const u32 in_row2, const u32 in_col2, const u32 in_slice2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 >  in_row2) || (in_col1 >  in_col2) || (in_slice1 >  in_slice2) ||
    (in_row2 >= n_rows)  || (in_col2 >= n_cols)  || (in_slice2 >= n_slices),
    "Cube::subcube(): indices out of bounds or incorrectly used"
    );
    
  return subview_cube<eT>(*this, in_row1, in_col1, in_slice1, in_row2, in_col2, in_slice2);
  }



//! create a cube from OpCube, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename op_type>
inline
Cube<eT>::Cube(const OpCube<T1, op_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  op_type::apply(*this, X);
  }



//! create a cube from OpCube, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator=(const OpCube<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  op_type::apply(*this, X);
  
  return *this;
  }



//! in-place cube addition, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator+=(const OpCube<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const Cube<eT> m(X);
  
  return (*this).operator+=(m);
  }



//! in-place cube subtraction, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator-=(const OpCube<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const Cube<eT> m(X);
  
  return (*this).operator-=(m);
  }



//! in-place cube element-wise multiplication, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator%=(const OpCube<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const Cube<eT> m(X);
  
  return (*this).operator%=(m);
  }



//! in-place cube element-wise division, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator/=(const OpCube<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const Cube<eT> m(X);
  
  return (*this).operator/=(m);
  }



//! create a cube from eOpCube, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename eop_type>
inline
Cube<eT>::Cube(const eOpCube<T1, eop_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply(*this, X);
  }



//! create a cube from eOpCube, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename eop_type>
inline
const Cube<eT>&
Cube<eT>::operator=(const eOpCube<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply(*this, X);
  
  return *this;
  }



//! in-place cube addition, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename eop_type>
inline
const Cube<eT>&
Cube<eT>::operator+=(const eOpCube<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply_inplace_plus(*this, X);
  
  return *this;
  }



//! in-place cube subtraction, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename eop_type>
inline
const Cube<eT>&
Cube<eT>::operator-=(const eOpCube<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  eop_type::apply_inplace_minus(*this, X);
  
  return *this;
  }



//! in-place cube element-wise multiplication, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename eop_type>
inline
const Cube<eT>&
Cube<eT>::operator%=(const eOpCube<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();

  eop_type::apply_inplace_schur(*this, X);
  
  return *this;
  }



//! in-place cube element-wise division, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename eop_type>
inline
const Cube<eT>&
Cube<eT>::operator/=(const eOpCube<T1, eop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();

  eop_type::apply_inplace_div(*this, X);
  
  return *this;
  }



//! create a cube from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Cube<eT>::Cube(const GlueCube<T1, T2, glue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  this->operator=(X);
  }



//! create a cube from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator=(const GlueCube<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_type::apply(*this, X);
  
  return *this;
  }


//! in-place cube addition, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator+=(const GlueCube<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  const Cube<eT> m(X);
  
  return (*this).operator+=(m);
  }



//! in-place cube subtraction, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator-=(const GlueCube<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  const Cube<eT> m(X);
  
  return (*this).operator-=(m);
  }



//! in-place cube element-wise multiplication, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator%=(const GlueCube<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  const Cube<eT> m(X);
  
  return (*this).operator%=(m);
  }



//! in-place cube element-wise division, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator/=(const GlueCube<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  const Cube<eT> m(X);
  
  return (*this).operator/=(m);
  }



//! create a cube from eGlue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
Cube<eT>::Cube(const eGlueCube<T1, T2, eglue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , use_aux_mem(false)
  , mat_ptrs(mat_ptrs)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  this->operator=(X);
  }



//! create a cube from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Cube<eT>&
Cube<eT>::operator=(const eGlueCube<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply(*this, X);
  
  return *this;
  }


//! in-place cube addition, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Cube<eT>&
Cube<eT>::operator+=(const eGlueCube<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply_inplace_plus(*this, X);
  
  return *this;
  }



//! in-place cube subtraction, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Cube<eT>&
Cube<eT>::operator-=(const eGlueCube<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply_inplace_minus(*this, X);
  
  return *this;
  }



//! in-place cube element-wise multiplication, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Cube<eT>&
Cube<eT>::operator%=(const eGlueCube<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply_inplace_schur(*this, X);
  
  return *this;
  }



//! in-place cube element-wise division, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Cube<eT>&
Cube<eT>::operator/=(const eGlueCube<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  eglue_type::apply_inplace_div(*this, X);
  
  return *this;
  }



//! linear element accessor (treats the cube as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT&
Cube<eT>::operator() (const u32 i)
  {
  arma_debug_check( (i >= n_elem), "Cube::operator(): index out of bounds");
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the cube as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT
Cube<eT>::operator() (const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "Cube::operator(): index out of bounds");
  return mem[i];
  }


//! linear element accessor (treats the cube as a vector); no bounds check.  
template<typename eT>
arma_inline
eT&
Cube<eT>::operator[] (const u32 i)
  {
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the cube as a vector); no bounds check
template<typename eT>
arma_inline
eT
Cube<eT>::operator[] (const u32 i) const
  {
  return mem[i];
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT&
Cube<eT>::operator() (const u32 in_row, const u32 in_col, const u32 in_slice)
  {
  arma_debug_check
    (
    (in_row >= n_rows) ||
    (in_col >= n_cols) ||
    (in_slice >= n_slices)
    ,
    "Cube::operator(): index out of bounds"
    );

  return access::rw(mem[in_slice*n_elem_slice + in_col*n_rows + in_row]);
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT
Cube<eT>::operator() (const u32 in_row, const u32 in_col, const u32 in_slice) const
  {
  arma_debug_check
    (
    (in_row >= n_rows) ||
    (in_col >= n_cols) ||
    (in_slice >= n_slices)
    ,
    "Cube::operator(): index out of bounds"
    );

  return mem[in_slice*n_elem_slice + in_col*n_rows + in_row];
  }



//! element accessor; no bounds check
template<typename eT>
arma_inline
eT&
Cube<eT>::at(const u32 in_row, const u32 in_col, const u32 in_slice)
  {
  return access::rw( mem[in_slice*n_elem_slice + in_col*n_rows + in_row] );
  }



//! element accessor; no bounds check
template<typename eT>
arma_inline
eT
Cube<eT>::at(const u32 in_row, const u32 in_col, const u32 in_slice) const
  {
  return mem[in_slice*n_elem_slice + in_col*n_rows + in_row];
  }



//! prefix ++
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator++()
  {
  Cube_aux::prefix_pp(*this);
  return *this;
  }



//! postfix ++  (must not return the object by reference)
template<typename eT>
arma_inline
void
Cube<eT>::operator++(int)
  {
  Cube_aux::postfix_pp(*this);
  }



//! prefix --
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator--()
  {
  Cube_aux::prefix_mm(*this);
  return *this;
  }



//! postfix --  (must not return the object by reference)
template<typename eT>
arma_inline
void
Cube<eT>::operator--(int)
  {
  Cube_aux::postfix_mm(*this);
  }



//! returns true if all of the elements are finite
template<typename eT>
arma_inline
bool
Cube<eT>::is_finite() const
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



//! returns a pointer to array of eTs used by the cube
template<typename eT>
arma_inline
eT*
Cube<eT>::memptr()
  {
  return const_cast<eT*>(mem);
  }



//! returns a pointer to array of eTs used by the cube
template<typename eT>
arma_inline
const eT*
Cube<eT>::memptr() const
  {
  return mem;
  }



//! returns a pointer to array of eTs used by the specified slice in the cube
template<typename eT>
arma_inline
eT*
Cube<eT>::slice_memptr(const u32 slice)
  {
  return const_cast<eT*>( &mem[ slice*n_elem_slice ] );
  }



//! returns a pointer to array of eTs used by the specified slice in the cube
template<typename eT>
arma_inline
const eT*
Cube<eT>::slice_memptr(const u32 slice) const
  {
  return &mem[ slice*n_elem_slice ];
  }



//! returns a pointer to array of eTs used by the specified slice in the cube
template<typename eT>
arma_inline
eT*
Cube<eT>::slice_colptr(const u32 slice, const u32 col)
  {
  return const_cast<eT*>( &mem[ slice*n_elem_slice + col*n_rows] );
  }



//! returns a pointer to array of eTs used by the specified slice in the cube
template<typename eT>
arma_inline
const eT*
Cube<eT>::slice_colptr(const u32 slice, const u32 col) const
  {
  return &mem[ slice*n_elem_slice + col*n_rows ];
  }



//! print contents of the cube (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's flags are restored to their original values.
template<typename eT>
inline
void
Cube<eT>::print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    cout << extra_text << '\n';
    }
  
  arma_ostream::print(cout, *this, true);
  }


//! print contents of the cube to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's flags are restored to their original values.
template<typename eT>
inline
void
Cube<eT>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    user_stream << extra_text << '\n';
    }
  
  arma_ostream::print(user_stream, *this, true);
  }



//! print contents of the cube (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the stream's flags are used as is and are not modified
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Cube<eT>::raw_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    cout << extra_text << '\n';
    }
  
  arma_ostream::print(cout, *this, false);
  }



//! print contents of the cube to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the stream's flags are used as is and are not modified.
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Cube<eT>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    user_stream << extra_text << '\n';
    }
  
  arma_ostream::print(user_stream, *this, false);
  }



//! change the cube to have user specified dimensions (data is not preserved)
template<typename eT>
inline
void
Cube<eT>::set_size(const u32 in_n_rows, const u32 in_n_cols, const u32 in_n_slices)
  {
  arma_extra_debug_sigprint();
  
  init(in_n_rows, in_n_cols, in_n_slices);
  }



//! change the cube (without preserving data) to have the same dimensions as the given cube 
template<typename eT>
template<typename eT2>
inline
void
Cube<eT>::copy_size(const Cube<eT2>& m)
  {
  arma_extra_debug_sigprint();
  
  init(m.n_rows, m.n_cols, m.n_slices);
  }



//! fill the cube with the specified value
template<typename eT>
inline
void
Cube<eT>::fill(const eT val)
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
Cube<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(0));
  }



template<typename eT>
inline
void
Cube<eT>::zeros(const u32 in_rows, const u32 in_cols, const u32 in_slices)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_rows = %d, in_cols = %d, in_slices = %d") % in_rows % in_cols % in_slices );

  set_size(in_rows, in_cols, in_slices);
  fill(eT(0));
  }



template<typename eT>
inline
void
Cube<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(1));
  }



template<typename eT>
inline
void
Cube<eT>::ones(const u32 in_rows, const u32 in_cols, const u32 in_slices)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_rows = %d, in_cols = %d, in_slices = %d") % in_rows % in_cols % in_slices );

  set_size(in_rows, in_cols, in_slices);
  fill(eT(1));
  }



template<typename eT>
inline
void
Cube<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  init(0,0,0);
  }



//! save the cube to a file
template<typename eT>
inline
void
Cube<eT>::save(const std::string name, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case raw_ascii:
      diskio::save_raw_ascii(*this, name);
      break;
    
    case arma_ascii:
      diskio::save_arma_ascii(*this, name);
      break;
    
    case arma_binary:
      diskio::save_arma_binary(*this, name);
      break;
      
    case ppm_binary:
      diskio::save_ppm_binary(*this, name);
      break;

    default:
      arma_stop("Cube::save(): unsupported file type");
    }
  
  }



//! save the cube to a stream
template<typename eT>
inline
void
Cube<eT>::save(std::ostream& os, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case raw_ascii:
      diskio::save_raw_ascii(*this, "[ostream]", os);
      break;
    
    case arma_ascii:
      diskio::save_arma_ascii(*this, "[ostream]", os);
      break;
    
    case arma_binary:
      diskio::save_arma_binary(*this, "[ostream]", os);
      break;
      
    case ppm_binary:
      diskio::save_ppm_binary(*this, "[ostream]", os);
      break;

    default:
      arma_stop("Cube::save(): unsupported file type");
    }
  
  }



//! load a cube from a file
template<typename eT>
inline
void
Cube<eT>::load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      diskio::load_auto_detect(*this, name);
      break;
    
    case raw_ascii:
      diskio::load_raw_ascii(*this, name);
      break;
    
    case arma_ascii:
      diskio::load_arma_ascii(*this, name);
      break;
    
    case arma_binary:
      diskio::load_arma_binary(*this, name);
      break;
      
    case ppm_binary:
      diskio::load_ppm_binary(*this, name);
      break;

    default:
      arma_stop("Cube::load(): unsupported file type");
    }
  
  }



//! load a cube from a stream
template<typename eT>
inline
void
Cube<eT>::load(std::istream& is, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      diskio::load_auto_detect(*this, "[istream]", is);
      break;
    
    case raw_ascii:
      diskio::load_raw_ascii(*this, "[istream]", is);
      break;
    
    case arma_ascii:
      diskio::load_arma_ascii(*this, "[istream]", is);
      break;
    
    case arma_binary:
      diskio::load_arma_binary(*this, "[istream]", is);
      break;
      
    case ppm_binary:
      diskio::load_ppm_binary(*this, "[istream]", is);
      break;
    
    default:
      arma_stop("Cube::load(): unsupported file type");
    }
  
  }



//! prefix ++
template<typename eT>
arma_inline
void
Cube_aux::prefix_pp(Cube<eT>& x)
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
Cube_aux::prefix_pp(Cube< std::complex<T> >& x)
  {
  x += T(1);
  }



//! postfix ++
template<typename eT>
arma_inline
void
Cube_aux::postfix_pp(Cube<eT>& x)
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
Cube_aux::postfix_pp(Cube< std::complex<T> >& x)
  {
  x += T(1);
  }



//! prefix --
template<typename eT>
arma_inline
void
Cube_aux::prefix_mm(Cube<eT>& x)
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
Cube_aux::prefix_mm(Cube< std::complex<T> >& x)
  {
  x -= T(1);
  }



//! postfix --
template<typename eT>
arma_inline
void
Cube_aux::postfix_mm(Cube<eT>& x)
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
Cube_aux::postfix_mm(Cube< std::complex<T> >& x)
  {
  x -= T(1);
  }



//! @}
