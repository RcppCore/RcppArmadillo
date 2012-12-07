// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
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
  
  if(mem_state == 0)
    {
    if(n_elem > Cube_prealloc::mem_n_elem)
      {
      memory::release( access::rw(mem) );
      }
    }
  
  if(arma_config::debug == true)
    {
    // try to expose buggy user code that accesses deleted objects
    access::rw(mat_ptrs) = 0;
    access::rw(mem)      = 0;
    }
  
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  }



template<typename eT>
inline
Cube<eT>::Cube()
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  }



//! construct the cube to have user specified dimensions
template<typename eT>
inline
Cube<eT>::Cube(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices)
  : n_rows(in_n_rows)
  , n_cols(in_n_cols)
  , n_elem_slice(in_n_rows*in_n_cols)
  , n_slices(in_n_slices)
  , n_elem(in_n_rows*in_n_cols*in_n_slices)
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold();
  }



template<typename eT>
inline
void
Cube<eT>::init_cold()
  {
  arma_extra_debug_sigprint( arma_boost::format("n_rows = %d, n_cols = %d, n_slices = %d") % n_rows % n_cols % n_slices );
  
  arma_debug_check
    (
      (
      ( (n_rows > 0x0FFF) || (n_cols > 0x0FFF) || (n_slices > 0xFF) )
        ? ( (float(n_rows) * float(n_cols) * float(n_slices)) > float(ARMA_MAX_UWORD) )
        : false
      ),
    "Cube::init(): requested size is too large"
    );
  
  if(n_elem <= Cube_prealloc::mem_n_elem)
    {
    access::rw(mem) = mem_local;
    }
  else
    {
    arma_extra_debug_print("Cube::init(): allocating memory");
    
    access::rw(mem) = memory::acquire<eT>(n_elem);
    
    arma_check_bad_alloc( (mem == 0), "Cube::init(): out of memory" );
    }
  
  
  if(n_elem == 0)
    {
    access::rw(n_rows)       = 0;
    access::rw(n_cols)       = 0;
    access::rw(n_elem_slice) = 0;
    access::rw(n_slices)     = 0;
    }
  else
    {
    create_mat();
    }
  }



//! internal cube construction; if the requested size is small enough, memory from the stack is used.
//! otherwise memory is allocated via 'new'
template<typename eT>
inline
void
Cube<eT>::init_warm(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_n_rows = %d, in_n_cols = %d, in_n_slices = %d") % in_n_rows % in_n_cols % in_n_slices );
  
  if( (n_rows == in_n_rows) && (n_cols == in_n_cols) && (n_slices == in_n_slices) )
    {
    return;
    }
  
  const uword t_mem_state = mem_state;
  
  bool  err_state = false;
  char* err_msg   = 0;
  
  arma_debug_set_error
    (
    err_state,
    err_msg,
    (t_mem_state == 3),
    "Cube::init(): size is fixed and hence cannot be changed"
    );
  
  arma_debug_set_error
    (
    err_state,
    err_msg,
      (
      ( (in_n_rows > 0x0FFF) || (in_n_cols > 0x0FFF) || (in_n_slices > 0xFF) )
        ? ( (float(in_n_rows) * float(in_n_cols) * float(in_n_slices)) > float(ARMA_MAX_UWORD) )
        : false
      ),
    "Cube::init(): requested size is too large"
    );
  
  arma_debug_check(err_state, err_msg);
  
  const uword old_n_elem = n_elem;
  const uword new_n_elem = in_n_rows * in_n_cols * in_n_slices;
  
  if(old_n_elem == new_n_elem)
    {
    delete_mat();
    
    if(new_n_elem > 0)
      {
      access::rw(n_rows)       = in_n_rows;
      access::rw(n_cols)       = in_n_cols;
      access::rw(n_elem_slice) = in_n_rows*in_n_cols;
      access::rw(n_slices)     = in_n_slices;
      
      create_mat();
      }
    }
  else
    {
    arma_debug_check( (t_mem_state == 2), "Cube::init(): requested size is not compatible with the size of auxiliary memory" );
    
    delete_mat();
    
    if(t_mem_state == 0)
      {
      if(n_elem > Cube_prealloc::mem_n_elem )
        {
        arma_extra_debug_print("Cube::init(): freeing memory");
        
        memory::release( access::rw(mem) );
        }
      }
    
    access::rw(mem_state) = 0;
    
    if(new_n_elem <= Cube_prealloc::mem_n_elem)
      {
      access::rw(mem) = mem_local;
      }
    else
      {
      arma_extra_debug_print("Cube::init(): allocating memory");
      
      access::rw(mem) = memory::acquire<eT>(new_n_elem);
      
      arma_check_bad_alloc( (mem == 0), "Cube::init(): out of memory" );
      }
    
    if(new_n_elem > 0)
      {
      access::rw(n_rows)       = in_n_rows;
      access::rw(n_cols)       = in_n_cols;
      access::rw(n_elem_slice) = in_n_rows*in_n_cols;
      access::rw(n_slices)     = in_n_slices;
      access::rw(n_elem)       = new_n_elem;
      
      create_mat();
      }
    }
  
  
  if(new_n_elem == 0)
    {
    access::rw(n_rows)       = 0;
    access::rw(n_cols)       = 0;
    access::rw(n_elem_slice) = 0;
    access::rw(n_slices)     = 0;
    access::rw(n_elem)       = 0;
    }
  }



//! for constructing a complex cube out of two non-complex cubes
template<typename eT>
template<typename T1, typename T2>
inline
void
Cube<eT>::init
  (
  const BaseCube<typename Cube<eT>::pod_type,T1>& X,
  const BaseCube<typename Cube<eT>::pod_type,T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          T;
  
  arma_type_check(( is_complex<eT>::value == false ));   //!< compile-time abort if eT isn't std::complex
  arma_type_check(( is_complex< T>::value == true  ));   //!< compile-time abort if T is std::complex
  
  arma_type_check(( is_same_type< std::complex<T>, eT >::value == false ));   //!< compile-time abort if types are not compatible
  
  const ProxyCube<T1> PX(X.get_ref());
  const ProxyCube<T2> PY(Y.get_ref());
  
  arma_assert_same_size(PX, PY, "Cube()");
  
  const uword local_n_rows   = PX.get_n_rows();
  const uword local_n_cols   = PX.get_n_cols();
  const uword local_n_slices = PX.get_n_slices();
  
  init_warm(local_n_rows, local_n_cols, local_n_slices);
  
  eT* out_mem = (*this).memptr();
  
  const bool prefer_at_accessor = ( ProxyCube<T1>::prefer_at_accessor || ProxyCube<T2>::prefer_at_accessor );
  
  if(prefer_at_accessor == false)
    {
    typedef typename ProxyCube<T1>::ea_type ea_type1;
    typedef typename ProxyCube<T2>::ea_type ea_type2;
    
    const uword N = n_elem;
    
    ea_type1 A = PX.get_ea();
    ea_type2 B = PY.get_ea();
        
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = std::complex<T>(A[i], B[i]);
      }
    }
  else
    {
    for(uword uslice = 0; uslice < local_n_slices; ++uslice)
    for(uword ucol   = 0;   ucol < local_n_cols;   ++ucol  )
    for(uword urow   = 0;   urow < local_n_rows;   ++urow  )
      {
      *out_mem = std::complex<T>( PX.at(urow,ucol,uslice), PY.at(urow,ucol,uslice) );
      out_mem++;
      }
    }
  }



//! try to steal the memory from a given cube; 
//! if memory can't be stolen, copy the given cube
template<typename eT>
inline
void
Cube<eT>::steal_mem(Cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    if( (x.mem_state == 0) && (x.n_elem > Cube_prealloc::mem_n_elem) )
      {
      reset();
      
      const uword x_n_slices = x.n_slices;
      
      access::rw(n_rows)       = x.n_rows;
      access::rw(n_cols)       = x.n_cols;
      access::rw(n_elem_slice) = x.n_elem_slice;
      access::rw(n_slices)     = x_n_slices;
      access::rw(n_elem)       = x.n_elem;
      access::rw(mem)          = x.mem;
      
      if(x_n_slices > Cube_prealloc::mat_ptrs_size)
        {
        access::rw(  mat_ptrs) = x.mat_ptrs;
        access::rw(x.mat_ptrs) = 0;
        }
      else
        {
        access::rw(mat_ptrs) = const_cast< const Mat<eT>** >(mat_ptrs_local);
        
        for(uword i=0; i < x_n_slices; ++i)
          {
            mat_ptrs[i] = x.mat_ptrs[i];
          x.mat_ptrs[i] = 0;
          }
        }
      
      access::rw(x.n_rows)       = 0;
      access::rw(x.n_cols)       = 0;
      access::rw(x.n_elem_slice) = 0;
      access::rw(x.n_slices)     = 0;
      access::rw(x.n_elem)       = 0;
      access::rw(x.mem)          = 0;
      }
    else
      {
      (*this).operator=(x);
      }
    }
  }



template<typename eT>
inline
void
Cube<eT>::delete_mat()
  {
  arma_extra_debug_sigprint();
  
  for(uword uslice = 0; uslice < n_slices; ++uslice)
    {
    delete access::rw(mat_ptrs[uslice]);
    }
  
  if(mem_state <= 2)
    {
    if(n_slices > Cube_prealloc::mat_ptrs_size)
      {
      delete [] mat_ptrs;
      }
    }
  }



template<typename eT>
inline
void
Cube<eT>::create_mat()
  {
  arma_extra_debug_sigprint();
  
  if(mem_state <= 2)
    {
    if(n_slices <= Cube_prealloc::mat_ptrs_size)
      {
      access::rw(mat_ptrs) = const_cast< const Mat<eT>** >(mat_ptrs_local);
      }
    else
      {
      access::rw(mat_ptrs) = new(std::nothrow) const Mat<eT>*[n_slices];
      
      arma_check_bad_alloc( (mat_ptrs == 0), "Cube::create_mat(): out of memory" );
      }
    }
  
  for(uword uslice = 0; uslice < n_slices; ++uslice)
    {
    mat_ptrs[uslice] = new Mat<eT>('j', slice_memptr(uslice), n_rows, n_cols);
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
  
  init_warm(1,1,1);
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
  
  arrayops::inplace_plus( memptr(), val, n_elem );
  
  return *this;
  }



//! In-place subtraction of a scalar from all elements of the cube
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_minus( memptr(), val, n_elem );
  
  return *this;
  }



//! In-place multiplication of all elements of the cube with a scalar
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_mul( memptr(), val, n_elem );
  
  return *this;
  }



//! In-place division of all elements of the cube with a scalar
template<typename eT>
arma_inline
const Cube<eT>&
Cube<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_div( memptr(), val, n_elem );
  
  return *this;
  }



//! construct a cube from a given cube
template<typename eT>
inline
Cube<eT>::Cube(const Cube<eT>& x)
  : n_rows(x.n_rows)
  , n_cols(x.n_cols)
  , n_elem_slice(x.n_elem_slice)
  , n_slices(x.n_slices)
  , n_elem(x.n_elem)
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  arma_extra_debug_sigprint(arma_boost::format("this = %x   in_cube = %x") % this % &x);
  
  init_cold();
  
  arrayops::copy( memptr(), x.mem, n_elem );
  }



//! construct a cube from a given cube
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator=(const Cube<eT>& x)
  {
  arma_extra_debug_sigprint(arma_boost::format("this = %x   in_cube = %x") % this % &x);
  
  if(this != &x)
    {
    init_warm(x.n_rows, x.n_cols, x.n_slices);
    
    arrayops::copy( memptr(), x.mem, n_elem );
    }
  
  return *this;
  }



//! construct a cube from a given auxiliary array of eTs.
//! if copy_aux_mem is true, new memory is allocated and the array is copied.
//! if copy_aux_mem is false, the auxiliary array is used directly (without allocating memory and copying).
//! note that in the latter case 
//! the default is to copy the array.

template<typename eT>
inline
Cube<eT>::Cube(eT* aux_mem, const uword aux_n_rows, const uword aux_n_cols, const uword aux_n_slices, const bool copy_aux_mem, const bool strict)
  : n_rows      ( aux_n_rows                          )
  , n_cols      ( aux_n_cols                          )
  , n_elem_slice( aux_n_rows*aux_n_cols               )
  , n_slices    ( aux_n_slices                        )
  , n_elem      ( aux_n_rows*aux_n_cols*aux_n_slices  )
  , mem_state   ( copy_aux_mem ? 0 : (strict ? 2 : 1) )
  , mat_ptrs    ( 0                                   )
  , mem         ( copy_aux_mem ? 0 : aux_mem          )
  {
  arma_extra_debug_sigprint_this(this);
  
  if(copy_aux_mem == true)
    {
    init_cold();
    
    arrayops::copy( memptr(), aux_mem, n_elem );
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
Cube<eT>::Cube(const eT* aux_mem, const uword aux_n_rows, const uword aux_n_cols, const uword aux_n_slices)
  : n_rows(aux_n_rows)
  , n_cols(aux_n_cols)
  , n_elem_slice(aux_n_rows*aux_n_cols)
  , n_slices(aux_n_slices)
  , n_elem(aux_n_rows*aux_n_cols*aux_n_slices)
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold();
  
  arrayops::copy( memptr(), aux_mem, n_elem );
  }



//! in-place cube addition
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator+=(const Cube<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(*this, m, "cube addition");
  
  arrayops::inplace_plus( memptr(), m.memptr(), n_elem );
  
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
  
  arrayops::inplace_minus( memptr(), m.memptr(), n_elem );
  
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
  
  arrayops::inplace_mul( memptr(), m.memptr(), n_elem );
  
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
  
  arrayops::inplace_div( memptr(), m.memptr(), n_elem );
  
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
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  init(A,B);
  }



//! construct a cube from a subview_cube instance (e.g. construct a cube from a delayed subcube operation)
template<typename eT>
inline
Cube<eT>::Cube(const subview_cube<eT>& X)
  : n_rows(X.n_rows)
  , n_cols(X.n_cols)
  , n_elem_slice(X.n_elem_slice)
  , n_slices(X.n_slices)
  , n_elem(X.n_elem)
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold();
  
  subview_cube<eT>::extract(*this, X);
  }



//! construct a cube from a subview_cube instance (e.g. construct a cube from a delayed subcube operation)
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::operator=(const subview_cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool alias = (this == &(X.m));
  
  if(alias == false)
    {
    init_warm(X.n_rows, X.n_cols, X.n_slices);
    
    subview_cube<eT>::extract(*this, X);
    }
  else
    {
    Cube<eT> tmp(X);
    
    steal_mem(tmp);
    }
  
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
Cube<eT>::slice(const uword in_slice)
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
Cube<eT>::slice(const uword in_slice) const
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
Cube<eT>::slices(const uword in_slice1, const uword in_slice2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_slice1 > in_slice2) || (in_slice2 >= n_slices),
    "Cube::slices(): indices out of bounds or incorrectly used"
    );
  
  const uword subcube_n_slices = in_slice2 - in_slice1 + 1;
  
  return subview_cube<eT>(*this, 0, 0, in_slice1, n_rows, n_cols, subcube_n_slices);
  }



//! creation of subview_cube (subcube comprised of specified slices)
template<typename eT>
arma_inline
const subview_cube<eT>
Cube<eT>::slices(const uword in_slice1, const uword in_slice2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_slice1 > in_slice2) || (in_slice2 >= n_slices),
    "Cube::rows(): indices out of bounds or incorrectly used"
    );
  
  const uword subcube_n_slices = in_slice2 - in_slice1 + 1;
  
  return subview_cube<eT>(*this, 0, 0, in_slice1, n_rows, n_cols, subcube_n_slices);
  }



//! creation of subview_cube (generic subcube)
template<typename eT>
arma_inline
subview_cube<eT>
Cube<eT>::subcube(const uword in_row1, const uword in_col1, const uword in_slice1, const uword in_row2, const uword in_col2, const uword in_slice2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 >  in_row2) || (in_col1 >  in_col2) || (in_slice1 >  in_slice2) ||
    (in_row2 >= n_rows)  || (in_col2 >= n_cols)  || (in_slice2 >= n_slices),
    "Cube::subcube(): indices out of bounds or incorrectly used"
    );
  
  const uword subcube_n_rows   = in_row2   - in_row1   + 1;
  const uword subcube_n_cols   = in_col2   - in_col1   + 1;
  const uword subcube_n_slices = in_slice2 - in_slice1 + 1;
  
  return subview_cube<eT>(*this, in_row1, in_col1, in_slice1, subcube_n_rows, subcube_n_cols, subcube_n_slices);
  }



//! creation of subview_cube (generic subcube)
template<typename eT>
arma_inline
const subview_cube<eT>
Cube<eT>::subcube(const uword in_row1, const uword in_col1, const uword in_slice1, const uword in_row2, const uword in_col2, const uword in_slice2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 >  in_row2) || (in_col1 >  in_col2) || (in_slice1 >  in_slice2) ||
    (in_row2 >= n_rows)  || (in_col2 >= n_cols)  || (in_slice2 >= n_slices),
    "Cube::subcube(): indices out of bounds or incorrectly used"
    );
    
  const uword subcube_n_rows   = in_row2   - in_row1   + 1;
  const uword subcube_n_cols   = in_col2   - in_col1   + 1;
  const uword subcube_n_slices = in_slice2 - in_slice1 + 1;
  
  return subview_cube<eT>(*this, in_row1, in_col1, in_slice1, subcube_n_rows, subcube_n_cols, subcube_n_slices);
  }



//! creation of subview_cube (generic subcube)
template<typename eT>
inline
subview_cube<eT>
Cube<eT>::subcube(const span& row_span, const span& col_span, const span& slice_span)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all   = row_span.whole;
  const bool col_all   = col_span.whole;
  const bool slice_all = slice_span.whole;
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  const uword in_row1          = row_all   ? 0              : row_span.a;
  const uword in_row2          =                              row_span.b;
  const uword subcube_n_rows   = row_all   ? local_n_rows   : in_row2 - in_row1 + 1;
  
  const uword in_col1          = col_all   ? 0              : col_span.a;
  const uword in_col2          =                              col_span.b;
  const uword subcube_n_cols   = col_all   ? local_n_cols   : in_col2 - in_col1 + 1;
  
  const uword in_slice1        = slice_all ? 0              : slice_span.a;
  const uword in_slice2        =                              slice_span.b;
  const uword subcube_n_slices = slice_all ? local_n_slices : in_slice2 - in_slice1 + 1;
  
  arma_debug_check
    (
    ( row_all   ? false : ((in_row1   >  in_row2)   || (in_row2   >= local_n_rows))   )
    ||
    ( col_all   ? false : ((in_col1   >  in_col2)   || (in_col2   >= local_n_cols))   )
    ||
    ( slice_all ? false : ((in_slice1 >  in_slice2) || (in_slice2 >= local_n_slices)) )
    ,
    "Cube::subcube(): indices out of bounds or incorrectly used"
    );
  
  return subview_cube<eT>(*this, in_row1, in_col1, in_slice1, subcube_n_rows, subcube_n_cols, subcube_n_slices);
  }



//! creation of subview_cube (generic subcube)
template<typename eT>
inline
const subview_cube<eT>
Cube<eT>::subcube(const span& row_span, const span& col_span, const span& slice_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool row_all   = row_span.whole;
  const bool col_all   = col_span.whole;
  const bool slice_all = slice_span.whole;
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  const uword in_row1          = row_all   ? 0              : row_span.a;
  const uword in_row2          =                              row_span.b;
  const uword subcube_n_rows   = row_all   ? local_n_rows   : in_row2 - in_row1 + 1;
  
  const uword in_col1          = col_all   ? 0              : col_span.a;
  const uword in_col2          =                              col_span.b;
  const uword subcube_n_cols   = col_all   ? local_n_cols   : in_col2 - in_col1 + 1;
  
  const uword in_slice1        = slice_all ? 0              : slice_span.a;
  const uword in_slice2        =                              slice_span.b;
  const uword subcube_n_slices = slice_all ? local_n_slices : in_slice2 - in_slice1 + 1;
  
  arma_debug_check
    (
    ( row_all   ? false : ((in_row1   >  in_row2)   || (in_row2   >= local_n_rows))   )
    ||
    ( col_all   ? false : ((in_col1   >  in_col2)   || (in_col2   >= local_n_cols))   )
    ||
    ( slice_all ? false : ((in_slice1 >  in_slice2) || (in_slice2 >= local_n_slices)) )
    ,
    "Cube::subcube(): indices out of bounds or incorrectly used"
    );
  
  return subview_cube<eT>(*this, in_row1, in_col1, in_slice1, subcube_n_rows, subcube_n_cols, subcube_n_slices);
  }



template<typename eT>
inline
subview_cube<eT>
Cube<eT>::operator()(const span& row_span, const span& col_span, const span& slice_span)
  {
  arma_extra_debug_sigprint();
  
  return (*this).subcube(row_span, col_span, slice_span);
  }



template<typename eT>
inline
const subview_cube<eT>
Cube<eT>::operator()(const span& row_span, const span& col_span, const span& slice_span) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).subcube(row_span, col_span, slice_span);
  }



//! remove specified slice
template<typename eT>
inline
void
Cube<eT>::shed_slice(const uword slice_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( slice_num >= n_slices, "Cube::shed_slice(): out of bounds");
  
  shed_slices(slice_num, slice_num);
  }



//! remove specified slices
template<typename eT>
inline
void
Cube<eT>::shed_slices(const uword in_slice1, const uword in_slice2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_slice1 > in_slice2) || (in_slice2 >= n_slices),
    "Cube::shed_slices(): indices out of bounds or incorrectly used"
    );
  
  const uword n_keep_front = in_slice1;
  const uword n_keep_back  = n_slices - (in_slice2 + 1);
  
  Cube<eT> X(n_rows, n_cols, n_keep_front + n_keep_back);
  
  if(n_keep_front > 0)
    {
    X.slices( 0, (n_keep_front-1) ) = slices( 0, (in_slice1-1) );
    }
  
  if(n_keep_back > 0)
    {
    X.slices( n_keep_front,  (n_keep_front+n_keep_back-1) ) = slices( (in_slice2+1), (n_slices-1) );
    }
  
  steal_mem(X);
  }



//! insert N slices at the specified slice position,
//! optionally setting the elements of the inserted slices to zero
template<typename eT>
inline
void
Cube<eT>::insert_slices(const uword slice_num, const uword N, const bool set_to_zero)
  {
  arma_extra_debug_sigprint();
  
  const uword t_n_slices = n_slices;
  
  const uword A_n_slices = slice_num;
  const uword B_n_slices = t_n_slices - slice_num;
  
  // insertion at slice_num == n_slices is in effect an append operation
  arma_debug_check( (slice_num > t_n_slices), "Cube::insert_slices(): out of bounds");
  
  if(N > 0)
    {
    Cube<eT> out(n_rows, n_cols, t_n_slices + N);
    
    if(A_n_slices > 0)
      {
      out.slices(0, A_n_slices-1) = slices(0, A_n_slices-1);
      }
    
    if(B_n_slices > 0)
      {
      out.slices(slice_num + N, t_n_slices + N - 1) = slices(slice_num, t_n_slices-1);
      }
    
    if(set_to_zero == true)
      {
      //out.slices(slice_num, slice_num + N - 1).zeros();
      
      for(uword i=slice_num; i < (slice_num + N); ++i)
        {
        out.slice(i).zeros();
        }
      }
    
    steal_mem(out);
    }
  }



//! insert the given object at the specified slice position; 
//! the given object must have the same number of rows and columns as the cube
template<typename eT>
template<typename T1>
inline
void
Cube<eT>::insert_slices(const uword slice_num, const BaseCube<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(X.get_ref());
  const Cube<eT>& C   = tmp.M;
  
  const uword N = C.n_slices;
  
  const uword t_n_slices = n_slices;
  
  const uword A_n_slices = slice_num;
  const uword B_n_slices = t_n_slices - slice_num;
  
  // insertion at slice_num == n_slices is in effect an append operation
  arma_debug_check( (slice_num  >  t_n_slices), "Cube::insert_slices(): out of bounds");
  
  arma_debug_check
    (
    ( (C.n_rows != n_rows) || (C.n_cols != n_cols) ),
    "Cube::insert_slices(): given object has an incompatible dimensions"
    );
  
  if(N > 0)
    {
    Cube<eT> out(n_rows, n_cols, t_n_slices + N);
    
    if(A_n_slices > 0)
      {
      out.slices(0, A_n_slices-1) = slices(0, A_n_slices-1);
      }
    
    if(B_n_slices > 0)
      {
      out.slices(slice_num + N, t_n_slices + N - 1) = slices(slice_num, t_n_slices - 1);
      }
    
    out.slices(slice_num, slice_num + N - 1) = C;
    
    steal_mem(out);
    }
  }



//! create a cube from OpCube, i.e. run the previously delayed unary operations
template<typename eT>
template<typename gen_type>
inline
Cube<eT>::Cube(const GenCube<eT, gen_type>& X)
  : n_rows(X.n_rows)
  , n_cols(X.n_cols)
  , n_elem_slice(X.n_rows*X.n_cols)
  , n_slices(X.n_slices)
  , n_elem(X.n_rows*X.n_cols*X.n_slices)
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold();
  
  X.apply(*this);
  }



template<typename eT>
template<typename gen_type>
inline
const Cube<eT>&
Cube<eT>::operator=(const GenCube<eT, gen_type>& X)
  {
  arma_extra_debug_sigprint();
  
  init_warm(X.n_rows, X.n_cols, X.n_slices);
  
  X.apply(*this);
  
  return *this;
  }



template<typename eT>
template<typename gen_type>
inline
const Cube<eT>&
Cube<eT>::operator+=(const GenCube<eT, gen_type>& X)
  {
  arma_extra_debug_sigprint();
  
  X.apply_inplace_plus(*this);
  
  return *this;
  }



template<typename eT>
template<typename gen_type>
inline
const Cube<eT>&
Cube<eT>::operator-=(const GenCube<eT, gen_type>& X)
  {
  arma_extra_debug_sigprint();
  
  X.apply_inplace_minus(*this);
  
  return *this;
  }



template<typename eT>
template<typename gen_type>
inline
const Cube<eT>&
Cube<eT>::operator%=(const GenCube<eT, gen_type>& X)
  {
  arma_extra_debug_sigprint();
  
  X.apply_inplace_schur(*this);
  
  return *this;
  }



template<typename eT>
template<typename gen_type>
inline
const Cube<eT>&
Cube<eT>::operator/=(const GenCube<eT, gen_type>& X)
  {
  arma_extra_debug_sigprint();
  
  X.apply_inplace_div(*this);
  
  return *this;
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
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);

  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
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

  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
  const Cube<eT> m(X);
  
  return (*this).operator/=(m);
  }



//! create a cube from eOpCube, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename eop_type>
inline
Cube<eT>::Cube(const eOpCube<T1, eop_type>& X)
  : n_rows(X.get_n_rows())
  , n_cols(X.get_n_cols())
  , n_elem_slice(X.get_n_elem_slice())
  , n_slices(X.get_n_slices())
  , n_elem(X.get_n_elem())
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
  init_cold();
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
  const bool bad_alias = ( X.P.has_subview  &&  X.P.is_alias(*this) );
  
  if(bad_alias == false)
    {
    init_warm(X.get_n_rows(), X.get_n_cols(), X.get_n_slices());
    
    eop_type::apply(*this, X);
    }
  else
    {
    Cube<eT> tmp(X);
    
    steal_mem(tmp);
    }
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));

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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));

  eop_type::apply_inplace_div(*this, X);
  
  return *this;
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
Cube<eT>::Cube(const mtOpCube<eT, T1, op_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  op_type::apply(*this, X);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator=(const mtOpCube<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  op_type::apply(*this, X);
  
  return *this;
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator+=(const mtOpCube<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> m(X);
  
  return (*this).operator+=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator-=(const mtOpCube<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> m(X);
  
  return (*this).operator-=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator%=(const mtOpCube<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> m(X);
  
  return (*this).operator%=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename op_type>
inline
const Cube<eT>&
Cube<eT>::operator/=(const mtOpCube<eT, T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> m(X);
  
  return (*this).operator/=(m);
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
  , mem_state(0)
  , mat_ptrs()
  , mem()
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
  const Cube<eT> m(X);
  
  return (*this).operator/=(m);
  }



//! create a cube from eGlue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
Cube<eT>::Cube(const eGlueCube<T1, T2, eglue_type>& X)
  : n_rows(X.get_n_rows())
  , n_cols(X.get_n_cols())
  , n_elem_slice(X.get_n_elem_slice())
  , n_slices(X.get_n_slices())
  , n_elem(X.get_n_elem())
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
  init_cold();
  
  eglue_type::apply(*this, X);
  }



//! create a cube from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename eglue_type>
inline
const Cube<eT>&
Cube<eT>::operator=(const eGlueCube<T1, T2, eglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
  const bool bad_alias = ( (X.P1.has_subview  &&  X.P1.is_alias(*this))  ||  (X.P2.has_subview  &&  X.P2.is_alias(*this)) );
  
  if(bad_alias == false)
    {
    init_warm(X.get_n_rows(), X.get_n_cols(), X.get_n_slices());
    
    eglue_type::apply(*this, X);
    }
  else
    {
    Cube<eT> tmp(X);
    
    steal_mem(tmp);
    }
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
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
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::value == false ));
  arma_type_check(( is_same_type< eT, typename T2::elem_type >::value == false ));
  
  eglue_type::apply_inplace_div(*this, X);
  
  return *this;
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Cube<eT>::Cube(const mtGlueCube<eT, T1, T2, glue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem_slice(0)
  , n_slices(0)
  , n_elem(0)
  , mem_state(0)
  , mat_ptrs()
  , mem()
  {
  arma_extra_debug_sigprint_this(this);
  
  glue_type::apply(*this, X);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator=(const mtGlueCube<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  glue_type::apply(*this, X);
  
  return *this;
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator+=(const mtGlueCube<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> m(X);
  
  return (*this).operator+=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator-=(const mtGlueCube<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> m(X);
  
  return (*this).operator-=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator%=(const mtGlueCube<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> m(X);
  
  return (*this).operator%=(m);
  }



//! EXPERIMENTAL
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Cube<eT>&
Cube<eT>::operator/=(const mtGlueCube<eT, T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> m(X);
  
  return (*this).operator/=(m);
  }



//! linear element accessor (treats the cube as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
arma_warn_unused
eT&
Cube<eT>::operator() (const uword i)
  {
  arma_debug_check( (i >= n_elem), "Cube::operator(): index out of bounds");
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the cube as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
arma_warn_unused
eT
Cube<eT>::operator() (const uword i) const
  {
  arma_debug_check( (i >= n_elem), "Cube::operator(): index out of bounds");
  return mem[i];
  }


//! linear element accessor (treats the cube as a vector); no bounds check.  
template<typename eT>
arma_inline
arma_warn_unused
eT&
Cube<eT>::operator[] (const uword i)
  {
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the cube as a vector); no bounds check
template<typename eT>
arma_inline
arma_warn_unused
eT
Cube<eT>::operator[] (const uword i) const
  {
  return mem[i];
  }



//! linear element accessor (treats the cube as a vector); no bounds check.  
template<typename eT>
arma_inline
arma_warn_unused
eT&
Cube<eT>::at(const uword i)
  {
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the cube as a vector); no bounds check
template<typename eT>
arma_inline
arma_warn_unused
eT
Cube<eT>::at(const uword i) const
  {
  return mem[i];
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
arma_warn_unused
eT&
Cube<eT>::operator() (const uword in_row, const uword in_col, const uword in_slice)
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
arma_warn_unused
eT
Cube<eT>::operator() (const uword in_row, const uword in_col, const uword in_slice) const
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
arma_warn_unused
eT&
Cube<eT>::at(const uword in_row, const uword in_col, const uword in_slice)
  {
  return access::rw( mem[in_slice*n_elem_slice + in_col*n_rows + in_row] );
  }



//! element accessor; no bounds check
template<typename eT>
arma_inline
arma_warn_unused
eT
Cube<eT>::at(const uword in_row, const uword in_col, const uword in_slice) const
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
arma_warn_unused
bool
Cube<eT>::is_finite() const
  {
  return arrayops::is_finite( memptr(), n_elem );
  }



//! returns true if the cube has no elements
template<typename eT>
arma_inline
arma_warn_unused
bool
Cube<eT>::is_empty() const
  {
  return (n_elem == 0);
  }



//! returns true if the given index is currently in range
template<typename eT>
arma_inline
arma_warn_unused
bool
Cube<eT>::in_range(const uword i) const
  {
  return (i < n_elem);
  }



//! returns true if the given start and end indices are currently in range
template<typename eT>
arma_inline
arma_warn_unused
bool
Cube<eT>::in_range(const span& x) const
  {
  arma_extra_debug_sigprint();
  
  if(x.whole == true)
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
template<typename eT>
arma_inline
arma_warn_unused
bool
Cube<eT>::in_range(const uword in_row, const uword in_col, const uword in_slice) const
  {
  return ( (in_row < n_rows) && (in_col < n_cols) && (in_slice < n_slices) );
  }



template<typename eT>
inline
arma_warn_unused
bool
Cube<eT>::in_range(const span& row_span, const span& col_span, const span& slice_span) const
  {
  arma_extra_debug_sigprint();
  
  const uword in_row1   = row_span.a;
  const uword in_row2   = row_span.b;
  
  const uword in_col1   = col_span.a;
  const uword in_col2   = col_span.b;
  
  const uword in_slice1 = slice_span.a;
  const uword in_slice2 = slice_span.b;
  
  
  const bool rows_ok   = row_span.whole   ? true : ( (in_row1   <= in_row2)   && (in_row2   < n_rows)   );
  const bool cols_ok   = col_span.whole   ? true : ( (in_col1   <= in_col2)   && (in_col2   < n_cols)   );
  const bool slices_ok = slice_span.whole ? true : ( (in_slice1 <= in_slice2) && (in_slice2 < n_slices) );
  
  
  return ( (rows_ok == true) && (cols_ok == true) && (slices_ok == true) );
  }



//! returns a pointer to array of eTs used by the cube
template<typename eT>
arma_inline
arma_warn_unused
eT*
Cube<eT>::memptr()
  {
  return const_cast<eT*>(mem);
  }



//! returns a pointer to array of eTs used by the cube
template<typename eT>
arma_inline
arma_warn_unused
const eT*
Cube<eT>::memptr() const
  {
  return mem;
  }



//! returns a pointer to array of eTs used by the specified slice in the cube
template<typename eT>
arma_inline
arma_warn_unused
eT*
Cube<eT>::slice_memptr(const uword uslice)
  {
  return const_cast<eT*>( &mem[ uslice*n_elem_slice ] );
  }



//! returns a pointer to array of eTs used by the specified slice in the cube
template<typename eT>
arma_inline
arma_warn_unused
const eT*
Cube<eT>::slice_memptr(const uword uslice) const
  {
  return &mem[ uslice*n_elem_slice ];
  }



//! returns a pointer to array of eTs used by the specified slice in the cube
template<typename eT>
arma_inline
arma_warn_unused
eT*
Cube<eT>::slice_colptr(const uword uslice, const uword col)
  {
  return const_cast<eT*>( &mem[ uslice*n_elem_slice + col*n_rows] );
  }



//! returns a pointer to array of eTs used by the specified slice in the cube
template<typename eT>
arma_inline
arma_warn_unused
const eT*
Cube<eT>::slice_colptr(const uword uslice, const uword col) const
  {
  return &mem[ uslice*n_elem_slice + col*n_rows ];
  }



//! print contents of the cube (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's state are restored to their original values.
template<typename eT>
inline
void
Cube<eT>::impl_print(const std::string& extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    ARMA_DEFAULT_OSTREAM << extra_text << '\n';
    }
  
  arma_ostream::print(ARMA_DEFAULT_OSTREAM, *this, true);
  }


//! print contents of the cube to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's state are restored to their original values.
template<typename eT>
inline
void
Cube<eT>::impl_print(std::ostream& user_stream, const std::string& extra_text) const
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
//! the stream's state are used as is and are not modified
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Cube<eT>::impl_raw_print(const std::string& extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    ARMA_DEFAULT_OSTREAM << extra_text << '\n';
    }
  
  arma_ostream::print(ARMA_DEFAULT_OSTREAM, *this, false);
  }



//! print contents of the cube to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the stream's state are used as is and are not modified.
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Cube<eT>::impl_raw_print(std::ostream& user_stream, const std::string& extra_text) const
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
Cube<eT>::set_size(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices)
  {
  arma_extra_debug_sigprint();
  
  init_warm(in_n_rows, in_n_cols, in_n_slices);
  }



//! change the cube to have user specified dimensions (data is preserved)
template<typename eT>
inline
void
Cube<eT>::reshape(const uword in_rows, const uword in_cols, const uword in_slices, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  *this = arma::reshape(*this, in_rows, in_cols, in_slices, dim);
  }



//! change the cube to have user specified dimensions (data is preserved)
template<typename eT>
inline
void
Cube<eT>::resize(const uword in_rows, const uword in_cols, const uword in_slices)
  {
  arma_extra_debug_sigprint();
  
  *this = arma::resize(*this, in_rows, in_cols, in_slices);
  }



//! change the cube (without preserving data) to have the same dimensions as the given cube 
template<typename eT>
template<typename eT2>
inline
void
Cube<eT>::copy_size(const Cube<eT2>& m)
  {
  arma_extra_debug_sigprint();
  
  init_warm(m.n_rows, m.n_cols, m.n_slices);
  }



//! fill the cube with the specified value
template<typename eT>
inline
const Cube<eT>&
Cube<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set( memptr(), val, n_elem );
  
  return *this;
  }



template<typename eT>
inline
const Cube<eT>&
Cube<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  return (*this).fill(eT(0));
  }



template<typename eT>
inline
const Cube<eT>&
Cube<eT>::zeros(const uword in_rows, const uword in_cols, const uword in_slices)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_rows, in_cols, in_slices);
  
  return (*this).fill(eT(0));
  }



template<typename eT>
inline
const Cube<eT>&
Cube<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  return (*this).fill(eT(1));
  }



template<typename eT>
inline
const Cube<eT>&
Cube<eT>::ones(const uword in_rows, const uword in_cols, const uword in_slices)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_rows, in_cols, in_slices);
  
  return (*this).fill(eT(1));
  }



template<typename eT>
inline
const Cube<eT>&
Cube<eT>::randu()
  {
  arma_extra_debug_sigprint();
  
  eop_aux_randu<eT>::fill( memptr(), n_elem );
  
  return *this;
  }



template<typename eT>
inline
const Cube<eT>&
Cube<eT>::randu(const uword in_rows, const uword in_cols, const uword in_slices)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_rows, in_cols, in_slices);
  
  return (*this).randu();
  }



template<typename eT>
inline
const Cube<eT>&
Cube<eT>::randn()
  {
  arma_extra_debug_sigprint();
  
  eop_aux_randn<eT>::fill( memptr(), n_elem );
  
  return *this;
  }



template<typename eT>
inline
const Cube<eT>&
Cube<eT>::randn(const uword in_rows, const uword in_cols, const uword in_slices)
  {
  arma_extra_debug_sigprint();
  
  set_size(in_rows, in_cols, in_slices);
  
  return (*this).randn();
  }



template<typename eT>
inline
void
Cube<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  init_warm(0,0,0);
  }



template<typename eT>
template<typename T1>
inline
void
Cube<eT>::set_real(const BaseCube<typename Cube<eT>::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Cube_aux::set_real(*this, X);
  }



template<typename eT>
template<typename T1>
inline
void
Cube<eT>::set_imag(const BaseCube<typename Cube<eT>::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Cube_aux::set_imag(*this, X);
  }



template<typename eT>
inline
arma_warn_unused
eT
Cube<eT>::min() const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (n_elem == 0), "min(): object has no elements" );
  
  return op_min::direct_min(memptr(), n_elem);
  }



template<typename eT>
inline
arma_warn_unused
eT
Cube<eT>::max() const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (n_elem == 0), "max(): object has no elements" );
  
  return op_max::direct_max(memptr(), n_elem);
  }



template<typename eT>
inline
eT
Cube<eT>::min(uword& index_of_min_val) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (n_elem == 0), "min(): object has no elements" );
  
  return op_min::direct_min(memptr(), n_elem, index_of_min_val);
  }



template<typename eT>
inline
eT
Cube<eT>::max(uword& index_of_max_val) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (n_elem == 0), "max(): object has no elements" );
  
  return op_max::direct_max(memptr(), n_elem, index_of_max_val);
  }



template<typename eT>
inline
eT
Cube<eT>::min(uword& row_of_min_val, uword& col_of_min_val, uword& slice_of_min_val) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (n_elem == 0), "min(): object has no elements" );
  
  uword i;
  
  eT val = op_min::direct_min(memptr(), n_elem, i);
  
  const uword in_slice = i / n_elem_slice;
  const uword offset   = in_slice * n_elem_slice;
  const uword j        = i - offset;
  
    row_of_min_val = j % n_rows;
    col_of_min_val = j / n_rows;
  slice_of_min_val = in_slice;
  
  return val;
  }



template<typename eT>
inline
eT
Cube<eT>::max(uword& row_of_max_val, uword& col_of_max_val, uword& slice_of_max_val) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (n_elem == 0), "max(): object has no elements" );
  
  uword i;
  
  eT val = op_max::direct_max(memptr(), n_elem, i);
  
  const uword in_slice = i / n_elem_slice;
  const uword offset   = in_slice * n_elem_slice;
  const uword j        = i - offset;
  
    row_of_max_val = j % n_rows;
    col_of_max_val = j / n_rows;
  slice_of_max_val = in_slice;
  
  return val;
  }



//! save the cube to a file
template<typename eT>
inline
bool
Cube<eT>::save(const std::string name, const file_type type, const bool print_status) const
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
    
    case raw_binary:
      save_okay = diskio::save_raw_binary(*this, name);
      break;
    
    case arma_binary:
      save_okay = diskio::save_arma_binary(*this, name);
      break;
    
    case ppm_binary:
      save_okay = diskio::save_ppm_binary(*this, name);
      break;
    
    default:
      arma_warn(print_status, "Cube::save(): unsupported file type");
      save_okay = false;
    }
  
  arma_warn( (print_status && (save_okay == false)), "Cube::save(): couldn't write to ", name);
  
  return save_okay;
  }



//! save the cube to a stream
template<typename eT>
inline
bool
Cube<eT>::save(std::ostream& os, const file_type type, const bool print_status) const
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
    
    case raw_binary:
      save_okay = diskio::save_raw_binary(*this, os);
      break;
    
    case arma_binary:
      save_okay = diskio::save_arma_binary(*this, os);
      break;
    
    case ppm_binary:
      save_okay = diskio::save_ppm_binary(*this, os);
      break;
    
    default:
      arma_warn(print_status, "Cube::save(): unsupported file type");
      save_okay = false;
    }
  
  arma_warn( (print_status && (save_okay == false)), "Cube::save(): couldn't write to given stream");
  
  return save_okay;
  }



//! load a cube from a file
template<typename eT>
inline
bool
Cube<eT>::load(const std::string name, const file_type type, const bool print_status)
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
    
    case raw_binary:
      load_okay = diskio::load_raw_binary(*this, name, err_msg);
      break;
    
    case arma_binary:
      load_okay = diskio::load_arma_binary(*this, name, err_msg);
      break;
    
    case ppm_binary:
      load_okay = diskio::load_ppm_binary(*this, name, err_msg);
      break;
    
    default:
      arma_warn(print_status, "Cube::load(): unsupported file type");
      load_okay = false;
    }
  
  if( (print_status == true) && (load_okay == false) )
    {
    if(err_msg.length() > 0)
      {
      arma_warn(true, "Cube::load(): ", err_msg, name);
      }
    else
      {
      arma_warn(true, "Cube::load(): couldn't read ", name);
      }
    }
  
  if(load_okay == false)
    {
    (*this).reset();
    }
    
  return load_okay;
  }



//! load a cube from a stream
template<typename eT>
inline
bool
Cube<eT>::load(std::istream& is, const file_type type, const bool print_status)
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
    
    case raw_binary:
      load_okay = diskio::load_raw_binary(*this, is, err_msg);
      break;
    
    case arma_binary:
      load_okay = diskio::load_arma_binary(*this, is, err_msg);
      break;
    
    case ppm_binary:
      load_okay = diskio::load_ppm_binary(*this, is, err_msg);
      break;
    
    default:
      arma_warn(print_status, "Cube::load(): unsupported file type");
      load_okay = false;
    }
  
  
  if( (print_status == true) && (load_okay == false) )
    {
    if(err_msg.length() > 0)
      {
      arma_warn(true, "Cube::load(): ", err_msg, "the given stream");
      }
    else
      {
      arma_warn(true, "Cube::load(): couldn't load from the given stream");
      }
    }
  
  if(load_okay == false)
    {
    (*this).reset();
    }
    
  return load_okay;
  }



//! save the cube to a file, without printing any error messages
template<typename eT>
inline
bool
Cube<eT>::quiet_save(const std::string name, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).save(name, type, false);
  }



//! save the cube to a stream, without printing any error messages
template<typename eT>
inline
bool
Cube<eT>::quiet_save(std::ostream& os, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).save(os, type, false);
  }



//! load a cube from a file, without printing any error messages
template<typename eT>
inline
bool
Cube<eT>::quiet_load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  return (*this).load(name, type, false);
  }



//! load a cube from a stream, without printing any error messages
template<typename eT>
inline
bool
Cube<eT>::quiet_load(std::istream& is, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  return (*this).load(is, type, false);
  }



template<typename eT>
inline
typename Cube<eT>::iterator
Cube<eT>::begin()
  {
  arma_extra_debug_sigprint();
  
  return memptr();
  }



template<typename eT>
inline
typename Cube<eT>::const_iterator
Cube<eT>::begin() const
  {
  arma_extra_debug_sigprint();
  
  return memptr();
  }



template<typename eT>
inline
typename Cube<eT>::iterator
Cube<eT>::end()
  {
  arma_extra_debug_sigprint();
  
  return memptr() + n_elem;
  }



template<typename eT>
inline
typename Cube<eT>::const_iterator
Cube<eT>::end() const
  {
  arma_extra_debug_sigprint();
  
  return memptr() + n_elem;
  }
  


template<typename eT>
inline
typename Cube<eT>::slice_iterator
Cube<eT>::begin_slice(const uword slice_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (slice_num >= n_slices), "begin_slice(): index out of bounds");
  
  return slice_memptr(slice_num);
  }



template<typename eT>
inline
typename Cube<eT>::const_slice_iterator
Cube<eT>::begin_slice(const uword slice_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (slice_num >= n_slices), "begin_slice(): index out of bounds");
  
  return slice_memptr(slice_num);
  }



template<typename eT>
inline
typename Cube<eT>::slice_iterator
Cube<eT>::end_slice(const uword slice_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (slice_num >= n_slices), "end_slice(): index out of bounds");
  
  return slice_memptr(slice_num) + n_elem_slice;
  }



template<typename eT>
inline
typename Cube<eT>::const_slice_iterator
Cube<eT>::end_slice(const uword slice_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (slice_num >= n_slices), "end_slice(): index out of bounds");
  
  return slice_memptr(slice_num) + n_elem_slice;
  }



template<typename eT>
template<uword fixed_n_rows, uword fixed_n_cols, uword fixed_n_slices>
arma_inline
void
Cube<eT>::fixed<fixed_n_rows, fixed_n_cols, fixed_n_slices>::mem_setup()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(fixed_n_elem > 0)
    {
    access::rw(Cube<eT>::n_rows)       = fixed_n_rows;
    access::rw(Cube<eT>::n_cols)       = fixed_n_cols;
    access::rw(Cube<eT>::n_elem_slice) = fixed_n_rows * fixed_n_cols;
    access::rw(Cube<eT>::n_slices)     = fixed_n_slices;
    access::rw(Cube<eT>::n_elem)       = fixed_n_elem;
    access::rw(Cube<eT>::mem_state)    = 3;
    access::rw(Cube<eT>::mat_ptrs)     = const_cast< const Mat<eT>** >( \
                                         (fixed_n_slices > Cube_prealloc::mat_ptrs_size) ? mat_ptrs_local_extra : mat_ptrs_local );
    access::rw(Cube<eT>::mem)          = (fixed_n_elem   > Cube_prealloc::mem_n_elem)    ? mem_local_extra      : mem_local;
    
    create_mat();
    }
  else
    {
    access::rw(Cube<eT>::n_rows)       = 0;
    access::rw(Cube<eT>::n_cols)       = 0;
    access::rw(Cube<eT>::n_elem_slice) = 0;
    access::rw(Cube<eT>::n_slices)     = 0;
    access::rw(Cube<eT>::n_elem)       = 0;
    access::rw(Cube<eT>::mem_state)    = 3;
    access::rw(Cube<eT>::mat_ptrs)     = 0;
    access::rw(Cube<eT>::mem)          = 0;
    }
  }



//! prefix ++
template<typename eT>
arma_inline
void
Cube_aux::prefix_pp(Cube<eT>& x)
  {
        eT*   memptr = x.memptr();
  const uword n_elem = x.n_elem;
  
  uword i,j;

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
  const uword n_elem = x.n_elem;
  
  uword i,j;
  
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
  const uword n_elem = x.n_elem;

  uword i,j;

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
  const uword n_elem = x.n_elem;

  uword i,j;

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



template<typename eT, typename T1>
inline
void
Cube_aux::set_real(Cube<eT>& out, const BaseCube<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(X.get_ref());
  const Cube<eT>& A   = tmp.M;
  
  arma_debug_assert_same_size( out, A, "Cube::set_real()" );
  
  out = A;
  }



template<typename eT, typename T1>
inline
void
Cube_aux::set_imag(Cube<eT>&, const BaseCube<eT,T1>&)
  {
  arma_extra_debug_sigprint();
  }



template<typename T, typename T1>
inline
void
Cube_aux::set_real(Cube< std::complex<T> >& out, const BaseCube<T,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const ProxyCube<T1> P(X.get_ref());
  
  const uword local_n_rows   = P.get_n_rows();
  const uword local_n_cols   = P.get_n_cols();
  const uword local_n_slices = P.get_n_slices();
  
  arma_debug_assert_same_size
    (
    out.n_rows,   out.n_cols,   out.n_slices,
    local_n_rows, local_n_cols, local_n_slices,
    "Cube::set_real()"
    );
  
  eT* out_mem = out.memptr();
  
  if(ProxyCube<T1>::prefer_at_accessor == false)
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    ea_type A = P.get_ea();
    
    const uword N = out.n_elem;
    
    for(uword i=0; i<N; ++i)
      {
      //out_mem[i].real() = PA[i];
      out_mem[i] = std::complex<T>( A[i], out_mem[i].imag() );
      }
    }
  else
    {
    for(uword slice = 0; slice < local_n_slices; ++slice)
    for(uword col   = 0; col   < local_n_cols;   ++col  )
    for(uword row   = 0; row   < local_n_rows;   ++row  )
      {
      (*out_mem) = std::complex<T>( P.at(row,col,slice), (*out_mem).imag() );
      out_mem++;
      }
    }
  }



template<typename T, typename T1>
inline
void
Cube_aux::set_imag(Cube< std::complex<T> >& out, const BaseCube<T,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const ProxyCube<T1> P(X.get_ref());
  
  const uword local_n_rows   = P.get_n_rows();
  const uword local_n_cols   = P.get_n_cols();
  const uword local_n_slices = P.get_n_slices();
  
  arma_debug_assert_same_size
    (
    out.n_rows,   out.n_cols,   out.n_slices,
    local_n_rows, local_n_cols, local_n_slices,
    "Cube::set_imag()"
    );
  
  eT* out_mem = out.memptr();
  
  if(ProxyCube<T1>::prefer_at_accessor == false)
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    ea_type A = P.get_ea();
    
    const uword N = out.n_elem;
    
    for(uword i=0; i<N; ++i)
      {
      //out_mem[i].imag() = PA[i];
      out_mem[i] = std::complex<T>( out_mem[i].real(), A[i] );
      }
    }
  else
    {
    for(uword slice = 0; slice < local_n_slices; ++slice)
    for(uword col   = 0; col   < local_n_cols;   ++col  )
    for(uword row   = 0; row   < local_n_rows;   ++row  )
      {
      (*out_mem) = std::complex<T>( (*out_mem).real(), P.at(row,col,slice) );
      out_mem++;
      }
    }
  }



#ifdef ARMA_EXTRA_CUBE_MEAT
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_CUBE_MEAT)
#endif



//! @}
