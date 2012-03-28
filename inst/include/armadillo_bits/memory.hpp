// Copyright (C) 2012 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup memory
//! @{


class memory
  {
  public:
  
  template<typename eT> arma_inline static eT*  acquire(const uword n_elem);
  
  template<typename eT> arma_inline static void release(eT* mem);
  };



template<typename eT>
arma_inline
eT*
memory::acquire(const uword n_elem)
  {
  #if   defined(ARMA_USE_TBB_ALLOC)
    {
    return ( (eT *) scalable_malloc( sizeof(eT)*n_elem) );
    }
  #elif defined(ARMA_USE_MKL_ALLOC)
    {
    return ( (eT *) mkl_malloc( sizeof(eT)*n_elem, 128 ) );
    }
  #else
    {
    return ( new(std::nothrow) eT[n_elem] );
    }
  #endif
  }



template<typename eT>
arma_inline
void
memory::release(eT* mem)
  {
  #if   defined(ARMA_USE_TBB_ALLOC)
    {
    scalable_free( (void *)(mem) );
    }
  #elif defined(ARMA_USE_MKL_ALLOC)
    {
    mkl_free( (void *)(mem) );
    }
  #else
    {
    delete [] mem;
    }
  #endif
  }



//! @}
