// Copyright (C) 2012 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup memory
//! @{


class memory
  {
  public:
  
                        arma_inline             static uword enlarge_to_mult_of_chunksize(const uword n_elem);
  
  template<typename eT> arma_inline arma_malloc static eT*   acquire(const uword n_elem);
  
  template<typename eT> arma_inline arma_malloc static eT*   acquire_chunked(const uword n_elem);
  
  template<typename eT> arma_inline             static void  release(eT* mem);
  };



arma_inline
uword
memory::enlarge_to_mult_of_chunksize(const uword n_elem)
  {
  const uword chunksize = arma_config::spmat_chunksize;
  
  // this relies on integer division
  const uword n_elem_mod = ((n_elem % chunksize) != 0) ? ((n_elem / chunksize) + 1) * chunksize : n_elem;
  
  return n_elem_mod;
  }



template<typename eT>
arma_inline
arma_malloc
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



//! get memory in multiples of chunks, holding at least n_elem
template<typename eT>
arma_inline
arma_malloc
eT*
memory::acquire_chunked(const uword n_elem)
  {
  const uword n_elem_mod = memory::enlarge_to_mult_of_chunksize(n_elem);
  
  return memory::acquire<eT>(n_elem_mod);
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
