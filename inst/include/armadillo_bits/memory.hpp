// Copyright (C) 2012-2013 NICTA (www.nicta.com.au)
// Copyright (C) 2012-2013 Conrad Sanderson
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
  
  
  template<typename eT> arma_inline static bool      is_aligned(const eT*  mem);
  template<typename eT> arma_inline static void mark_as_aligned(      eT*& mem);
  template<typename eT> arma_inline static void mark_as_aligned(const eT*& mem);
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
    return ( (eT *) scalable_malloc(sizeof(eT)*n_elem) );
    }
  #elif defined(ARMA_USE_MKL_ALLOC)
    {
    return ( (eT *) mkl_malloc( sizeof(eT)*n_elem, 128 ) );
    }
  #elif defined(ARMA_HAVE_POSIX_MEMALIGN)
    {
    eT* memptr;
    
    const size_t alignment = 16;  // change the 16 to 64 if you wish to align to the cache line
    
    int status = posix_memalign((void **)&memptr, ( (alignment >= sizeof(void*)) ? alignment : sizeof(void*) ), sizeof(eT)*n_elem);
    
    return (status == 0) ? memptr : NULL;
    }
  #elif defined(_MSC_VER)
    {
    return ( (eT *) _aligned_malloc( sizeof(eT)*n_elem, 16 ) );  // lives in malloc.h
    }
  #else
    {
    //return ( new(std::nothrow) eT[n_elem] );
    return ( (eT *) malloc(sizeof(eT)*n_elem) );
    }
  #endif
  
  // TODO: for mingw, use __mingw_aligned_malloc
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
  #elif defined(ARMA_HAVE_POSIX_MEMALIGN)
    {
    free( (void *)(mem) );
    }
  #elif defined(_MSC_VER)
    {
    _aligned_free( (void *)(mem) );
    }
  #else
    {
    //delete [] mem;
    free( (void *)(mem) );
    }
  #endif
  
  // TODO: for mingw, use __mingw_aligned_free
  }



template<typename eT>
arma_inline
bool
memory::is_aligned(const eT* mem)
  {
  #if (defined(ARMA_HAVE_ICC_ASSUME_ALIGNED) || defined(ARMA_HAVE_GCC_ASSUME_ALIGNED)) && !defined(ARMA_DONT_CHECK_ALIGNMENT)
    {
    return ((std::ptrdiff_t(mem) & 0x0F) == 0);
    }
  #else
    {
    arma_ignore(mem);
    
    return false;
    }
  #endif
  }



template<typename eT>
arma_inline
void
memory::mark_as_aligned(eT*& mem)
  {
  #if defined(ARMA_HAVE_ICC_ASSUME_ALIGNED)
    {
    __assume_aligned(mem, 16);
    }
  #elif defined(ARMA_HAVE_GCC_ASSUME_ALIGNED)
    {
    mem = (eT*)__builtin_assume_aligned(mem, 16);
    }
  #else
    {
    arma_ignore(mem);
    }
  #endif
  
  // TODO: MSVC?  __assume( (mem & 0x0F) == 0 );
  //
  // http://comments.gmane.org/gmane.comp.gcc.patches/239430
  // GCC __builtin_assume_aligned is similar to ICC's __assume_aligned,
  // so for lvalue first argument ICC's __assume_aligned can be emulated using
  // #define __assume_aligned(lvalueptr, align) lvalueptr = __builtin_assume_aligned (lvalueptr, align) 
  //
  // http://www.inf.ethz.ch/personal/markusp/teaching/263-2300-ETH-spring11/slides/class19.pdf
  // http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/cpp/lin/index.htm
  // http://d3f8ykwhia686p.cloudfront.net/1live/intel/CompilerAutovectorizationGuide.pdf
  }



template<typename eT>
arma_inline
void
memory::mark_as_aligned(const eT*& mem)
  {
  #if defined(ARMA_HAVE_ICC_ASSUME_ALIGNED)
    {
    __assume_aligned(mem, 16);
    }
  #elif defined(ARMA_HAVE_GCC_ASSUME_ALIGNED)
    {
    mem = (const eT*)__builtin_assume_aligned(mem, 16);
    }
  #else
    {
    arma_ignore(mem);
    }
  #endif
  }



//! @}
