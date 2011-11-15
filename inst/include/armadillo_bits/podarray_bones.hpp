// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup podarray
//! @{



struct podarray_prealloc_n_elem
  {
  static const uword val = 16;
  };



//! A lightweight array for POD types. If the amount of memory requested is small, the stack is used.

template<typename eT>
class podarray
  {
  public:
  
  arma_aligned const uword       n_elem; //!< number of elements held
  arma_aligned const eT* const mem;    //!< pointer to memory used by the object
  
  
  protected:
  //! Internal memory, to avoid calling the 'new' operator for small amounts of memory.
  arma_aligned eT mem_local[ podarray_prealloc_n_elem::val ];
  
  
  public:
  
  inline ~podarray();
  inline  podarray();
  
  inline                 podarray (const podarray& x);
  inline const podarray& operator=(const podarray& x);
  
  arma_inline explicit podarray(const uword new_N);
  
  arma_inline explicit podarray(const eT* X, const uword new_N);
  
  arma_inline eT& operator[] (const uword i);
  arma_inline eT  operator[] (const uword i) const;
  
  arma_inline eT& operator() (const uword i);
  arma_inline eT  operator() (const uword i) const;
  
  inline void set_size(const uword new_n_elem);
  inline void reset();
  
  inline void fill(const eT val);
  
  inline void zeros();
  inline void zeros(const uword new_n_elem);
  
  arma_inline       eT* memptr();
  arma_inline const eT* memptr() const;
  
  arma_hot inline void copy_row(const Mat<eT>& A, const uword row);
  
  protected:
  
  inline void init(const uword new_n_elem);
  
  };



//! @}
