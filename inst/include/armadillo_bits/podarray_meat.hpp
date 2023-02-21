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


//! \addtogroup podarray
//! @{


template<typename eT>
inline
podarray<eT>::~podarray()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(n_elem > podarray_prealloc_n_elem::val )
    {
    memory::release( mem );
    }
  }



template<typename eT>
inline
podarray<eT>::podarray()
  : n_elem(0)
  , mem   (0)
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
inline
podarray<eT>::podarray(const podarray& x)
  : n_elem(x.n_elem)
  {
  arma_extra_debug_sigprint();
  
  const uword x_n_elem = x.n_elem;
  
  init_cold(x_n_elem);
  
  arrayops::copy( memptr(), x.memptr(), x_n_elem );
  }



template<typename eT>
inline
const podarray<eT>&
podarray<eT>::operator=(const podarray& x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    const uword x_n_elem = x.n_elem;
    
    init_warm(x_n_elem);
    
    arrayops::copy( memptr(), x.memptr(), x_n_elem );
    }
  
  return *this;
  }



template<typename eT>
arma_inline
podarray<eT>::podarray(const uword new_n_elem)
  : n_elem(new_n_elem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold(new_n_elem);
  }



template<typename eT>
template<bool do_zeros>
inline
podarray<eT>::podarray(const uword new_n_elem, const arma_initmode_indicator<do_zeros>&)
  : n_elem(new_n_elem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold(new_n_elem);
  
  if(do_zeros)
    {
    arma_extra_debug_print("podarray::constructor: zeroing memory");
    arrayops::fill_zeros(memptr(), n_elem);
    }
  }



template<typename eT>
arma_inline
eT
podarray<eT>::operator[] (const uword i) const
  {
  return mem[i];
  }



template<typename eT>
arma_inline
eT&
podarray<eT>::operator[] (const uword i)
  {
  return access::rw(mem[i]);
  }



template<typename eT>
arma_inline
eT
podarray<eT>::operator() (const uword i) const
  {
  arma_debug_check_bounds( (i >= n_elem), "podarray::operator(): index out of bounds" );
  
  return mem[i];
  }



template<typename eT>
arma_inline
eT&
podarray<eT>::operator() (const uword i)
  {
  arma_debug_check_bounds( (i >= n_elem), "podarray::operator(): index out of bounds" );
  
  return access::rw(mem[i]);
  }



template<typename eT>
inline
void
podarray<eT>::set_min_size(const uword min_n_elem)
  {
  arma_extra_debug_sigprint();
  
  if(min_n_elem > n_elem)  { init_warm(min_n_elem); }
  }



template<typename eT>
inline
void
podarray<eT>::set_size(const uword new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  init_warm(new_n_elem);
  }



template<typename eT>
inline
void
podarray<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  init_warm(0);
  }



template<typename eT>
inline
void
podarray<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arrayops::inplace_set(memptr(), val, n_elem);
  }



template<typename eT>
inline
void
podarray<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  arrayops::fill_zeros(memptr(), n_elem);
  }



template<typename eT>
inline
void
podarray<eT>::zeros(const uword new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  init_warm(new_n_elem);
  
  arrayops::fill_zeros(memptr(), n_elem);
  }



template<typename eT>
arma_inline
eT*
podarray<eT>::memptr()
  {
  return mem;
  }



template<typename eT>
arma_inline
const eT*
podarray<eT>::memptr() const
  {
  return mem;
  }



template<typename eT>
inline
void
podarray<eT>::copy_row(const Mat<eT>& A, const uword row)
  {
  arma_extra_debug_sigprint();
  
  // note: this function assumes that the podarray has been set to the correct size beforehand
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  const eT*   A_mem = &(A.at(row,0));
        eT* out_mem = memptr();
  
  for(uword i=0; i < n_cols; ++i)
    {
    out_mem[i] = (*A_mem);
    
    A_mem += n_rows;
    }
  }



template<typename eT>
inline
void
podarray<eT>::init_cold(const uword new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  mem = (new_n_elem <= podarray_prealloc_n_elem::val) ? mem_local : memory::acquire<eT>(new_n_elem);
  }



template<typename eT>
inline
void
podarray<eT>::init_warm(const uword new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  if(n_elem == new_n_elem)  { return; }
    
  if(n_elem > podarray_prealloc_n_elem::val)  { memory::release( mem ); }
  
  mem = (new_n_elem <= podarray_prealloc_n_elem::val) ? mem_local : memory::acquire<eT>(new_n_elem);
  
  access::rw(n_elem) = new_n_elem;
  }



//! @}
