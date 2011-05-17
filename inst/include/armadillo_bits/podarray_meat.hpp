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


template<typename eT>
inline
podarray<eT>::~podarray()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(n_elem > sizeof(mem_local)/sizeof(eT) )
    {
    delete [] mem;
    }
  
  if(arma_config::debug == true)
    {
    access::rw(n_elem) = 0;
    access::rw(mem)    = 0;
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
  : n_elem(0)
  , mem   (0)
  {
  arma_extra_debug_sigprint();
  
  this->operator=(x);
  }
  
  
  
template<typename eT>
inline
const podarray<eT>&
podarray<eT>::operator=(const podarray& x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    init(x.n_elem);
    
    arrayops::copy( memptr(), x.memptr(), n_elem );
    }
  
  return *this;
  }



template<typename eT>
arma_inline
podarray<eT>::podarray(const u32 new_n_elem)
  : n_elem(0)
  , mem   (0)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(new_n_elem);
  }



template<typename eT>
arma_inline
podarray<eT>::podarray(const eT* X, const u32 new_n_elem)
  : n_elem(0)
  , mem   (0)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(new_n_elem);
  
  arrayops::copy( memptr(), X, new_n_elem );
  }



template<typename eT>
arma_inline
eT
podarray<eT>::operator[] (const u32 i) const
  {
  return mem[i];
  }



template<typename eT>
arma_inline
eT&
podarray<eT>::operator[] (const u32 i)
  {
  return access::rw(mem[i]);
  }
  
  
  
template<typename eT>
arma_inline
eT
podarray<eT>::operator() (const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "podarray::operator(): index out of bounds");
  return mem[i];
  }



template<typename eT>
arma_inline
eT&
podarray<eT>::operator() (const u32 i)
  {
  arma_debug_check( (i >= n_elem), "podarray::operator(): index out of bounds");
  return access::rw(mem[i]);
  }



template<typename eT>
inline
void
podarray<eT>::set_size(const u32 new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  init(new_n_elem);
  }



template<typename eT>
inline
void
podarray<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  init(0);
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
  
  fill(eT(0));
  }



template<typename eT>
inline
void
podarray<eT>::zeros(const u32 new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  init(new_n_elem);
  fill(eT(0));
  }



template<typename eT>
arma_inline
eT*
podarray<eT>::memptr()
  {
  return const_cast<eT*>(mem);
  }
  
  

template<typename eT>
arma_inline
const eT*
podarray<eT>::memptr() const
  {
  return mem;
  }



template<typename eT>
arma_hot
inline
void
podarray<eT>::copy_row(const Mat<eT>& A, const u32 row)
  {
  const u32 cols = A.n_cols;
  
  // note: this function assumes that the podarray has been set to the correct size beforehand
  eT* out = memptr();
  
  switch(cols)
    {
    default:
      {
      u32 i,j;
      for(i=0, j=1; j < cols; i+=2, j+=2)
        {
        out[i] = A.at(row, i);
        out[j] = A.at(row, j);
        }
      
      if(i < cols)
        {
        out[i] = A.at(row, i);
        }
      }
      break;
    
    case 8:
      out[7] = A.at(row, 7);

    case 7:
      out[6] = A.at(row, 6);

    case 6:
      out[5] = A.at(row, 5);

    case 5:
      out[4] = A.at(row, 4);

    case 4:
      out[3] = A.at(row, 3);

    case 3:
      out[2] = A.at(row, 2);

    case 2:
      out[1] = A.at(row, 1);

    case 1:
      out[0] = A.at(row, 0);
    }
  }
  


template<typename eT>
inline
void
podarray<eT>::init(const u32 new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  if(n_elem == new_n_elem)
    {
    return;
    }
    
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
    access::rw(mem) = new eT[new_n_elem];
    }
  
  access::rw(n_elem) = new_n_elem;
  }



//! @}
