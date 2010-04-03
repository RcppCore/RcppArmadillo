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


//! \addtogroup syslib
//! @{


class syslib
  {
  public:
  
  template<typename eT>
  arma_hot
  inline
  static
  void
  copy_elem(eT* dest, const eT* src, const u32 n_elem)
    {
    if( n_elem <= (128/sizeof(eT)) )
      {
      u32 i,j;
      
      for(i=0, j=1; j<n_elem; i+=2, j+=2)
        {
        dest[i] = src[i];
        dest[j] = src[j];
        }
      
      if(i < n_elem)
        {
        dest[i] = src[i];
        }
      }
    else
      {
      std::memcpy(dest, src, n_elem*sizeof(eT));
      }
    
    }
  
  };



//! @}
