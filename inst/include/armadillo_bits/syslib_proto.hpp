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
    #if !defined(__OPTIMIZE__)
      {
      std::memcpy(dest, src, n_elem*sizeof(eT));
      }
    #else
      {
      switch(n_elem)
        {
        case 0:
          break;
        
        case 1:
          *dest = *src;
          break;
          
        case 2:
          dest[0] = src[0];
          dest[1] = src[1];
          break;
        
        case 3:
          dest[0] = src[0];
          dest[1] = src[1];
          dest[2] = src[2];
          break;
        
        case 4:
          dest[0] = src[0];
          dest[1] = src[1];
          dest[2] = src[2];
          dest[3] = src[3];
          break;
        
        default:
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
      }
    #endif
    }
  
  
  
  template<typename out_eT, typename in_eT>
  arma_hot
  inline
  static
  void
  copy_and_convert_elem(out_eT* dest, const in_eT* src, const u32 n_elem)
    {
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      dest[i] = out_eT( src[i] );
      dest[j] = out_eT( src[j] );
      }
    
    if(i < n_elem)
      {
      dest[i] = out_eT( src[i] );
      }
    }
  
  
  
  //
  // TODO: the functions below will need more work
  
  
  template<typename out_eT, typename in_eT>
  arma_hot
  arma_inline
  static
  void
  convert_cx_scalar(out_eT& out, const in_eT& in)
    {
    out = out_eT(in);
    }
  
  
  
  template<typename out_eT, typename in_T>
  arma_hot
  arma_inline
  static
  void
  convert_cx_scalar(out_eT& out, const std::complex<in_T>& in)
    {
    out = out_eT( in.real() );
    }
  
  
  
  template<typename out_T, typename in_T>
  arma_hot
  arma_inline
  static
  void
  convert_cx_scalar(std::complex<out_T>& out, const std::complex<in_T>& in)
    {
    typedef std::complex<out_T> out_eT;
    
    out = out_eT(in);
    }
  
  
  
  template<typename out_eT, typename in_eT>
  arma_hot
  inline
  static
  void
  copy_and_convert_cx_elem(out_eT* dest, const in_eT* src, const u32 n_elem)
    {
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      convert_cx_scalar( dest[i], src[i] );
      convert_cx_scalar( dest[j], src[j] );
      }
    
    if(i < n_elem)
      {
      convert_cx_scalar( dest[i], src[i] );
      }
    }
  
  
  
  };



//! @}
