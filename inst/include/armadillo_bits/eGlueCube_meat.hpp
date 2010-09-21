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


//! \addtogroup eGlueCube
//! @{



template<typename T1, typename T2, typename eglue_type>
arma_inline
eGlueCube<T1,T2,eglue_type>::~eGlueCube()
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
eGlueCube<T1,T2,eglue_type>::eGlueCube(const T1& in_A, const T2& in_B)
  : P1(in_A)
  , P2(in_B)
  {
  arma_extra_debug_sigprint();
  
  arma_assert_same_size
    (
    P1.get_n_rows(), P1.get_n_cols(), P1.get_n_slices(),
    P2.get_n_rows(), P2.get_n_cols(), P2.get_n_slices(), 
    eglue_type::text()
    );
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
u32
eGlueCube<T1,T2,eglue_type>::get_n_rows() const
  {
  return P1.get_n_rows();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
u32
eGlueCube<T1,T2,eglue_type>::get_n_cols() const
  {
  return P1.get_n_cols();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
u32
eGlueCube<T1,T2,eglue_type>::get_n_slices() const
  {
  return P1.get_n_slices();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
u32
eGlueCube<T1,T2,eglue_type>::get_n_elem_slice() const
  {
  return P1.get_n_elem_slice();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
u32
eGlueCube<T1,T2,eglue_type>::get_n_elem() const
  {
  return P1.get_n_elem();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
typename T1::elem_type
eGlueCube<T1,T2,eglue_type>::operator[] (const u32 i) const
  {
  typedef typename T1::elem_type eT;
  
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return P1[i] * P2[i]; }
  else
    {
    arma_stop("eGlueCube::operator[]: unhandled eglue_type");
    return eT();
    }
  }


template<typename T1, typename T2, typename eglue_type>
arma_inline
typename T1::elem_type
eGlueCube<T1,T2,eglue_type>::at(const u32 row, const u32 col, const u32 slice) const
  {
  typedef typename T1::elem_type eT;
  
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return P1.at(row,col,slice) + P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return P1.at(row,col,slice) - P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return P1.at(row,col,slice) / P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return P1.at(row,col,slice) * P2.at(row,col,slice); }
  else
    {
    arma_stop("eGlueCube::at(): unhandled eglue_type");
    return eT();
    }
  }



//! @}
