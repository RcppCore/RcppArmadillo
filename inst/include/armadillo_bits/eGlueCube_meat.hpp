// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


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
uword
eGlueCube<T1,T2,eglue_type>::get_n_rows() const
  {
  return P1.get_n_rows();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
uword
eGlueCube<T1,T2,eglue_type>::get_n_cols() const
  {
  return P1.get_n_cols();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
uword
eGlueCube<T1,T2,eglue_type>::get_n_slices() const
  {
  return P1.get_n_slices();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
uword
eGlueCube<T1,T2,eglue_type>::get_n_elem_slice() const
  {
  return P1.get_n_elem_slice();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
uword
eGlueCube<T1,T2,eglue_type>::get_n_elem() const
  {
  return P1.get_n_elem();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
typename T1::elem_type
eGlueCube<T1,T2,eglue_type>::operator[] (const uword i) const
  {
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return P1[i] * P2[i]; }
  }


template<typename T1, typename T2, typename eglue_type>
arma_inline
typename T1::elem_type
eGlueCube<T1,T2,eglue_type>::at(const uword row, const uword col, const uword slice) const
  {
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return P1.at(row,col,slice) + P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return P1.at(row,col,slice) - P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return P1.at(row,col,slice) / P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return P1.at(row,col,slice) * P2.at(row,col,slice); }
  }



//! @}
