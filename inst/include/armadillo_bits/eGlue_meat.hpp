// Copyright (C) 2010-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup eGlue
//! @{



template<typename T1, typename T2, typename eglue_type>
arma_inline
eGlue<T1,T2,eglue_type>::~eGlue()
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
eGlue<T1,T2,eglue_type>::eGlue(const T1& in_A, const T2& in_B)
  : P1(in_A)
  , P2(in_B)
  {
  arma_extra_debug_sigprint();
  
  // arma_debug_assert_same_size( P1, P2, eglue_type::text() );
  arma_debug_assert_same_size
    (
    P1.get_n_rows(), P1.get_n_cols(),
    P2.get_n_rows(), P2.get_n_cols(),
    eglue_type::text()
    );
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
uword
eGlue<T1,T2,eglue_type>::get_n_rows() const
  {
  return is_row ? 1 : ( Proxy<T1>::is_fixed ? P1.get_n_rows() : ( Proxy<T2>::is_fixed ? P2.get_n_rows() : P1.get_n_rows() ) );
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
uword
eGlue<T1,T2,eglue_type>::get_n_cols() const
  {
  return is_col ? 1 : ( Proxy<T1>::is_fixed ? P1.get_n_cols() : ( Proxy<T2>::is_fixed ? P2.get_n_cols() : P1.get_n_cols() ) );
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
uword
eGlue<T1,T2,eglue_type>::get_n_elem() const
  {
  return Proxy<T1>::is_fixed ? P1.get_n_elem() : ( Proxy<T2>::is_fixed ? P2.get_n_elem() : P1.get_n_elem() ) ;
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
typename T1::elem_type
eGlue<T1,T2,eglue_type>::operator[] (const uword ii) const
  {
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return P1[ii] + P2[ii]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return P1[ii] - P2[ii]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return P1[ii] / P2[ii]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return P1[ii] * P2[ii]; }
  }


template<typename T1, typename T2, typename eglue_type>
arma_inline
typename T1::elem_type
eGlue<T1,T2,eglue_type>::at(const uword row, const uword col) const
  {
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return P1.at(row,col) + P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return P1.at(row,col) - P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return P1.at(row,col) / P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return P1.at(row,col) * P2.at(row,col); }
  }



//! @}
