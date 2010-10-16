// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


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
u32
eGlue<T1,T2,eglue_type>::get_n_rows() const
  {
  return P1.get_n_rows();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
u32
eGlue<T1,T2,eglue_type>::get_n_cols() const
  {
  return P1.get_n_cols();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
u32
eGlue<T1,T2,eglue_type>::get_n_elem() const
  {
  return P1.get_n_elem();
  }



template<typename T1, typename T2, typename eglue_type>
arma_inline
typename T1::elem_type
eGlue<T1,T2,eglue_type>::operator[] (const u32 i) const
  {
  typedef typename T1::elem_type eT;
  
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return P1[i] * P2[i]; }
  else
    {
    arma_stop("eGlue::operator[]: unhandled eglue_type");
    return eT();
    }
  }


template<typename T1, typename T2, typename eglue_type>
arma_inline
typename T1::elem_type
eGlue<T1,T2,eglue_type>::at(const u32 row, const u32 col) const
  {
  typedef typename T1::elem_type eT;
  
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return P1.at(row,col) + P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return P1.at(row,col) - P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return P1.at(row,col) / P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return P1.at(row,col) * P2.at(row,col); }
  else
    {
    arma_stop("eGlue::at(): unhandled eglue_type");
    return eT();
    }
  }




//! @}
