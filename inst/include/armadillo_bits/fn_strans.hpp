// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_strans
//! @{



template<typename T1>
arma_inline
const Op<T1, op_strans>
strans(const Base<typename T1::elem_type,T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_strans>(X.get_ref());
  }



template<typename T1>
arma_inline
const Op<T1, op_htrans>
strans(const Base<typename T1::elem_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_htrans>(X.get_ref());
  }



//! two consecutive transpose operations cancel each other
template<typename T1>
arma_inline
const T1&
strans(const Op<T1, op_strans>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("strans(): removing op_strans");
  
  return X.m;
  }



//! @}
