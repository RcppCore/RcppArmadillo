// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_trans
//! @{


template<typename T1>
arma_inline
const Op<T1, op_htrans>
trans
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename enable_if<is_basevec<T1>::value == false>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_htrans>(X.get_ref());
  }



template<typename T1>
arma_inline
arma_warn_unused
const Op<T1, op_htrans>
trans
  (
  const T1& X,
  const typename enable_if<is_basevec<T1>::value == true>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_htrans>(X);
  }



template<typename T1>
arma_inline
const Op<T1, op_htrans>
htrans
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename enable_if<is_basevec<T1>::value == false>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_htrans>(X.get_ref());
  }



template<typename T1>
arma_inline
arma_warn_unused
const Op<T1, op_htrans>
htrans
  (
  const T1& X,
  const typename enable_if<is_basevec<T1>::value == true>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_htrans>(X);
  }



//! two consecutive transpose operations cancel each other
template<typename T1>
arma_inline
const T1&
trans(const Op<T1, op_htrans>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("trans(): removing op_htrans");
  
  return X.m;
  }



template<typename T1>
arma_inline
const T1&
htrans(const Op<T1, op_htrans>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("htrans(): removing op_htrans");
  
  return X.m;
  }



//! @}
