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


//! \addtogroup fn_htrans
//! @{



template<typename T1>
arma_inline
const Op<T1, op_htrans>
htrans(const Base<std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_htrans>(X.get_ref());
  }



template<typename T1>
arma_inline
const Op<T1, op_trans>
htrans(const Base<typename T1::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("htrans(): non-complex object given -- using regular transpose");
  
  return Op<T1, op_trans>(X.get_ref());
  }



//! two consecutive hermitian transpose operations cancel each other
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
