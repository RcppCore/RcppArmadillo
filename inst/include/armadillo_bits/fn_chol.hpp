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


//! \addtogroup fn_chol
//! @{



template<typename T1>
inline
const Op<T1, op_chol>
chol(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_chol>(X.get_ref());
  }



template<typename T1>
inline
bool
chol(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  out = chol(X);
  
  return (out.n_elem == 0) ? false : true;
  }



//! @}
