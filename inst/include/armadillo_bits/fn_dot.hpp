// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_dot
//! @{


template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename T1::elem_type
dot
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B
  )
  {
  arma_extra_debug_sigprint();
  
  return op_dot::apply(A,B);
  }



template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename T1::elem_type
norm_dot
  (
  const Base<typename T1::elem_type,T1>& A, 
  const Base<typename T1::elem_type,T2>& B,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return op_norm_dot::apply(A,B);
  }



//
// cdot



template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename T1::elem_type
cdot
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const typename arma_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return op_cdot::apply(A,B);
  }



template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename T1::elem_type
cdot
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return op_dot::apply(A,B);
  }




// convert dot(htrans(x), y) to cdot(x,y)

template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename T1::elem_type
dot
  (
  const Op<T1, op_htrans>& A,
  const Base<typename T1::elem_type,T2>& B,
  const typename arma_cx_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return cdot(A.m, B);
  }



//! @}
