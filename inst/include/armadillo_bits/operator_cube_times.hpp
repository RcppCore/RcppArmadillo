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




//! \addtogroup operator_cube_times
//! @{



//! BaseCube * scalar
template<typename T1>
arma_inline
const eOpCube<T1, eop_scalar_times>
operator*
  (
  const BaseCube<typename T1::elem_type,T1>& X,
  const typename T1::elem_type               k
  )
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_scalar_times>(X.get_ref(), k);
  }



//! scalar * BaseCube
template<typename T1>
arma_inline
const eOpCube<T1, eop_scalar_times>
operator*
  (
  const typename T1::elem_type               k,
  const BaseCube<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_scalar_times>(X.get_ref(), k);
  }



//! non-complex BaseCube * complex scalar (experimental)
template<typename T1>
arma_inline
const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>
operator*
  (
  const BaseCube<typename T1::pod_type, T1>& X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>('j', X.get_ref(), k);
  }



//! complex scalar * non-complex BaseCube (experimental)
template<typename T1>
arma_inline
const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>
operator*
  (
  const std::complex<typename T1::pod_type>& k,
  const BaseCube<typename T1::pod_type, T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>('j', X.get_ref(), k);
  }



//! @}
