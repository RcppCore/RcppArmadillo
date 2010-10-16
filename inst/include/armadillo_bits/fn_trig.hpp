// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_trig
//! @{

//
// trigonometric functions:
// cos family: cos, acos, cosh, acosh
// sin family: sin, asin, sinh, asinh
// tan family: tan, atan, tanh, atanh


//
// cos

template<typename T1>
arma_inline
const eOp<T1, eop_cos>
cos(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_cos>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cos>
cos(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cos>(A.get_ref());
  }



//
// acos

template<typename T1>
arma_inline
const eOp<T1, eop_acos>
acos(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_acos>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_acos>
acos(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_acos>(A.get_ref());
  }



//
// cosh

template<typename T1>
arma_inline
const eOp<T1, eop_cosh>
cosh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_cosh>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cosh>
cosh(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cosh>(A.get_ref());
  }



//
// acosh

template<typename T1>
arma_inline
const eOp<T1, eop_acosh>
acosh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_acosh>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_acosh>
acosh(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_acosh>(A.get_ref());
  }



//
// sin

template<typename T1>
arma_inline
const eOp<T1, eop_sin>
sin(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_sin>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_sin>
sin(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_sin>(A.get_ref());
  }



//
// asin

template<typename T1>
arma_inline
const eOp<T1, eop_asin>
asin(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_asin>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_asin>
asin(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_asin>(A.get_ref());
  }



//
// sinh

template<typename T1>
arma_inline
const eOp<T1, eop_sinh>
sinh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_sinh>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_sinh>
sinh(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_sinh>(A.get_ref());
  }



//
// asinh

template<typename T1>
arma_inline
const eOp<T1, eop_asinh>
asinh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_asinh>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_asinh>
asinh(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_asinh>(A.get_ref());
  }



//
// tan

template<typename T1>
arma_inline
const eOp<T1, eop_tan>
tan(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_tan>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_tan>
tan(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_tan>(A.get_ref());
  }



//
// atan

template<typename T1>
arma_inline
const eOp<T1, eop_atan>
atan(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_atan>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_atan>
atan(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_atan>(A.get_ref());
  }



//
// tanh

template<typename T1>
arma_inline
const eOp<T1, eop_tanh>
tanh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_tanh>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_tanh>
tanh(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_tanh>(A.get_ref());
  }



//
// atanh

template<typename T1>
arma_inline
const eOp<T1, eop_atanh>
atanh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_atanh>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_atanh>
atanh(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_atanh>(A.get_ref());
  }



//! @}
