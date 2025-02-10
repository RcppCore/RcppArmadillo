// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup fn_elem
//! @{


//
// real

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const T1& >::result
real(const T1& X)
  {
  arma_debug_sigprint();
  
  return X;
  }



template<typename T1>
arma_warn_unused
arma_inline
const T1&
real(const BaseCube<typename T1::pod_type, T1>& X)
  {
  arma_debug_sigprint();
  
  return X.get_ref();
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::no), const T1& >::result
real(const T1& X)
  {
  arma_debug_sigprint();
  
  return X;
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::yes), const mtOp<typename T1::pod_type, T1, op_real> >::result
real(const T1& X)
  {
  arma_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_real>( X );
  }



template<typename T1>
arma_warn_unused
inline
const mtOpCube<typename T1::pod_type, T1, op_real>
real(const BaseCube<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_debug_sigprint();
  
  return mtOpCube<typename T1::pod_type, T1, op_real>( X.get_ref() );
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::yes), const mtSpOp<typename T1::pod_type, T1, spop_real> >::result
real(const T1& X)
  {
  arma_debug_sigprint();
  
  return mtSpOp<typename T1::pod_type, T1, spop_real>(X);
  }



//
// imag

template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_arma_type<T1>::value, const mtOp<typename T1::pod_type, T1, op_imag> >::result
imag(const T1& X)
  {
  arma_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_imag>( X );
  }



template<typename T1>
arma_warn_unused
inline
const mtOpCube<typename T1::pod_type, T1, op_imag>
imag(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_debug_sigprint();
  
  return mtOpCube<typename T1::pod_type, T1, op_imag>( X.get_ref() );
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_sparse_type<T1>::value, const mtSpOp<typename T1::pod_type, T1, spop_imag> >::result
imag(const T1& X)
  {
  arma_debug_sigprint();
  
  return mtSpOp<typename T1::pod_type, T1, spop_imag>(X);
  }



//
// log

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_log> >::result
log(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_log>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_log>
log(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_log>(A.get_ref());
  }



//
// log2

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_log2> >::result
log2(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_log2>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_log2>
log2(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_log2>(A.get_ref());
  }



//
// log10

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_log10> >::result
log10(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_log10>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_log10>
log10(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_log10>(A.get_ref());
  }



//
// log1p

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_log1p> >::result
log1p(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_log1p>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const eOpCube<T1, eop_log1p> >::result
log1p(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_log1p>(A.get_ref());
  }



//
// exp

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_exp> >::result
exp(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_exp>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_exp>
exp(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_exp>(A.get_ref());
  }



// exp2

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_exp2> >::result
exp2(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_exp2>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_exp2>
exp2(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_exp2>(A.get_ref());
  }



// exp10

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_exp10> >::result
exp10(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_exp10>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_exp10>
exp10(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_exp10>(A.get_ref());
  }



// expm1

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_expm1> >::result
expm1(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_expm1>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const eOpCube<T1, eop_expm1> >::result
expm1(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_expm1>(A.get_ref());
  }



//
// abs


template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_abs> >::result
abs(const T1& X)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_abs>(X);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_abs>
abs(const BaseCube<typename T1::elem_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = nullptr)
  {
  arma_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOpCube<T1, eop_abs>(X.get_ref());
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::yes), const mtOp<typename T1::pod_type, T1, op_abs> >::result
abs(const T1& X)
  {
  arma_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_abs>(X);
  }



template<typename T1>
arma_warn_unused
inline
const mtOpCube<typename T1::pod_type, T1, op_abs>
abs(const BaseCube< std::complex<typename T1::pod_type>,T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr)
  {
  arma_debug_sigprint();
  
  arma_ignore(junk);
  
  return mtOpCube<typename T1::pod_type, T1, op_abs>( X.get_ref() );
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::no), const SpOp<T1, spop_abs> >::result
abs(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_abs>(X);
  }




template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::yes), const mtSpOp<typename T1::pod_type, T1, spop_cx_abs> >::result
abs(const T1& X)
  {
  arma_debug_sigprint();
  
  return mtSpOp<typename T1::pod_type, T1, spop_cx_abs>(X);
  }



//
// arg


template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_arg> >::result
arg(const T1& X)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_arg>(X);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_arg>
arg(const BaseCube<typename T1::elem_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = nullptr)
  {
  arma_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOpCube<T1, eop_arg>(X.get_ref());
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::yes), const mtOp<typename T1::pod_type, T1, op_arg> >::result
arg(const T1& X)
  {
  arma_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_arg>(X);
  }



template<typename T1>
arma_warn_unused
inline
const mtOpCube<typename T1::pod_type, T1, op_arg>
arg(const BaseCube< std::complex<typename T1::pod_type>,T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr)
  {
  arma_debug_sigprint();
  
  arma_ignore(junk);
  
  return mtOpCube<typename T1::pod_type, T1, op_arg>( X.get_ref() );
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::no), const SpOp<T1, spop_arg> >::result
arg(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_arg>(X);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::yes), const mtSpOp<typename T1::pod_type, T1, spop_cx_arg> >::result
arg(const T1& X)
  {
  arma_debug_sigprint();
  
  return mtSpOp<typename T1::pod_type, T1, spop_cx_arg>(X);
  }



//
// square

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_square> >::result
square(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_square>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_square>
square(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_square>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpOp<T1, spop_square> >::result
square(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_square>(X);
  }



//
// sqrt

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_sqrt> >::result
sqrt(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_sqrt>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_sqrt>
sqrt(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_sqrt>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpOp<T1, spop_sqrt> >::result
sqrt(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_sqrt>(X);
  }



//
// cbrt

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_cbrt> >::result
cbrt(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_cbrt>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const eOpCube<T1, eop_cbrt> >::result
cbrt(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_cbrt>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::no), const SpOp<T1, spop_cbrt> >::result
cbrt(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_cbrt>(X);
  }



//
// conj

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const T1& >::result
conj(const T1& X)
  {
  arma_debug_sigprint();
  
  return X;
  }



template<typename T1>
arma_warn_unused
arma_inline
const T1&
conj(const BaseCube<typename T1::pod_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return A.get_ref();
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::no), const T1& >::result
conj(const T1& X)
  {
  arma_debug_sigprint();
  
  return X;
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::yes), const eOp<T1, eop_conj> >::result
conj(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_conj>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_conj>
conj(const BaseCube<std::complex<typename T1::pod_type>,T1>& A)
  {
  arma_debug_sigprint();

  return eOpCube<T1, eop_conj>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::yes), const SpOp<T1, spop_conj> >::result
conj(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_conj>(X);
  }



// pow

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_pow> >::result
pow(const T1& A, const typename T1::elem_type exponent)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_pow>(A, exponent);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_pow>
pow(const BaseCube<typename T1::elem_type,T1>& A, const typename T1::elem_type exponent)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_pow>(A.get_ref(), exponent);
  }



// pow, specialised handling (non-complex exponent for complex matrices)

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::yes), const eOp<T1, eop_pow> >::result
pow(const T1& A, const typename T1::elem_type::value_type exponent)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  return eOp<T1, eop_pow>(A, eT(exponent));
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_pow>
pow(const BaseCube<std::complex<typename T1::pod_type>,T1>& A, const typename T1::elem_type::value_type exponent)
  {
  arma_debug_sigprint();
  
  typedef std::complex<typename T1::pod_type> eT;
  
  return eOpCube<T1, eop_pow>(A.get_ref(), eT(exponent));
  }



//
// floor

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_floor> >::result
floor(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_floor>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_floor>
floor(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_floor>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpOp<T1, spop_floor> >::result
floor(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_floor>(X);
  }



//
// ceil

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_ceil> >::result
ceil(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_ceil>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_ceil>
ceil(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_ceil>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpOp<T1, spop_ceil> >::result
ceil(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_ceil>(X);
  }



//
// round

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_round> >::result
round(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_round>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_round>
round(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_round>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpOp<T1, spop_round> >::result
round(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_round>(X);
  }



//
// trunc

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_trunc> >::result
trunc(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_trunc>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_trunc>
trunc(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_trunc>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpOp<T1, spop_trunc> >::result
trunc(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_trunc>(X);
  }



//
// sign

template<typename eT>
arma_warn_unused
arma_inline
typename arma_scalar_only<eT>::result
sign(const eT x)
  {
  arma_debug_sigprint();
  
  return arma_sign(x);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_sign> >::result
sign(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_sign>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
const eOpCube<T1, eop_sign>
sign(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_sign>(A.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpOp<T1, spop_sign> >::result
sign(const T1& X)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_sign>(X);
  }



//
// erf

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_erf> >::result
erf(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_erf>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const eOpCube<T1, eop_erf> >::result
erf(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_erf>(A.get_ref());
  }



//
// erfc

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_erfc> >::result
erfc(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_erfc>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const eOpCube<T1, eop_erfc> >::result
erfc(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_erfc>(A.get_ref());
  }



//
// lgamma

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_lgamma> >::result
lgamma(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_lgamma>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const eOpCube<T1, eop_lgamma> >::result
lgamma(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_lgamma>(A.get_ref());
  }



//
// tgamma

template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no), const eOp<T1, eop_tgamma> >::result
tgamma(const T1& A)
  {
  arma_debug_sigprint();
  
  return eOp<T1, eop_tgamma>(A);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const eOpCube<T1, eop_tgamma> >::result
tgamma(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_debug_sigprint();
  
  return eOpCube<T1, eop_tgamma>(A.get_ref());
  }



// the functions below are currently unused; reserved for potential future use

template<typename T1> void exp_approx(const T1&) { arma_stop_logic_error("unimplemented"); }
template<typename T1> void log_approx(const T1&) { arma_stop_logic_error("unimplemented"); }
template<typename T1> void approx_exp(const T1&) { arma_stop_logic_error("unimplemented"); }
template<typename T1> void approx_log(const T1&) { arma_stop_logic_error("unimplemented"); }

//! @}
