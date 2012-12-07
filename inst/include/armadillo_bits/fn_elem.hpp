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


//! \addtogroup fn_elem
//! @{



//
// find

template<typename eT, typename T1>
inline
const mtOp<uword, T1, op_find>
find(const Base<eT,T1>& X, const uword k = 0, const char* direction = "first")
  {
  arma_extra_debug_sigprint();
  
  const char sig = direction[0];
  
  arma_debug_check
    (
    (sig != 'f' && sig != 'F' && sig != 'l' && sig != 'L'),
    "find(): 3rd input argument must be \"first\" or \"last\""
    );
  
  const uword type = (sig == 'f' || sig == 'F') ? 0 : 1;
  
  return mtOp<uword, T1, op_find>(X.get_ref(), k, type);
  }



//
// real

template<typename T1>
arma_inline
const T1&
real(const Base<typename T1::pod_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.get_ref();
  }



template<typename T1>
arma_inline
const T1&
real(const BaseCube<typename T1::pod_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.get_ref();
  }



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_real>
real(const Base<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_real>( X.get_ref() );
  }



template<typename T1>
inline
const mtOpCube<typename T1::pod_type, T1, op_real>
real(const BaseCube<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<typename T1::pod_type, T1, op_real>( X.get_ref() );
  }



//
// imag

template<typename T1>
inline
const Gen< Mat<typename T1::pod_type>, gen_zeros >
imag(const Base<typename T1::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.get_ref());
  
  return Gen< Mat<typename T1::pod_type>, gen_zeros>(A.get_n_rows(), A.get_n_cols());
  }



template<typename T1>
inline
const GenCube<typename T1::pod_type, gen_zeros>
imag(const BaseCube<typename T1::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<T1> A(X.get_ref());
  
  return GenCube<typename T1::pod_type, gen_zeros>(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  }



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_imag>
imag(const Base<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_imag>( X.get_ref() );
  }



template<typename T1>
inline
const mtOpCube<typename T1::pod_type, T1, op_imag>
imag(const BaseCube<std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<typename T1::pod_type, T1, op_imag>( X.get_ref() );
  }



//
// log

template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_log> >::result
log(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_log>(A);
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_log>
log(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_log>(A.get_ref());
  }



//
// log2

template<typename T1>
arma_inline
const eOp<T1, eop_log2>
log2(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_log2>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_log2>
log2(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_log2>(A.get_ref());
  }



//
// log10

template<typename T1>
arma_inline
const eOp<T1, eop_log10>
log10(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_log10>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_log10>
log10(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_log10>(A.get_ref());
  }



//
// exp

template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_exp> >::result
exp(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_exp>(A);
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_exp>
exp(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_exp>(A.get_ref());
  }



// exp2

template<typename T1>
arma_inline
const eOp<T1, eop_exp2>
exp2(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_exp2>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_exp2>
exp2(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_exp2>(A.get_ref());
  }



// exp10

template<typename T1>
arma_inline
const eOp<T1, eop_exp10>
exp10(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_exp10>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_exp10>
exp10(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_exp10>(A.get_ref());
  }



//
// abs


template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_abs> >::result
abs(const T1& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOp<T1, eop_abs>(X);
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_abs>
abs(const BaseCube<typename T1::elem_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOpCube<T1, eop_abs>(X.get_ref());
  }



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_abs>
abs(const Base<std::complex<typename T1::pod_type>, T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return mtOp<typename T1::pod_type, T1, op_abs>( X.get_ref() );
  }



template<typename T1>
inline
const mtOpCube<typename T1::pod_type, T1, op_abs>
abs(const BaseCube< std::complex<typename T1::pod_type>,T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return mtOpCube<typename T1::pod_type, T1, op_abs>( X.get_ref() );
  }



template<typename T1>
arma_inline
const SpOp<T1, spop_abs>
abs(const SpBase<typename T1::elem_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return SpOp<T1, spop_abs>(X.get_ref());
  }



template<typename T1>
arma_inline
const mtSpOp<typename T1::pod_type, T1, spop_cx_abs>
abs(const SpBase< std::complex<typename T1::pod_type>, T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return mtSpOp<typename T1::pod_type, T1, spop_cx_abs>(X.get_ref());
  }



//
// square

template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_square> >::result
square(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_square>(A);
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_square>
square(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_square>(A.get_ref());
  }



template<typename T1>
arma_inline
const SpOp<T1, spop_square>
square(const SpBase<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_square>(A.get_ref());
  }



//
// sqrt

template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_sqrt> >::result
sqrt(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_sqrt>(A);
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_sqrt>
sqrt(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_sqrt>(A.get_ref());
  }



template<typename T1>
arma_inline
const SpOp<T1, spop_sqrt>
sqrt(const SpBase<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_sqrt>(A.get_ref());
  }



//
// conj

template<typename T1>
arma_inline
const T1&
conj(const Base<typename T1::pod_type,T1>& A)
  {
  arma_extra_debug_sigprint();

  return A.get_ref();
  }



template<typename T1>
arma_inline
const T1&
conj(const BaseCube<typename T1::pod_type,T1>& A)
  {
  arma_extra_debug_sigprint();

  return A.get_ref();
  }



template<typename T1>
arma_inline
const eOp<T1, eop_conj>
conj(const Base<std::complex<typename T1::pod_type>,T1>& A)
  {
  arma_extra_debug_sigprint();

  return eOp<T1, eop_conj>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_conj>
conj(const BaseCube<std::complex<typename T1::pod_type>,T1>& A)
  {
  arma_extra_debug_sigprint();

  return eOpCube<T1, eop_conj>(A.get_ref());
  }



template<typename T1>
arma_inline
const T1&
conj(const eOp<T1, eop_conj>& A)
  {
  arma_extra_debug_sigprint();
  
  return A.m;
  }



template<typename T1>
arma_inline
const T1&
conj(const eOpCube<T1, eop_conj>& A)
  {
  arma_extra_debug_sigprint();
  
  return A.m;
  }



// TODO: this needs a more elaborate template restriction mechanism to work properly,
//       i.e. an overloaded version of thus function should do nothing if the input type is non-complex
// 
// //! the conjugate of the transpose of a complex matrix is the same as the hermitian transpose
// template<typename T1>
// arma_inline
// const Op<T1, op_htrans>
// conj(const Op<T1, op_strans>& A)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<T1, op_htrans>(A.m);
//   }



// pow

template<typename T1>
arma_inline
const eOp<T1, eop_pow>
pow(const Base<typename T1::elem_type,T1>& A, const typename T1::elem_type exponent)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_pow>(A.get_ref(), exponent);
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_pow>
pow(const BaseCube<typename T1::elem_type,T1>& A, const typename T1::elem_type exponent)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_pow>(A.get_ref(), exponent);
  }



// pow, specialised handling (non-complex exponent for complex matrices)

template<typename T1>
arma_inline
const eOp<T1, eop_pow>
pow(const Base<typename T1::elem_type,T1>& A, const typename T1::elem_type::value_type exponent)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  return eOp<T1, eop_pow>(A.get_ref(), eT(exponent));
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_pow>
pow(const BaseCube<typename T1::elem_type,T1>& A, const typename T1::elem_type::value_type exponent)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  return eOpCube<T1, eop_pow>(A.get_ref(), eT(exponent));
  }



//
// floor

template<typename T1>
arma_inline
const eOp<T1, eop_floor>
floor(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_floor>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_floor>
floor(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_floor>(A.get_ref());
  }



//
// ceil

template<typename T1>
arma_inline
const eOp<T1, eop_ceil>
ceil(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_ceil>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_ceil>
ceil(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_ceil>(A.get_ref());
  }



//! @}
