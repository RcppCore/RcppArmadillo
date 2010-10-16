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


//! \addtogroup operator_minus
//! @{



//! unary -
template<typename T1>
arma_inline
const eOp<T1, eop_neg>
operator-
(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1,eop_neg>(X.get_ref());
  }



//! cancellation of two consecutive negations: -(-T1)
template<typename T1>
arma_inline
const T1&
operator-
(const eOp<T1, eop_neg>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.m;
  }



//! Base - scalar
template<typename T1>
arma_inline
const eOp<T1, eop_scalar_minus_post>
operator-
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename T1::elem_type           k
  )
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_minus_post>(X.get_ref(), k);
  }



//! scalar - Base
template<typename T1>
arma_inline
const eOp<T1, eop_scalar_minus_pre>
operator-
  (
  const typename T1::elem_type           k,
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_minus_pre>(X.get_ref(), k);
  }



//! complex scalar - non-complex Base (experimental)
template<typename T1>
arma_inline
const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_pre>
operator-
  (
  const std::complex<typename T1::pod_type>& k,
  const Base<typename T1::pod_type, T1>&     X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_pre>('j', X.get_ref(), k);
  }



//! non-complex Base - complex scalar (experimental)
template<typename T1>
arma_inline
const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_post>
operator-
  (
  const Base<typename T1::pod_type, T1>&     X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_post>('j', X.get_ref(), k);
  }



//! subtraction of Base objects with same element type
template<typename T1, typename T2>
arma_inline
const eGlue<T1, T2, eglue_minus>
operator-
  (
  const Base<typename T1::elem_type,T1>& X,
  const Base<typename T1::elem_type,T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return eGlue<T1, T2, eglue_minus>(X.get_ref(), Y.get_ref());
  }



//! subtraction of Base objects with different element types
template<typename T1, typename T2>
inline
const mtGlue<typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_minus>
operator-
  (
  const Base< typename force_different_type<typename T1::elem_type, typename T2::elem_type>::T1_result, T1>& X,
  const Base< typename force_different_type<typename T1::elem_type, typename T2::elem_type>::T2_result, T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  return mtGlue<out_eT, T1, T2, glue_mixed_minus>( X.get_ref(), Y.get_ref() );
  }



//! @}
