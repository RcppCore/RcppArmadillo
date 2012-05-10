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



//! \addtogroup operator_times
//! @{



//! Base * scalar
template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_scalar_times> >::result
operator*
(const T1& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_times>(X,k);
  }



//! scalar * Base
template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_scalar_times> >::result
operator*
(const typename T1::elem_type k, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_times>(X,k);  // NOTE: order is swapped
  }



//! non-complex Base * complex scalar
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_complex<typename T1::elem_type>::value == false),
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>
  >::result
operator*
  (
  const T1&                                  X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>('j', X, k);
  }



//! complex scalar * non-complex Base
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_complex<typename T1::elem_type>::value == false),
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>
  >::result
operator*
  (
  const std::complex<typename T1::pod_type>& k,
  const T1&                                  X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>('j', X, k);
  }



//! scalar * trans(T1)
template<typename T1>
arma_inline
const Op<T1, op_htrans2>
operator*
(const typename T1::elem_type k, const Op<T1, op_htrans>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_htrans2>(X.m, k);
  }



//! trans(T1) * scalar
template<typename T1>
arma_inline
const Op<T1, op_htrans2>
operator*
(const Op<T1, op_htrans>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_htrans2>(X.m, k);
  }



//! Base * diagmat
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const Glue<T1, Op<T2, op_diagmat>, glue_times_diag>
  >::result
operator*
(const T1& X, const Op<T2, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Op<T2, op_diagmat>, glue_times_diag>(X, Y);
  }



//! diagmat * Base
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const Glue<Op<T1, op_diagmat>, T2, glue_times_diag>
  >::result
operator*
(const Op<T1, op_diagmat>& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1, op_diagmat>, T2, glue_times_diag>(X, Y);
  }



//! diagmat * diagmat
template<typename T1, typename T2>
inline
Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
operator*
(const Op<T1, op_diagmat>& X, const Op<T2, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const diagmat_proxy<T1> A(X.m);
  const diagmat_proxy<T2> B(Y.m);
  
  arma_debug_assert_mul_size(A.n_elem, A.n_elem, B.n_elem, B.n_elem, "matrix multiply");
  
  const uword N = A.n_elem;
  
  Mat<out_eT> out(N,N);
  
  out.zeros();
  
  for(uword i=0; i<N; ++i)
    {
    out.at(i,i) = upgrade_val<eT1,eT2>::apply( A[i] ) * upgrade_val<eT1,eT2>::apply( B[i] );
    }
  
  return out;
  }



//! multiplication of Base objects with same element type
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value,
  const Glue<T1, T2, glue_times>
  >::result
operator*
(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X, Y);
  }



//! multiplication of Base objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value == false)),
  const mtGlue< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_times >
  >::result
operator*
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  return mtGlue<out_eT, T1, T2, glue_mixed_times>( X, Y );
  }



//! @}
