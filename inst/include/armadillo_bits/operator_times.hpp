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



//! \addtogroup operator_times
//! @{



//! Base * scalar
template<typename T1>
arma_inline
const eOp<T1, eop_scalar_times>
operator*
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_times>(X.get_ref(),k);
  }



//! scalar * Base
template<typename T1>
arma_inline
const eOp<T1, eop_scalar_times>
operator*
(const typename T1::elem_type k, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_times>(X.get_ref(),k);  // NOTE: order is swapped
  }



//! non-complex Base * complex scalar (experimental)
template<typename T1>
arma_inline
Mat<typename std::complex<typename T1::pod_type> >
operator*
  (
  const Base<typename T1::pod_type, T1>&     X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.get_ref());
  
  Mat<eT> out(A.n_rows, A.n_cols);
  
  const u32 n_elem  = A.n_elem;
        eT* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] * k;
    }
  
  return out;
  }



//! complex scalar * non-complex Base (experimental)
template<typename T1>
arma_inline
Mat<typename std::complex<typename T1::pod_type> >
operator*
  (
  const std::complex<typename T1::pod_type>& k,
  const Base<typename T1::pod_type, T1>&     X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.get_ref());
  
  Mat<eT> out(A.n_rows, A.n_cols);
  
  const u32 n_elem  = A.n_elem;
        eT* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = k * A[i];
    }
  
  return out;
  }



//! scalar * trans(T1)
template<typename T1>
arma_inline
const Op<T1, op_trans2>
operator*
(const typename T1::elem_type k, const Op<T1, op_trans>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_trans2>(X.m, k);
  }



//! trans(T1) * scalar
template<typename T1>
arma_inline
const Op<T1, op_trans2>
operator*
(const Op<T1, op_trans>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_trans2>(X.m, k);
  }



//! Base * diagmat
template<typename T1, typename T2>
arma_inline
const Glue<T1, Op<T2, op_diagmat>, glue_times_diag>
operator*
(const Base<typename T2::elem_type,T1>& X, const Op<T2, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Op<T2, op_diagmat>, glue_times_diag>(X.get_ref(), Y);
  }



//! diagmat * Base
template<typename T1, typename T2>
arma_inline
const Glue<Op<T1, op_diagmat>, T2, glue_times_diag>
operator*
(const Op<T1, op_diagmat>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1, op_diagmat>, T2, glue_times_diag>(X, Y.get_ref());
  }



//! diagmat * diagmat
template<typename T1, typename T2>
arma_inline
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
  
  const u32 N = A.n_elem;
  
  Mat<out_eT> out(N,N);
  
  out.zeros();
  
  for(u32 i=0; i<N; ++i)
    {
    out.at(i,i) = upgrade_val<eT1,eT2>::apply( A[i] ) * upgrade_val<eT1,eT2>::apply( B[i] );
    }
  
  return out;
  }



//! multiplication of Base objects with same element type
template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



//! multiplication of Base objects with different element types
template<typename T1, typename T2>
arma_inline
Mat<typename promote_type<typename T1::elem_type, typename T2::elem_type>::result>
operator*
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
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT1>& A = tmp1.M;
  const Mat<eT2>& B = tmp2.M;
  
  Mat<out_eT> out;
  
  glue_times::apply_mixed(out, A, B);
  
  return out;
  }



//! @}
