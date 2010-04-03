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


//! \addtogroup operator_cube_minus
//! @{



//! unary -
template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_neg>
operator-
  (
  const BaseCube<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_neg>(X.get_ref());
  }



//! cancellation of two consecutive negations: -(-T1)
template<typename T1>
arma_inline
const T1&
operator-
  (
  const eOpCube<T1, eop_cube_neg>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return X.m;
  }



//! BaseCube - scalar
template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_scalar_minus_post>
operator-
  (
  const BaseCube<typename T1::elem_type,T1>& X,
  const typename T1::elem_type               k
  )
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_scalar_minus_post>(X.get_ref(), k);
  }



//! scalar - BaseCube
template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_scalar_minus_pre>
operator-
  (
  const typename T1::elem_type               k,
  const BaseCube<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_scalar_minus_pre>(X.get_ref(), k);
  }



//! subtraction of BaseCube objects with same element type
template<typename T1, typename T2>
arma_inline
const eGlueCube<T1, T2, eglue_cube_minus>
operator-
  (
  const BaseCube<typename T1::elem_type,T1>& X,
  const BaseCube<typename T1::elem_type,T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return eGlueCube<T1, T2, eglue_cube_minus>(X.get_ref(), Y.get_ref());
  }



//! subtraction of Base objects with different element types
template<typename T1, typename T2>
arma_inline
Cube<typename promote_type<typename T1::elem_type, typename T2::elem_type>::result>
operator-
  (
  const BaseCube< typename force_different_type<typename T1::elem_type, typename T2::elem_type>::T1_result, T1>& X,
  const BaseCube< typename force_different_type<typename T1::elem_type, typename T2::elem_type>::T2_result, T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.get_ref());
  const ProxyCube<T2> B(Y.get_ref());
  
  arma_debug_assert_same_size(A, B, "cube subtraction");
  
  Cube<out_eT> out(A.n_rows, A.n_cols, A.n_slices);

        out_eT* out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) - upgrade_val<eT1,eT2>::apply(B[i]);
    }
  
  return out;
  }



//! @}
