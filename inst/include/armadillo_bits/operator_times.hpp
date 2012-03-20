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
const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>
operator*
  (
  const Base<typename T1::pod_type, T1>&     X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>('j', X.get_ref(), k);
  }



//! complex scalar * non-complex Base (experimental)
template<typename T1>
arma_inline
const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>
operator*
  (
  const std::complex<typename T1::pod_type>& k,
  const Base<typename T1::pod_type, T1>&     X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>('j', X.get_ref(), k);
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
const Glue<T1, T2, glue_times>
operator*
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



//! multiplication of Base objects with different element types
template<typename T1, typename T2>
inline
const mtGlue< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_times >
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
  
  return mtGlue<out_eT, T1, T2, glue_mixed_times>( X.get_ref(), Y.get_ref() );
  }






//
//
//

#undef armaObj
#undef armaObjA
#undef armaObjB

// TODO: add handling of scalar*trans(object)
// TODO: add handling of trans(object)*scalar

#undef arma_operator_obj_scalar
#undef arma_operator_scalar_obj
#undef arma_operator_obj_cx_scalar
#undef arma_operator_cx_scalar_obj
#undef arma_operator_obj_base
#undef arma_operator_base_obj
#undef arma_operator_obj_diagmat
#undef arma_operator_diagmat_obj
#undef arma_operator_obj1_obj2
#undef arma_operator_obj1_obj2_mixed



#define arma_operator_obj_scalar(armaObj) \
template<typename eT> \
arma_inline \
const eOp<armaObj<eT>, eop_scalar_times> \
operator* \
(const armaObj<eT>& X, const typename armaObj<eT>::elem_type k) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return eOp<armaObj<eT>, eop_scalar_times>(X, k); \
  }



#define arma_operator_scalar_obj(armaObj) \
template<typename eT> \
arma_inline \
const eOp<armaObj<eT>, eop_scalar_times> \
operator* \
(const typename armaObj<eT>::elem_type k, const armaObj<eT>& X) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return eOp<armaObj<eT>, eop_scalar_times>(X, k); \
  }



#define arma_operator_obj_cx_scalar(armaObj) \
template<typename eT> \
arma_inline \
const mtOp<typename std::complex<typename armaObj<eT>::pod_type>, armaObj<eT>, op_cx_scalar_times> \
operator* \
  ( \
  const armaObj<eT>&                                  X, \
  const std::complex<typename armaObj<eT>::pod_type>& k  \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return mtOp<typename std::complex<typename armaObj<eT>::pod_type>, armaObj<eT>, op_cx_scalar_times>('j', X, k); \
  }



#define arma_operator_cx_scalar_obj(armaObj) \
template<typename eT> \
arma_inline \
const mtOp<typename std::complex<typename armaObj<eT>::pod_type>, armaObj<eT>, op_cx_scalar_times> \
operator* \
  ( \
  const std::complex<typename armaObj<eT>::pod_type>& k, \
  const armaObj<eT>&                                  X \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return mtOp<typename std::complex<typename armaObj<eT>::pod_type>, armaObj<eT>, op_cx_scalar_times>('j', X, k); \
  }



#define arma_operator_obj_base(armaObj) \
template<typename BT> \
arma_inline \
const Glue<armaObj<typename BT::elem_type>, BT, glue_times> \
operator* \
  ( \
  const armaObj<typename BT::elem_type>&     X, \
  const    Base<typename BT::elem_type, BT>& Y \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return Glue<armaObj<typename BT::elem_type>, BT, glue_times>(X, Y.get_ref()); \
  }



#define arma_operator_base_obj(armaObj) \
template<typename BT> \
arma_inline \
const Glue<BT, armaObj<typename BT::elem_type>, glue_times> \
operator* \
  ( \
  const    Base<typename BT::elem_type, BT>& X, \
  const armaObj<typename BT::elem_type>    & Y  \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return Glue<BT, armaObj<typename BT::elem_type>, glue_times>(X.get_ref(), Y); \
  }



#define arma_operator_obj_diagmat(armaObj) \
template<typename T1> \
arma_inline \
const Glue< armaObj<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag > \
operator* \
  ( \
  const armaObj<typename T1::elem_type>& X, \
  const Op<T1, op_diagmat>&              Y  \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return Glue< armaObj<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag >(X, Y); \
  }



#define arma_operator_diagmat_obj(armaObj) \
template<typename T1> \
arma_inline \
const Glue< Op<T1, op_diagmat>, armaObj<typename T1::elem_type>, glue_times_diag > \
operator* \
  ( \
  const Op<T1, op_diagmat>&              X, \
  const armaObj<typename T1::elem_type>& Y  \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return Glue< Op<T1, op_diagmat>, armaObj<typename T1::elem_type>, glue_times_diag >(X, Y); \
  }



#define arma_operator_obj1_obj2(armaObjA, armaObjB) \
template<typename eT> \
arma_inline \
const Glue<armaObjA<eT>, armaObjB<eT>, glue_times> \
operator* \
  ( \
  const armaObjA<eT>& X, \
  const armaObjB<eT>& Y \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return Glue<armaObjA<eT>, armaObjB<eT>, glue_times>(X, Y); \
  }



#define arma_operator_obj1_obj2_mixed(armaObjA, armaObjB) \
template<typename eT1, typename eT2> \
inline \
const mtGlue<typename promote_type<eT1, eT2>::result, armaObjA<eT1>, armaObjB<eT2>, glue_mixed_times> \
operator* \
  ( \
  const armaObjA<eT1>& X, \
  const armaObjB<eT2>& Y  \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  typedef typename promote_type<eT1,eT2>::result out_eT; \
  \
  promote_type<eT1,eT2>::check(); \
  \
  return mtGlue<out_eT, armaObjA<eT1>, armaObjB<eT2>, glue_mixed_times>( X, Y ); \
  }



arma_operator_obj_scalar(Col)
arma_operator_obj_scalar(Row)
arma_operator_obj_scalar(diagview)
arma_operator_obj_scalar(subview_col)
arma_operator_obj_scalar(subview_row)
                         
arma_operator_scalar_obj(Col)
arma_operator_scalar_obj(Row)
arma_operator_scalar_obj(diagview)
arma_operator_scalar_obj(subview_col)
arma_operator_scalar_obj(subview_row)

arma_operator_obj_cx_scalar(Col)
arma_operator_obj_cx_scalar(Row)
arma_operator_obj_cx_scalar(diagview)
arma_operator_obj_cx_scalar(subview_col)
arma_operator_obj_cx_scalar(subview_row)
                            
arma_operator_cx_scalar_obj(Col)
arma_operator_cx_scalar_obj(Row)
arma_operator_cx_scalar_obj(diagview)
arma_operator_cx_scalar_obj(subview_col)
arma_operator_cx_scalar_obj(subview_row)

arma_operator_obj_base(Col)
arma_operator_obj_base(Row)
arma_operator_obj_base(diagview)
arma_operator_obj_base(subview_col)
arma_operator_obj_base(subview_row)

arma_operator_base_obj(Col)
arma_operator_base_obj(Row)
arma_operator_base_obj(diagview)
arma_operator_base_obj(subview_col)
arma_operator_base_obj(subview_row)

arma_operator_obj_diagmat(Col)
arma_operator_obj_diagmat(Row)
arma_operator_obj_diagmat(diagview)
arma_operator_obj_diagmat(subview_col)
arma_operator_obj_diagmat(subview_row)

arma_operator_diagmat_obj(Col)
arma_operator_diagmat_obj(Row)
arma_operator_diagmat_obj(diagview)
arma_operator_diagmat_obj(subview_col)
arma_operator_diagmat_obj(subview_row)


arma_operator_obj1_obj2(Col,Col)
arma_operator_obj1_obj2(Col,Row)
arma_operator_obj1_obj2(Col,diagview)
arma_operator_obj1_obj2(Col,subview_col)
arma_operator_obj1_obj2(Col,subview_row)

arma_operator_obj1_obj2(Row,Col)
arma_operator_obj1_obj2(Row,Row)
arma_operator_obj1_obj2(Row,diagview)
arma_operator_obj1_obj2(Row,subview_col)
arma_operator_obj1_obj2(Row,subview_row)

arma_operator_obj1_obj2(diagview,Col)
arma_operator_obj1_obj2(diagview,Row)
arma_operator_obj1_obj2(diagview,diagview)
arma_operator_obj1_obj2(diagview,subview_col)
arma_operator_obj1_obj2(diagview,subview_row)

arma_operator_obj1_obj2(subview_col,Col)
arma_operator_obj1_obj2(subview_col,Row)
arma_operator_obj1_obj2(subview_col,diagview)
arma_operator_obj1_obj2(subview_col,subview_col)
arma_operator_obj1_obj2(subview_col,subview_row)

arma_operator_obj1_obj2(subview_row,Col)
arma_operator_obj1_obj2(subview_row,Row)
arma_operator_obj1_obj2(subview_row,diagview)
arma_operator_obj1_obj2(subview_row,subview_col)
arma_operator_obj1_obj2(subview_row,subview_row)



arma_operator_obj1_obj2_mixed(Col,Col)
arma_operator_obj1_obj2_mixed(Col,Row)
arma_operator_obj1_obj2_mixed(Col,diagview)
arma_operator_obj1_obj2_mixed(Col,subview_col)
arma_operator_obj1_obj2_mixed(Col,subview_row)

arma_operator_obj1_obj2_mixed(Row,Col)
arma_operator_obj1_obj2_mixed(Row,Row)
arma_operator_obj1_obj2_mixed(Row,diagview)
arma_operator_obj1_obj2_mixed(Row,subview_col)
arma_operator_obj1_obj2_mixed(Row,subview_row)

arma_operator_obj1_obj2_mixed(diagview,Col)
arma_operator_obj1_obj2_mixed(diagview,Row)
arma_operator_obj1_obj2_mixed(diagview,diagview)
arma_operator_obj1_obj2_mixed(diagview,subview_col)
arma_operator_obj1_obj2_mixed(diagview,subview_row)

arma_operator_obj1_obj2_mixed(subview_col,Col)
arma_operator_obj1_obj2_mixed(subview_col,Row)
arma_operator_obj1_obj2_mixed(subview_col,diagview)
arma_operator_obj1_obj2_mixed(subview_col,subview_col)
arma_operator_obj1_obj2_mixed(subview_col,subview_row)
arma_operator_obj1_obj2_mixed(subview_row,Col)
arma_operator_obj1_obj2_mixed(subview_row,Row)
arma_operator_obj1_obj2_mixed(subview_row,diagview)
arma_operator_obj1_obj2_mixed(subview_row,subview_col)
arma_operator_obj1_obj2_mixed(subview_row,subview_row)



// TODO: explicit handling of subview_elem1, perhaps via expanding the above macros to take another paramater (for template<T1>)


//! @}
