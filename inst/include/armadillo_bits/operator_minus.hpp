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



//
//
//

#undef armaObj
#undef armaObjA
#undef armaObjB

#undef arma_operator_unary
#undef arma_operator_obj_scalar
#undef arma_operator_scalar_obj
#undef arma_operator_obj_cx_scalar
#undef arma_operator_cx_scalar_obj
#undef arma_operator_obj_base
#undef arma_operator_base_obj
#undef arma_operator_obj1_obj2
#undef arma_operator_obj1_obj2_mixed




#define arma_operator_unary(armaObj) \
template<typename eT> \
arma_inline \
const eOp<armaObj<eT>, eop_neg> \
operator- \
(const armaObj<eT>& X) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return eOp<armaObj<eT>,eop_neg>(X); \
  }



#define arma_operator_obj_scalar(armaObj) \
template<typename eT> \
arma_inline \
const eOp<armaObj<eT>, eop_scalar_minus_post> \
operator- \
(const armaObj<eT>& X, const typename armaObj<eT>::elem_type k) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return eOp<armaObj<eT>, eop_scalar_minus_post>(X, k); \
  }



#define arma_operator_scalar_obj(armaObj) \
template<typename eT> \
arma_inline \
const eOp<armaObj<eT>, eop_scalar_minus_pre> \
operator- \
(const typename armaObj<eT>::elem_type k, const armaObj<eT>& X) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return eOp<armaObj<eT>, eop_scalar_minus_pre>(X, k); \
  }



#define arma_operator_obj_cx_scalar(armaObj) \
template<typename eT> \
arma_inline \
const mtOp<typename std::complex<typename armaObj<eT>::pod_type>, armaObj<eT>, op_cx_scalar_minus_post> \
operator- \
  ( \
  const armaObj<eT>&                                  X, \
  const std::complex<typename armaObj<eT>::pod_type>& k  \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return mtOp<typename std::complex<typename armaObj<eT>::pod_type>, armaObj<eT>, op_cx_scalar_minus_post>('j', X, k); \
  }



#define arma_operator_cx_scalar_obj(armaObj) \
template<typename eT> \
arma_inline \
const mtOp<typename std::complex<typename armaObj<eT>::pod_type>, armaObj<eT>, op_cx_scalar_minus_pre> \
operator- \
  ( \
  const std::complex<typename armaObj<eT>::pod_type>& k, \
  const armaObj<eT>&                                  X \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return mtOp<typename std::complex<typename armaObj<eT>::pod_type>, armaObj<eT>, op_cx_scalar_minus_pre>('j', X, k); \
  }



#define arma_operator_obj_base(armaObj) \
template<typename BT> \
arma_inline \
const eGlue<armaObj<typename BT::elem_type>, BT, eglue_minus> \
operator- \
  ( \
  const armaObj<typename BT::elem_type>&     X, \
  const    Base<typename BT::elem_type, BT>& Y \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return eGlue<armaObj<typename BT::elem_type>, BT, eglue_minus>(X, Y.get_ref()); \
  }



#define arma_operator_base_obj(armaObj) \
template<typename BT> \
arma_inline \
const eGlue<BT, armaObj<typename BT::elem_type>, eglue_minus> \
operator- \
  ( \
  const    Base<typename BT::elem_type, BT>& X, \
  const armaObj<typename BT::elem_type>    & Y  \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return eGlue<BT, armaObj<typename BT::elem_type>, eglue_minus>(X.get_ref(), Y); \
  }



#define arma_operator_obj1_obj2(armaObjA, armaObjB) \
template<typename eT> \
arma_inline \
const eGlue<armaObjA<eT>, armaObjB<eT>, eglue_minus> \
operator- \
  ( \
  const armaObjA<eT>& X, \
  const armaObjB<eT>& Y \
  ) \
  { \
  arma_extra_debug_sigprint(); \
  \
  return eGlue<armaObjA<eT>, armaObjB<eT>, eglue_minus>(X, Y); \
  }



#define arma_operator_obj1_obj2_mixed(armaObjA, armaObjB) \
template<typename eT1, typename eT2> \
inline \
const mtGlue<typename promote_type<eT1, eT2>::result, armaObjA<eT1>, armaObjB<eT2>, glue_mixed_minus> \
operator- \
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
  return mtGlue<out_eT, armaObjA<eT1>, armaObjB<eT2>, glue_mixed_minus>( X, Y ); \
  }



arma_operator_unary(Col)
arma_operator_unary(Row)
arma_operator_unary(diagview)
arma_operator_unary(subview_col)
arma_operator_unary(subview_row)

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
