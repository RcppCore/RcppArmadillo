// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup operator_cube_relational
//! @{



// <  : lt
// >  : gt
// <= : lteq
// >= : gteq
// == : eq
// != : noteq



template<typename T1, typename T2>
inline
const mtGlueCube<uword, T1, T2, glue_rel_lt>
operator<
(const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlueCube<uword, T1, T2, glue_rel_lt>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlueCube<uword, T1, T2, glue_rel_gt>
operator>
(const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlueCube<uword, T1, T2, glue_rel_gt>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlueCube<uword, T1, T2, glue_rel_lteq>
operator<=
(const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlueCube<uword, T1, T2, glue_rel_lteq>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlueCube<uword, T1, T2, glue_rel_gteq>
operator>=
(const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlueCube<uword, T1, T2, glue_rel_gteq>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlueCube<uword, T1, T2, glue_rel_eq>
operator==
(const BaseCube<typename T1::elem_type,T1>& X, const BaseCube<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlueCube<uword, T1, T2, glue_rel_eq>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlueCube<uword, T1, T2, glue_rel_noteq>
operator!=
(const BaseCube<typename T1::elem_type,T1>& X, const BaseCube<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlueCube<uword, T1, T2, glue_rel_noteq>( X.get_ref(), Y.get_ref() );
  }



//
//
//



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_lt_pre>
operator<
(const typename arma_not_cx<typename T1::elem_type>::result val, const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_lt_pre>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_lt_post>
operator<
(const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result val)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_lt_post>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_gt_pre>
operator>
(const typename arma_not_cx<typename T1::elem_type>::result val, const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_gt_pre>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_gt_post>
operator>
(const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result val)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_gt_post>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_lteq_pre>
operator<=
(const typename arma_not_cx<typename T1::elem_type>::result val, const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_lteq_pre>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_lteq_post>
operator<=
(const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result val)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_lteq_post>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_gteq_pre>
operator>=
(const typename arma_not_cx<typename T1::elem_type>::result val, const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_gteq_pre>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_gteq_post>
operator>=
(const BaseCube<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result val)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_gteq_post>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_eq>
operator==
(const typename T1::elem_type val, const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_eq>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_eq>
operator==
(const BaseCube<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_eq>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_noteq>
operator!=
(const typename T1::elem_type val, const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_noteq>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOpCube<uword, T1, op_rel_noteq>
operator!=
(const BaseCube<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<uword, T1, op_rel_noteq>(X.get_ref(), val);
  }



//! @}
