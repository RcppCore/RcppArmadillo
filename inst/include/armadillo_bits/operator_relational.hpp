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


//! \addtogroup operator_relational
//! @{


// <  : lt
// >  : gt
// <= : lteq
// >= : gteq
// == : eq
// != : noteq



template<typename T1, typename T2>
inline
const mtGlue<uword, T1, T2, glue_rel_lt>
operator<
(const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const Base<typename arma_not_cx<typename T1::elem_type>::result,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_lt>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlue<uword, T1, T2, glue_rel_gt>
operator>
(const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const Base<typename arma_not_cx<typename T1::elem_type>::result,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_gt>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlue<uword, T1, T2, glue_rel_lteq>
operator<=
(const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const Base<typename arma_not_cx<typename T1::elem_type>::result,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_lteq>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlue<uword, T1, T2, glue_rel_gteq>
operator>=
(const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const Base<typename arma_not_cx<typename T1::elem_type>::result,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_gteq>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlue<uword, T1, T2, glue_rel_eq>
operator==
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_eq>( X.get_ref(), Y.get_ref() );
  }



template<typename T1, typename T2>
inline
const mtGlue<uword, T1, T2, glue_rel_noteq>
operator!=
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword, T1, T2, glue_rel_noteq>( X.get_ref(), Y.get_ref() );
  }



//
//
//



template<typename T1>
inline
const mtOp<uword, T1, op_rel_lt_pre>
operator<
(const typename arma_not_cx<typename T1::elem_type>::result val, const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_lt_pre>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_lt_post>
operator<
(const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_lt_post>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_gt_pre>
operator>
(const typename arma_not_cx<typename T1::elem_type>::result val, const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_gt_pre>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_gt_post>
operator>
(const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_gt_post>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_lteq_pre>
operator<=
(const typename arma_not_cx<typename T1::elem_type>::result val, const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_lteq_pre>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_lteq_post>
operator<=
(const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_lteq_post>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_gteq_pre>
operator>=
(const typename arma_not_cx<typename T1::elem_type>::result val, const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_gteq_pre>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_gteq_post>
operator>=
(const Base<typename arma_not_cx<typename T1::elem_type>::result,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_gteq_post>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_eq>
operator==
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_eq>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_eq>
operator==
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_eq>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_noteq>
operator!=
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_noteq>(X.get_ref(), val);
  }



template<typename T1>
inline
const mtOp<uword, T1, op_rel_noteq>
operator!=
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword, T1, op_rel_noteq>(X.get_ref(), val);
  }



//! @}
