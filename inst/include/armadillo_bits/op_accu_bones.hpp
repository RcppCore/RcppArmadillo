// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup op_accu
//! @{



struct  op_accu_mat
  : public traits_op_passthru
  {
  template<typename T1>
  static inline typename T1::elem_type apply_proxy_linear(const Proxy<T1>& P);
  
  template<typename T1>
  static inline typename T1::elem_type apply_proxy_at(const Proxy<T1>& P);
  
  template<typename T1>
  static inline typename T1::elem_type apply(const T1& X);
  
  template<typename T1, typename functor>
  static inline typename T1::elem_type apply_omit_helper(const Proxy<T1>& P, functor is_omitted);
  
  template<typename T1>
  static inline typename T1::elem_type apply(const Op<T1, op_omit>& in);
  
  template<typename T1>
  static inline typename T1::elem_type apply(const eOp<T1,eop_square>& expr);
  
  template<typename T1>
  static inline typename T1::elem_type apply(const eOp<T1,eop_pow>& expr);
  
  template<typename T1, typename T2>
  static inline typename T1::elem_type apply(const eGlue<T1,T2,eglue_schur>& expr);
  
  template<typename T1, typename op_type>
  static inline uword apply(const mtOp<uword,T1,op_type>& X, const typename arma_op_rel_only<op_type>::result* junk1 = nullptr, const typename arma_not_cx<typename T1::elem_type>::result* junk2 = nullptr);
  
  template<typename T1, typename op_type>
  static inline uword apply(const mtOp<uword,T1,op_type>& X, const typename arma_op_rel_only<op_type>::result* junk1 = nullptr, const typename arma_cx_only<typename T1::elem_type>::result* junk2 = nullptr);
  
  template<typename T1, typename T2>
  static inline uword apply(const mtGlue<uword,T1,T2,glue_rel_noteq>& X);
  
  template<typename T1, typename T2>
  static inline uword apply(const mtGlue<uword,T1,T2,glue_rel_eq>& X);
  
  template<typename eT>
  static inline eT apply(const subview<eT>& X);
  
  template<typename eT>
  static inline eT apply(const subview_col<eT>& X);
  
  template<typename eT>
  static inline eT apply(const subview_row<eT>& X);
  };



struct  op_accu_cube
  : public traits_op_passthru
  {
  template<typename T1>
  static inline typename T1::elem_type apply_proxy_linear(const ProxyCube<T1>& P);
  
  template<typename T1>
  static inline typename T1::elem_type apply_proxy_at(const ProxyCube<T1>& P);
  
  template<typename T1>
  static inline typename T1::elem_type apply(const BaseCube<typename T1::elem_type,T1>& X);
  
  template<typename T1>
  static inline typename T1::elem_type apply(const eOpCube<T1,eop_square>& expr);
  
  template<typename T1>
  static inline typename T1::elem_type apply(const eOpCube<T1,eop_pow>& expr);
  
  template<typename T1, typename T2>
  static inline typename T1::elem_type apply(const eGlueCube<T1,T2,eglue_schur>& expr);
  
  template<typename T1, typename functor>
  static inline typename T1::elem_type apply_omit_helper(const ProxyCube<T1>& P, functor is_omitted);
  
  template<typename T1>
  static inline typename T1::elem_type apply(const CubeToMatOp<T1, op_omit_cube>& in);
  
  template<typename eT>
  static inline eT apply(const subview_cube<eT>& sv);
  };



//! @}
