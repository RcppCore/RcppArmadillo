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


//! \addtogroup spop_accu
//! @{



struct  spop_accu
  : public traits_op_passthru
  {
  template<typename T1>
  static inline typename T1::elem_type apply(const SpBase<typename T1::elem_type,T1>& expr);

  template<typename T1, typename T2>
  static inline typename T1::elem_type apply(const SpGlue<T1,T2,spglue_plus>& expr);

  template<typename T1, typename T2>
  static inline typename T1::elem_type apply(const SpGlue<T1,T2,spglue_minus>& expr);

  template<typename T1, typename T2>
  static inline typename T1::elem_type apply(const SpGlue<T1,T2,spglue_schur>& expr);

  template<typename T1, typename spop_type>
  static inline typename T1::elem_type apply(const SpOp<T1, spop_type>& expr);

  template<typename T1>
  static inline typename T1::elem_type apply(const SpOp<T1, spop_square>& expr);

  template<typename T1, typename functor>
  static inline typename T1::elem_type apply_spop_omit_helper(const T1& expr, functor is_omitted);

  template<typename T1>
  static inline typename T1::elem_type apply(const SpOp<T1, spop_omit>& expr);

  template<typename T1, typename spop_type>
  static inline uword apply(const mtSpOp<uword,T1,spop_type>& X, const typename arma_spop_rel_only<spop_type>::result* junk1 = nullptr, const typename arma_not_cx<typename T1::elem_type>::result* junk2 = nullptr);

  template<typename T1, typename spop_type>
  static inline uword apply(const mtSpOp<uword,T1,spop_type>& X, const typename arma_spop_rel_only<spop_type>::result* junk1 = nullptr, const typename arma_cx_only<typename T1::elem_type>::result* junk2 = nullptr);
  };



//! @}
