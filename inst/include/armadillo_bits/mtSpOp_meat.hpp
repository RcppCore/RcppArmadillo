// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup mtSpOp
//! @{



template<typename out_eT, typename T1, typename op_type>
inline
mtSpOp<out_eT, T1, op_type>::mtSpOp(const T1& in_m)
  : m(in_m)
  {
  arma_extra_debug_sigprint();
  }



template<typename out_eT, typename T1, typename op_type>
inline
mtSpOp<out_eT, T1, op_type>::mtSpOp(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b)
  : m(in_m)
  , aux_uword_a(in_aux_uword_a)
  , aux_uword_b(in_aux_uword_b)
  {
  arma_extra_debug_sigprint();
  }



template<typename out_eT, typename T1, typename op_type>
inline
mtSpOp<out_eT, T1, op_type>::mtSpOp(const char junk, const T1& in_m, const out_eT in_aux)
  : m(in_m)
  , aux_out_eT(in_aux)
  {
  arma_ignore(junk);
  
  arma_extra_debug_sigprint();
  }



template<typename out_eT, typename T1, typename op_type>
inline
mtSpOp<out_eT, T1, op_type>::~mtSpOp()
  {
  arma_extra_debug_sigprint();
  }



template<typename out_eT, typename T1, typename op_type>
template<typename eT2>
arma_inline
bool
mtSpOp<out_eT, T1, op_type>::is_alias(const SpMat<eT2>& X) const
  {
  return (void_ptr(&X) == void_ptr(&m));
  }



//! @}
