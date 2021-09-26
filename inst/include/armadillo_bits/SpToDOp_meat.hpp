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


//! \addtogroup SpToDOp
//! @{



template<typename T1, typename op_type>
inline
SpToDOp<T1, op_type>::SpToDOp(const T1& in_m)
  : m(in_m)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename op_type>
inline
SpToDOp<T1, op_type>::SpToDOp(const T1& in_m, const typename T1::elem_type in_aux)
  : m(in_m)
  , aux(in_aux)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename op_type>
inline
SpToDOp<T1, op_type>::~SpToDOp()
  {
  arma_extra_debug_sigprint();
  }



//! @}
