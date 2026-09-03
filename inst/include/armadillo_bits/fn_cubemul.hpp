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


//! \addtogroup fn_cubemul
//! @{



template<typename T1, typename T2>
arma_warn_unused
inline
const GlueCube<T1, T2, glue_cubemul>
cubemul(const BaseCube<typename T1::elem_type,T1>& A, const BaseCube<typename T1::elem_type,T2>& B)
  {
  arma_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cubemul>(A.get_ref(), B.get_ref());
  }



template<typename T1, typename T2>
arma_warn_unused
inline
Cube<typename T1::elem_type>
cubemul(const BaseCube<typename T1::elem_type,T1>& A, const Base<typename T1::elem_type,T2>& B)
  {
  arma_debug_sigprint();
  
  return glue_cubemul::apply(A.get_ref(), B.get_ref());
  }



template<typename T1, typename T2>
arma_warn_unused
inline
Cube<typename T1::elem_type>
cubemul(const Base<typename T1::elem_type,T1>& A, const BaseCube<typename T1::elem_type,T2>& B)
  {
  arma_debug_sigprint();
  
  return glue_cubemul::apply(A.get_ref(), B.get_ref());
  }



//! @}
