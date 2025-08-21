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


//! \addtogroup cond_rel
//! @{


//
// for preventing pedantic compiler warnings

template<const bool do_eval>
class cond_rel
  {
  public:
  
  template<typename eT> arma_inline static bool lt(const eT A, const eT B);
  template<typename eT> arma_inline static bool gt(const eT A, const eT B);

  template<typename eT> arma_inline static bool leq(const eT A, const eT B);
  template<typename eT> arma_inline static bool geq(const eT A, const eT B);
  
  template<typename eT> arma_inline static eT make_neg(const eT val);
  };



//! @}
