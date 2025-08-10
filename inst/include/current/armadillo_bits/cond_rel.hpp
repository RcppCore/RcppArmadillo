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


//! \addtogroup cond_rel
//! @{


//
// for preventing pedantic compiler warnings

template<bool do_eval>
struct cond_rel
  {
  template<typename eT> static constexpr bool lt(const eT A, const eT B);
  template<typename eT> static constexpr bool gt(const eT A, const eT B);

  template<typename eT> static constexpr bool leq(const eT A, const eT B);
  template<typename eT> static constexpr bool geq(const eT A, const eT B);
  
  template<typename eT> static constexpr eT make_neg(const eT val);
  };



template<>
template<typename eT>
constexpr
bool
cond_rel<true>::lt(const eT A, const eT B)
  {
  return (A < B);
  }
  


template<>
template<typename eT>
constexpr
bool
cond_rel<false>::lt(const eT, const eT)
  {
  return false;
  }
  


template<>
template<typename eT>
constexpr
bool
cond_rel<true>::gt(const eT A, const eT B)
  {
  return (A > B);
  }
  


template<>
template<typename eT>
constexpr
bool
cond_rel<false>::gt(const eT, const eT)
  {
  return false;
  }
  


template<>
template<typename eT>
constexpr
bool
cond_rel<true>::leq(const eT A, const eT B)
  {
  return (A <= B);
  }
  


template<>
template<typename eT>
constexpr
bool
cond_rel<false>::leq(const eT, const eT)
  {
  return false;
  }
  


template<>
template<typename eT>
constexpr
bool
cond_rel<true>::geq(const eT A, const eT B)
  {
  return (A >= B);
  }
  


template<>
template<typename eT>
constexpr
bool
cond_rel<false>::geq(const eT, const eT)
  {
  return false;
  }



template<>
template<typename eT>
constexpr
eT
cond_rel<true>::make_neg(const eT val)
  {
  return -val;
  }
  


template<>
template<typename eT>
constexpr
eT
cond_rel<false>::make_neg(const eT)
  {
  return eT(0);
  }
  


//! @}
