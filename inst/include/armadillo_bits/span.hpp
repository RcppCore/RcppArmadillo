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



//! \addtogroup span
//! @{


template<typename T> struct span_base;


template<typename T>
struct span_base
  {
  const uword a;
  const uword b;
  const bool  whole;
  
  constexpr          span_base()                                   : a(0),    b(0),    whole(true)  {}
  constexpr explicit span_base(const uword in_a)                   : a(in_a), b(in_a), whole(false) {}
  constexpr explicit span_base(const uword in_a, const uword in_b) : a(in_a), b(in_b), whole(false) {}
  
  static const span_base<T> all;
  };


template<typename T>
const span_base<T> span_base<T>::all = span_base<T>();

typedef span_base<int> span;


//! @}
