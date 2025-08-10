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


//! \addtogroup op_norm2est
//! @{



template<typename eT>
struct norm2est_randu_filler
  {
  typedef typename promote_type<eT, float>::result eT_promoted;
  
  std::mt19937_64                             local_engine;
  std::uniform_real_distribution<eT_promoted> local_u_distr;
  
  inline norm2est_randu_filler();
  
  inline void fill(eT* mem, const uword N);
  };


template<typename T>
struct norm2est_randu_filler< std::complex<T> >
  {
  typedef typename promote_type<T, float>::result T_promoted;
  
  std::mt19937_64                            local_engine;
  std::uniform_real_distribution<T_promoted> local_u_distr;
  
  inline norm2est_randu_filler();
  
  inline void fill(std::complex<T>* mem, const uword N);
  };



struct op_norm2est
  : public traits_op_default
  {
  template<typename T1> inline static typename T1::pod_type norm2est(const   Base<typename T1::elem_type, T1>& X, const typename T1::pod_type tolerance, const uword max_iter);
  template<typename T1> inline static typename T1::pod_type norm2est(const SpBase<typename T1::elem_type, T1>& X, const typename T1::pod_type tolerance, const uword max_iter);
  };



//! @}
