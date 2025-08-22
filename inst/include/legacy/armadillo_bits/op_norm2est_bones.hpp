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


//! \addtogroup op_norm2est
//! @{



template<typename eT>
struct norm2est_randu_filler
  {
  std::mt19937_64                    local_engine;
  std::uniform_real_distribution<eT> local_u_distr;
  
  inline norm2est_randu_filler();
  
  inline void fill(eT* mem, const uword N);
  };


template<typename T>
struct norm2est_randu_filler< std::complex<T> >
  {
  std::mt19937_64                   local_engine;
  std::uniform_real_distribution<T> local_u_distr;
  
  inline norm2est_randu_filler();
  
  inline void fill(std::complex<T>* mem, const uword N);
  };



class op_norm2est
  : public traits_op_default
  {
  public:
  
  template<typename T1> inline static typename T1::pod_type norm2est(const   Base<typename T1::elem_type, T1>& X, const typename T1::pod_type tolerance, const uword max_iter);
  template<typename T1> inline static typename T1::pod_type norm2est(const SpBase<typename T1::elem_type, T1>& X, const typename T1::pod_type tolerance, const uword max_iter);
  };



//! @}
