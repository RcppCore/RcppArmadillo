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


//! \addtogroup op_sp_stddev
//! @{



class op_sp_stddev
  : public traits_op_xvec
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::pod_type>& out, const mtSpReduceOp<typename T1::pod_type, T1, op_sp_stddev>& in);
  
  template<typename T1>
  inline static void apply_slow(Mat<typename T1::pod_type>& out, const SpProxy<T1>& p, const uword norm_type, const uword dim);
  
  template<typename T1>
  inline static typename T1::pod_type stddev_vec(const T1& X, const uword norm_type = 0);
  };



//! @}
