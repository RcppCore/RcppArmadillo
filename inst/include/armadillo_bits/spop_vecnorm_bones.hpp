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


//! \addtogroup spop_vecnorm
//! @{


class spop_vecnorm
  : public traits_op_xvec
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type,T1,spop_vecnorm>& expr);
  
  template<typename eT>
  inline static void apply_direct(Mat< typename get_pod_type<eT>::result >& out, const SpMat<eT>& X, const uword k);
  };


//


class spop_vecnorm_ext
  : public traits_op_xvec
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type,T1,spop_vecnorm_ext>& expr);
  
  template<typename eT>
  inline static void apply_direct(Mat< typename get_pod_type<eT>::result >& out, const SpMat<eT>& X, const uword method_id);
  };


//! @}
