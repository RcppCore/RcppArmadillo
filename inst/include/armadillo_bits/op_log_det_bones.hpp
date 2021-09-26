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


//! \addtogroup op_log_det
//! @{



class op_log_det
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static bool apply_direct(typename T1::elem_type& out_val, typename T1::pod_type& out_sign, const Base<typename T1::elem_type,T1>& expr);
  
  template<typename T1>
  inline static bool apply_diagmat(typename T1::elem_type& out_val, typename T1::pod_type& out_sign, const Base<typename T1::elem_type,T1>& expr);
  
  template<typename T1>
  inline static bool apply_trimat(typename T1::elem_type& out_val, typename T1::pod_type& out_sign, const Base<typename T1::elem_type,T1>& expr);
  };



class op_log_det_sympd
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static bool apply_direct(typename T1::pod_type& out_val, const Base<typename T1::elem_type,T1>& expr);
  };



//! @}
