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


//! \addtogroup op_elem
//! @{



struct op_real
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_real>& X);
  
  template<typename T1>
  inline static void apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_real>& X);
  };



struct op_imag
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_imag>& X);
  
  template<typename T1>
  inline static void apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_imag>& X);
  };



struct op_abs
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_abs>& X);
  
  template<typename T1>
  inline static void apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_abs>& X);
  };



struct op_arg
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_arg>& X);
  
  template<typename T1>
  inline static void apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_arg>& X);
  };



struct op_replace
  : public traits_op_passthru
  {
  template<typename eT, typename T1> inline static void apply(Mat<eT>& out, const mtOp<eT,T1,op_replace>& in);
  
  template<typename eT, typename T1> inline static void apply(Cube<eT>& out, const mtOpCube<eT,T1,op_replace>& in);
  };



struct op_eps
  : public traits_op_passthru
  {
  template<typename eT> inline static typename get_pod_type<eT>::result direct_eps(const eT& x);
  
  //
  
  template<typename T1>
  inline static void apply(Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_eps>& in);
  
  template<typename T, typename eT>
  inline static void apply_mat_noalias(Mat<T>& out, const Mat<eT>& X);
  
  template<typename T, typename T1>
  inline static void apply_proxy_noalias(Mat<T>& out, const Proxy<T1>& P);
  //
  
  template<typename T1>
  inline static void apply(Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_eps>& in);
  
  template<typename T, typename eT>
  inline static void apply_noalias(Cube<T>& out, const Cube<eT>& X);
  };



//! @}
