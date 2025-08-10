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


//! \addtogroup spop_misc
//! @{


struct spop_scalar_times
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_scalar_times>& in);
  };



struct spop_cx_scalar_times
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat< std::complex<typename T1::pod_type> >& out, const mtSpOp< std::complex<typename T1::pod_type>, T1, spop_cx_scalar_times>& in);
  };



struct spop_square
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_square>& in);
  };



struct spop_sqrt
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_sqrt>& in);
  };



struct spop_cbrt
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_cbrt>& in);
  };



struct spop_abs
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_abs>& in);
  };



struct spop_cx_abs
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_cx_abs>& in);
  };



struct spop_arg
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_arg>& in);
  };



struct spop_cx_arg
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_cx_arg>& in);
  };



struct spop_real
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_real>& in);
  };



struct spop_imag
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_imag>& in);
  };



struct spop_conj
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_conj>& in);
  };



struct spop_repelem
  : public traits_op_default
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_repelem>& in);
  };



struct spop_reshape
  : public traits_op_default
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_reshape>& in);
  };



struct spop_resize
  : public traits_op_default
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_resize>& in);
  };



struct spop_floor
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_floor>& in);
  };



struct spop_ceil
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_ceil>& in);
  };



struct spop_round
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_round>& in);
  };



struct spop_trunc
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_trunc>& in);
  };



struct spop_sign
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_sign>& in);
  };



struct spop_flipud
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_flipud>& in);
  };



struct spop_fliplr
  : public traits_op_passthru
  {
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_fliplr>& in);
  };



struct spop_replace
  : public traits_op_passthru
  {
  template<typename eT, typename T1>
  inline static void apply(SpMat<eT>& out, const mtSpOp<eT, T1, spop_replace>& in);
  };



//! @}
