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



//! \addtogroup op_shift
//! @{



class op_shift_vec
  : public traits_op_passthru
  {
  public:
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_shift_vec>& in);
  };



class op_shift
  : public traits_op_default
  {
  public:
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_shift>& in);
  
  template<typename eT> inline static void apply_direct(Mat<eT>& out, const Mat<eT>& X, const uword len, const uword neg, const uword dim);
  
  template<typename eT> inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword len, const uword neg, const uword dim);
  
  template<typename eT> inline static void apply_alias(Mat<eT>& out, const uword len, const uword neg, const uword dim);
  };



//! @}
