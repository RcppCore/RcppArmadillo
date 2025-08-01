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


//! \addtogroup op_omit
//! @{


class op_omit
  : public traits_op_col
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_omit>& in);
  
  template<typename T1, typename functor>
  inline static void apply(Mat<typename T1::elem_type>& out, const T1& X, functor is_omitted);
  };



class op_omit_cube
  : public traits_op_col
  {
  public:
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const CubeToMatOp<T1, op_omit_cube>& in);
  
  template<typename T1, typename functor>
  inline static void apply(Mat<typename T1::elem_type>& out, const T1& X, functor is_omitted);
  };



//! @}
