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


//! \addtogroup spop_relational
//! @{



class spop_rel_lt_pre
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_lt_pre>& X);
  };



class spop_rel_lt_post
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_lt_post>& X);
  };



class spop_rel_gt_pre
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_gt_pre>& X);
  };



class spop_rel_gt_post
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_gt_post>& X);
  };



class spop_rel_lteq_pre
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_lteq_pre>& X);
  };



class spop_rel_lteq_post
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_lteq_post>& X);
  };



class spop_rel_gteq_pre
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_gteq_pre>& X);
  };



class spop_rel_gteq_post
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_gteq_post>& X);
  };



class spop_rel_eq
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_eq>& X);
  };



class spop_rel_noteq
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_noteq>& X);
  };



//! @}
