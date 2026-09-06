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



//! \addtogroup glue_cubemul
//! @{



struct glue_cubemul
  {
  template<typename T1, typename T2>
  inline static void apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_cubemul>& X);
  
  template<typename eT, bool do_strans_A, bool do_htrans_A, bool do_strans_B, bool do_htrans_B>
  inline static void apply_noalias(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B);
  
  //
  
  template<typename T1, typename T2>
  inline static Cube<typename T1::elem_type> apply(const BaseCube<typename T1::elem_type,T1>& expr_A, const Base<typename T1::elem_type,T2>& expr_B);
  
  template<typename eT, bool do_strans_A, bool do_htrans_A, bool do_strans_B, bool do_htrans_B>
  inline static void apply_noalias(Cube<eT>& out, const Cube<eT>& A, const Mat<eT>& B);
  
  //
  
  template<typename T1, typename T2>
  inline static Cube<typename T1::elem_type> apply(const Base<typename T1::elem_type,T1>& expr_A, const BaseCube<typename T1::elem_type,T2>& expr_B);
  
  template<typename eT, bool do_strans_A, bool do_htrans_A, bool do_strans_B, bool do_htrans_B>
  inline static void apply_noalias(Cube<eT>& out, const Mat<eT>& A, const Cube<eT>& B);
  };



//! @}
