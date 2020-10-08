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



//! \addtogroup op_find
//! @{



class op_find
  : public traits_op_col
  {
  public:
  
  template<typename T1>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const Base<typename T1::elem_type, T1>& X
    );
  
  template<typename T1, typename op_type>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const mtOp<uword, T1, op_type>& X,
    const typename arma_op_rel_only<op_type>::result* junk1 = nullptr,
    const typename arma_not_cx<typename T1::elem_type>::result* junk2 = nullptr
    );
  
  template<typename T1, typename op_type>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const mtOp<uword, T1, op_type>& X,
    const typename arma_op_rel_only<op_type>::result* junk1 = nullptr,
    const typename arma_cx_only<typename T1::elem_type>::result* junk2 = nullptr
    );
  
  template<typename T1, typename T2, typename glue_type>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const mtGlue<uword, T1, T2, glue_type>& X,
    const typename arma_glue_rel_only<glue_type>::result* junk1 = nullptr,
    const typename arma_not_cx<typename T1::elem_type>::result* junk2 = nullptr,
    const typename arma_not_cx<typename T2::elem_type>::result* junk3 = nullptr
    );
  
  template<typename T1, typename T2, typename glue_type>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const mtGlue<uword, T1, T2, glue_type>& X,
    const typename arma_glue_rel_only<glue_type>::result* junk1 = nullptr,
    const typename arma_cx_only<typename T1::elem_type>::result* junk2 = nullptr,
    const typename arma_cx_only<typename T2::elem_type>::result* junk3 = nullptr
    );
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const mtOp<uword, T1, op_find>& X);
  };



class op_find_simple
  : public traits_op_col
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const mtOp<uword, T1, op_find_simple>& X);
  };



class op_find_finite
  : public traits_op_col
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const mtOp<uword, T1, op_find_finite>& X);
  };



class op_find_nonfinite
  : public traits_op_col
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const mtOp<uword, T1, op_find_nonfinite>& X);
  };



//! @}
