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


//! \addtogroup op_diagmat
//! @{



class op_diagmat
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Proxy<T1>& P);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_times>, op_diagmat>& X);
  
  template<typename T1, typename T2>
  inline static void apply_times(Mat<typename T1::elem_type>& out, const T1& X, const T2& Y, const typename arma_not_cx<typename T1::elem_type>::result* junk = nullptr);
  
  template<typename T1, typename T2>
  inline static void apply_times(Mat<typename T1::elem_type>& out, const T1& X, const T2& Y, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr);
  };



class op_diagmat2
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat2>& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword row_offset, const uword col_offset);
  };



//! @}
