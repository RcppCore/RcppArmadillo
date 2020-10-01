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


//! \addtogroup glue_affmul
//! @{



class glue_affmul
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = T1::is_row;
    static constexpr bool is_col  = T2::is_col;
    static constexpr bool is_xvec = false;
    };
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_affmul>& X);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(Mat<typename T1::elem_type>& out, const T1& A, const T2& B);
  
  template<typename T1, typename T2>
  inline static void apply_noalias_square(Mat<typename T1::elem_type>& out, const T1& A, const T2& B);
  
  template<typename T1, typename T2>
  inline static void apply_noalias_rectangle(Mat<typename T1::elem_type>& out, const T1& A, const T2& B);
  
  template<typename T1, typename T2>
  inline static void apply_noalias_generic(Mat<typename T1::elem_type>& out, const T1& A, const T2& B);
  };



//! @}

