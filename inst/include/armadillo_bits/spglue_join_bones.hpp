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


//! \addtogroup spglue_join
//! @{



class spglue_join_cols
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = false;
    static constexpr bool is_col  = (T1::is_col && T2::is_col);
    static constexpr bool is_xvec = false;
    };
  
  template<typename T1, typename T2>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_join_cols>& X);
  
  template<typename eT>
  inline static void apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B);
  
  template<typename eT, typename T1, typename T2, typename T3>
  inline static void apply(SpMat<eT>& out, const SpBase<eT,T1>& A, const SpBase<eT,T2>& B, const SpBase<eT,T3>& C);
  
  template<typename eT, typename T1, typename T2, typename T3, typename T4>
  inline static void apply(SpMat<eT>& out, const SpBase<eT,T1>& A, const SpBase<eT,T2>& B, const SpBase<eT,T3>& C, const SpBase<eT,T4>& D);
  };



class spglue_join_rows
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = (T1::is_row && T2::is_row);
    static constexpr bool is_col  = false;
    static constexpr bool is_xvec = false;
    };
  
  template<typename T1, typename T2>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_join_rows>& X);
  
  template<typename eT>
  inline static void apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B);
  
  template<typename eT, typename T1, typename T2, typename T3>
  inline static void apply(SpMat<eT>& out, const SpBase<eT,T1>& A, const SpBase<eT,T2>& B, const SpBase<eT,T3>& C);
  
  template<typename eT, typename T1, typename T2, typename T3, typename T4>
  inline static void apply(SpMat<eT>& out, const SpBase<eT,T1>& A, const SpBase<eT,T2>& B, const SpBase<eT,T3>& C, const SpBase<eT,T4>& D);
  };



//! @}
