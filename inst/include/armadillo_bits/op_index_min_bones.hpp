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


//! \addtogroup op_index_min
//! @{


class op_index_min
  : public traits_op_xvec
  {
  public:
  
  // dense matrices
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const mtOp<uword,T1,op_index_min>& in);
  
  template<typename eT>
  inline static void apply_noalias(Mat<uword>& out, const Mat<eT>& X, const uword dim);
  
  
  // cubes
  
  template<typename T1>
  inline static void apply(Cube<uword>& out, const mtOpCube<uword, T1, op_index_min>& in);  
  
  template<typename eT>
  inline static void apply_noalias(Cube<uword>& out, const Cube<eT>& X, const uword dim, const typename arma_not_cx<eT>::result* junk = nullptr);
  
  template<typename eT>
  inline static void apply_noalias(Cube<uword>& out, const Cube<eT>& X, const uword dim, const typename arma_cx_only<eT>::result* junk = nullptr);
  
  
  // sparse matrices
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const SpBase<typename T1::elem_type,T1>& expr, const uword dim);
  };



//! @}
