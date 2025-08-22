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


//! \addtogroup fn_resize
//! @{



template<typename T1>
arma_warn_unused
inline
const Op<T1, op_resize>
resize(const Base<typename T1::elem_type,T1>& X, const uword in_n_rows, const uword in_n_cols)
  {
  arma_debug_sigprint();
  
  return Op<T1, op_resize>(X.get_ref(), in_n_rows, in_n_cols);
  }



template<typename T1>
arma_warn_unused
inline
const Op<T1, op_resize>
resize(const Base<typename T1::elem_type,T1>& X, const SizeMat& s)
  {
  arma_debug_sigprint();
  
  return Op<T1, op_resize>(X.get_ref(), s.n_rows, s.n_cols);
  }



template<typename T1>
arma_warn_unused
inline
const OpCube<T1, op_resize>
resize(const BaseCube<typename T1::elem_type,T1>& X, const uword in_n_rows, const uword in_n_cols, const uword in_n_slices)
  {
  arma_debug_sigprint();
  
  return OpCube<T1, op_resize>(X.get_ref(), in_n_rows, in_n_cols, in_n_slices);
  }



template<typename T1>
arma_warn_unused
inline
const OpCube<T1, op_resize>
resize(const BaseCube<typename T1::elem_type,T1>& X, const SizeCube& s)
  {
  arma_debug_sigprint();
  
  return OpCube<T1, op_resize>(X.get_ref(), s.n_rows, s.n_cols, s.n_slices);
  }



template<typename T1>
arma_warn_unused
inline
const SpOp<T1, spop_resize>
resize(const SpBase<typename T1::elem_type, T1>& X, const uword in_n_rows, const uword in_n_cols)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_resize>(X.get_ref(), in_n_rows, in_n_cols);
  }



template<typename T1>
arma_warn_unused
inline
const SpOp<T1, spop_resize>
resize(const SpBase<typename T1::elem_type, T1>& X, const SizeMat& s)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_resize>(X.get_ref(), s.n_rows, s.n_cols);
  }



template<typename oT>
arma_warn_unused
inline
field<oT>
resize(const field<oT>& A, const uword in_n_rows, const uword in_n_cols, const uword in_n_slices = uword(1))
  {
  arma_debug_sigprint();
  
  // better-than-nothing implementation
  
  field<oT> B(in_n_rows, in_n_cols, in_n_slices);
  
  if((B.n_elem > 0) && (A.n_elem > 0))
    {
    const uword end_row   = (std::min)(in_n_rows,   A.n_rows  ) - 1;
    const uword end_col   = (std::min)(in_n_cols,   A.n_cols  ) - 1;
    const uword end_slice = (std::min)(in_n_slices, A.n_slices) - 1;
    
    B.subfield(0, 0, 0, end_row, end_col, end_slice) = A.subfield(0, 0, 0, end_row, end_col, end_slice);
    }
  
  return B;
  }



template<typename oT>
arma_warn_unused
inline
field<oT>
resize(const field<oT>& A, const SizeMat& s)
  {
  arma_debug_sigprint();
  
  return resize(A, s.n_rows, s.n_cols);
  }



template<typename oT>
arma_warn_unused
inline
field<oT>
resize(const field<oT>& A, const SizeCube& s)
  {
  arma_debug_sigprint();
  
  return resize(A, s.n_rows, s.n_cols, s.n_slices);
  }



//! @}
