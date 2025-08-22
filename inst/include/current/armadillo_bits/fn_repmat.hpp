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



//! \addtogroup fn_repmat
//! @{



template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_arma_type<T1>::value, const Op<T1, op_repmat> >::result
repmat(const T1& A, const uword r, const uword c)
  {
  arma_debug_sigprint();
  
  return Op<T1, op_repmat>(A, r, c);
  }



template<typename T1>
arma_warn_unused
inline
const SpOp<T1, spop_repmat>
repmat(const SpBase<typename T1::elem_type,T1>& A, const uword r, const uword c)
  {
  arma_debug_sigprint();
  
  return SpOp<T1, spop_repmat>(A.get_ref(), r, c);
  }



template<typename T1>
arma_warn_unused
inline
const OpCube<T1, op_repcube>
repcube(const BaseCube<typename T1::elem_type,T1>& A, const uword r, const uword c, const uword s)
  {
  arma_debug_sigprint();
  
  return OpCube<T1, op_repcube>(A.get_ref(), r, c, s);
  }



template<typename T1>
arma_warn_unused
inline
Cube<typename T1::elem_type>
repcube(const Base<typename T1::elem_type,T1>& A, const uword copies_per_row, const uword copies_per_col, const uword copies_per_slice)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT> B = arma::repmat(A.get_ref(), copies_per_row, copies_per_col);
  
  const eT*   B_mem    = B.memptr();
  const uword B_n_elem = B.n_elem;
  
  Cube<eT> out(B.n_rows, B.n_cols, copies_per_slice);
  
  if( (B_n_elem > 0) && (out.n_elem > 0) )
    {
    for(uword s=0; s < copies_per_slice; ++s)
      {
      arrayops::copy(out.slice_memptr(s), B_mem, B_n_elem);
      }
    }
  
  return out;
  }



//! @}
