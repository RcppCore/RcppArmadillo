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


//! \addtogroup op_sp_diagvec
//! @{



template<typename T1>
inline
void
op_sp_diagvec::apply(Mat<typename T1::elem_type>& out, const mtSpReduceOp<typename T1::elem_type, T1,op_sp_diagvec>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(in.m);
  const SpMat<eT>& X =   U.M;
  
  const uword a = in.aux_uword_a;
  const uword b = in.aux_uword_b;
  
  const uword row_offset = (b >  0) ? a : 0;
  const uword col_offset = (b == 0) ? a : 0;
  
  arma_conform_check_bounds
    (
    ((row_offset > 0) && (row_offset >= X.n_rows)) || ((col_offset > 0) && (col_offset >= X.n_cols)),
    "diagvec(): requested diagonal out of bounds"
    );
  
  const uword len = (std::min)(X.n_rows - row_offset, X.n_cols - col_offset);
  
  out.set_size(len, 1);
  
  eT* out_mem = out.memptr();
  
  for(uword i=0; i < len; ++i)
    {
    out_mem[i] = X.at(i + row_offset, i + col_offset);
    }
  }



//! @}
