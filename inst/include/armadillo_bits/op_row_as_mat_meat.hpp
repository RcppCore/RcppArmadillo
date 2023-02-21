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


//! \addtogroup op_row_as_mat
//! @{



template<typename T1>
inline
void
op_row_as_mat::apply(Mat<typename T1::elem_type>& out, const CubeToMatOp<T1, op_row_as_mat>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> U(expr.m);
  const Cube<eT>&   A = U.M;
  
  const uword in_row = expr.aux_uword;
  
  arma_debug_check_bounds( (in_row >= A.n_rows), "Cube::row_as_mat(): index out of bounds" );
  
  const uword A_n_cols   = A.n_cols;
  const uword A_n_rows   = A.n_rows;
  const uword A_n_slices = A.n_slices;
  
  out.set_size(A_n_slices, A_n_cols);
  
  for(uword s=0; s < A_n_slices; ++s)
    {
    const eT*   A_mem = &(A.at(in_row, 0, s));
          eT* out_mem = &(out.at(s,0));
    
    for(uword c=0; c < A_n_cols; ++c)
      {
      (*out_mem) = (*A_mem);
      
        A_mem += A_n_rows;
      out_mem += A_n_slices;
      }
    }
  }



//! @}
