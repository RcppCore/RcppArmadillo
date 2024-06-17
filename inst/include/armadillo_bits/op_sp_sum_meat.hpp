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


//! \addtogroup op_sp_sum
//! @{



template<typename T1>
inline
void
op_sp_sum::apply(Mat<typename T1::elem_type>& out, const mtSpReduceOp<typename T1::elem_type, T1, op_sp_sum>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  arma_conform_check( (dim > 1), "sum(): parameter 'dim' must be 0 or 1" );
  
  const SpProxy<T1> p(in.m);
  
  const uword p_n_rows = p.get_n_rows();
  const uword p_n_cols = p.get_n_cols();
  
  if(dim == 0)  { out.zeros(1, p_n_cols); }
  if(dim == 1)  { out.zeros(p_n_rows, 1); }
  
  if(p.get_n_nonzero() == 0)  { return; }
  
  eT* out_mem = out.memptr();
    
  if(dim == 0) // find the sum in each column
    {
    if(SpProxy<T1>::use_iterator)
      {
      typename SpProxy<T1>::const_iterator_type it = p.begin();
      
      const uword N = p.get_n_nonzero();
      
      for(uword i=0; i < N; ++i)
        {
        out_mem[it.col()] += (*it);
        ++it;
        }
      }
    else
      {
      for(uword col = 0; col < p_n_cols; ++col)
        {
        out_mem[col] = arrayops::accumulate
          (
          &p.get_values()[p.get_col_ptrs()[col]],
          p.get_col_ptrs()[col + 1] - p.get_col_ptrs()[col]
          );
        }
      }
    }
  else
  if(dim == 1)  // find the sum in each row
    {
    typename SpProxy<T1>::const_iterator_type it = p.begin();
    
    const uword N = p.get_n_nonzero();
    
    for(uword i=0; i < N; ++i)
      {
      out_mem[it.row()] += (*it);
      ++it;
      }
    }
  }



//! @}
