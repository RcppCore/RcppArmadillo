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



//! \addtogroup op_sp_nonzeros
//! @{



template<typename T1>
inline
void
op_sp_nonzeros::apply(Mat<typename T1::elem_type>& out, const SpToDOp<T1, op_sp_nonzeros>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_SpMat<T1>::value || is_SpMat<typename SpProxy<T1>::stored_type>::value)
    {
    const unwrap_spmat<T1> U(X.m);
    
    out.set_size(U.M.n_nonzero,1);
    
    arrayops::copy(out.memptr(), U.M.values, U.M.n_nonzero);
    
    return;
    }
  
  if(is_SpSubview<T1>::value)
    {
    const SpSubview<eT>& sv = reinterpret_cast< const SpSubview<eT>& >(X.m);
    
    if(sv.n_rows == sv.m.n_rows)
      {
      arma_debug_print("op_sp_nonzeros::apply(): SpSubview optimisation");
      
      const SpMat<eT>& m   = sv.m;
      const uword      col = sv.aux_col1;
      const uword      N   = sv.n_nonzero;
      
      out.set_size(N, 1);
      
      arrayops::copy(out.memptr(), &(m.values[ m.col_ptrs[col] ]), N);
      
      return;
      }
    }
  
  const SpProxy<T1> P(X.m);
  
  const uword N = P.get_n_nonzero();
  
  out.set_size(N,1);
  
  if(N == 0)  { return; }
  
  eT* out_mem = out.memptr();
  
  typename SpProxy<T1>::const_iterator_type it = P.begin();
  
  for(uword i=0; i<N; ++i)  { out_mem[i] = (*it); ++it; }
  }



//! @}
