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



//! \addtogroup spop_omit
//! @{



template<typename T1>
inline
void
spop_omit::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_omit>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword omit_mode = in.aux_uword_a;
  
  if(arma_config::fast_math_warn)
    {
    if(omit_mode == 1)  { arma_warn(1, "omit_nan(): detection of NaN is not reliable in fast math mode"); }
    if(omit_mode == 2)  { arma_warn(1, "omit_nonfinite(): detection of non-finite values is not reliable in fast math mode"); }
    }
  
  auto is_omitted_1 = [](const eT& x) -> bool { return arma_isnan(x);       };
  auto is_omitted_2 = [](const eT& x) -> bool { return arma_isnonfinite(x); };
  
  const SpProxy<T1> P(in.m);
  
  if(P.is_alias(out))
    {
    SpMat<eT> tmp;
    
    if(omit_mode == 1)  { spop_omit::apply_noalias(tmp, P, is_omitted_1); }
    if(omit_mode == 2)  { spop_omit::apply_noalias(tmp, P, is_omitted_2); }
    
    out.steal_mem(tmp);
    }
  else
    {
    if(omit_mode == 1)  { spop_omit::apply_noalias(out, P, is_omitted_1); }
    if(omit_mode == 2)  { spop_omit::apply_noalias(out, P, is_omitted_2); }
    }
  }



template<typename T1, typename functor>
inline
void
spop_omit::apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& P, functor is_omitted)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_rows        = P.get_n_rows();
  const uword max_n_nonzero = P.get_n_nonzero();
  
  if(max_n_nonzero == 0)  { out.reset(); return; }
  
  out.reserve(P.get_n_elem(), 1, max_n_nonzero);
  
  typename SpProxy<T1>::const_iterator_type it     = P.begin();
  typename SpProxy<T1>::const_iterator_type it_end = P.end();
  
  uword count = 0;
  
  for(; it != it_end; ++it)
    {
    const eT val = (*it);
    
    if(is_omitted(val) == false)
      {
      const uword index = it.row() + it.col()*n_rows;
      
      access::rw(out.values[count])      = val;
      access::rw(out.row_indices[count]) = index;
      access::rw(out.col_ptrs[1])++;
      ++count;
      }
    }
  
  if(count < max_n_nonzero)
    {
    if(count <= (max_n_nonzero/2))
      {
      out.mem_resize(count);
      }
    else
      {
      // quick resize without reallocating memory and copying data
      access::rw(         out.n_nonzero) = count;
      access::rw(     out.values[count]) = eT(0);
      access::rw(out.row_indices[count]) = uword(0);
      }
    }
  }



//! @}
