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


//! \addtogroup spop_spones
//! @{



template<typename T1>
inline
void
spop_spones::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_spones>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(expr.m);
  
  if(U.is_alias(out))
    {
    arrayops::inplace_set( access::rwp(out.values), eT(1), out.n_nonzero );
    }
  else
    {
    SpMat<eT> tmp(arma_layout_indicator(), U.M);
    
    arrayops::inplace_set( access::rwp(tmp.values), eT(1), tmp.n_nonzero );
    
    out.steal_mem(tmp);
    }
  }



//! @}
