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


//! \addtogroup op_sp_plus
//! @{


template<typename T1>
inline
void
op_sp_plus::apply(Mat<typename T1::elem_type>& out, const SpToDOp<T1,op_sp_plus>& in)
  {
  arma_debug_sigprint();
  
  const SpProxy<T1> P(in.m);
  
  out.set_size(P.get_n_rows(), P.get_n_cols());
  out.fill(in.aux);
  
  out += P.Q;
  }



// used for the optimization of sparse % (sparse + scalar)
template<typename eT, typename T2, typename T3>
inline
void
op_sp_plus::apply_inside_schur(SpMat<eT>& out, const T2& x, const SpToDOp<T3, op_sp_plus>& y)
  {
  arma_debug_sigprint();
  
  const SpProxy<T2> P2(x);
  const SpProxy<T3> P3(y.m);
  
  arma_conform_assert_same_size(P2.get_n_rows(), P2.get_n_cols(), P3.get_n_rows(), P3.get_n_cols(), "element-wise multiplication");
  
  out.zeros(P2.get_n_rows(), P2.get_n_cols());
  
  typename SpProxy<T2>::const_iterator_type it     = P2.begin();
  typename SpProxy<T2>::const_iterator_type it_end = P2.end();
  
  const eT k = y.aux;
  
  for(; it != it_end; ++it)
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    out.at(it_row, it_col) = (*it) * (P3.at(it_row, it_col) + k);
    }
  }



// used for the optimization of sparse / (sparse + scalar)
template<typename eT, typename T2, typename T3>
inline
void
op_sp_plus::apply_inside_div(SpMat<eT>& out, const T2& x, const SpToDOp<T3, op_sp_plus>& y)
  {
  arma_debug_sigprint();
  
  const SpProxy<T2> P2(x);
  const SpProxy<T3> P3(y.m);
  
  arma_conform_assert_same_size(P2.get_n_rows(), P2.get_n_cols(), P3.get_n_rows(), P3.get_n_cols(), "element-wise division");
  
  out.zeros(P2.get_n_rows(), P2.get_n_cols());
  
  typename SpProxy<T2>::const_iterator_type it     = P2.begin();
  typename SpProxy<T2>::const_iterator_type it_end = P2.end();
  
  const eT k = y.aux;
  
  for(; it != it_end; ++it)
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    out.at(it_row, it_col) = (*it) / (P3.at(it_row, it_col) + k);
    }
  }



//! @}
