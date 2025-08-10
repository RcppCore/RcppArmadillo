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


//! \addtogroup op_sp_stddev
//! @{



template<typename T1>
inline
void
op_sp_stddev::apply(Mat<typename T1::pod_type>& out, const mtSpReduceOp<typename T1::pod_type, T1, op_sp_stddev>& in)
  {
  arma_debug_sigprint();
  
  const uword norm_type = in.aux_uword_a;
  const uword dim       = in.aux_uword_b;
  
  arma_conform_check( (norm_type > 1), "stddev(): parameter 'norm_type' must be 0 or 1" );
  arma_conform_check( (dim > 1),       "stddev(): parameter 'dim' must be 0 or 1"       );
  
  const SpProxy<T1> p(in.m);
  
  const uword p_n_rows = p.get_n_rows();
  const uword p_n_cols = p.get_n_cols();
  
  if( (p_n_rows == 0) || (p_n_cols == 0) || (p.get_n_nonzero() == 0) )
    {
    if(dim == 0)  { out.zeros((p_n_rows > 0) ? 1 : 0, p_n_cols); }
    if(dim == 1)  { out.zeros(p_n_rows, (p_n_cols > 0) ? 1 : 0); }
    
    return;
    }
  
  op_sp_stddev::apply_slow(out, p, norm_type, dim);
  }



template<typename T1>
inline
void
op_sp_stddev::apply_slow
  (
        Mat<typename T1::pod_type>& out,
  const SpProxy<T1>&                p,
  const uword                       norm_type,
  const uword                       dim
  )
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type  in_eT;
  //typedef typename T1::pod_type  out_eT;
  
  const uword p_n_rows = p.get_n_rows();
  const uword p_n_cols = p.get_n_cols();
  
  if(dim == 0)  // find variance in each column
    {
    arma_debug_print("op_sp_stddev::apply_slow(): dim = 0");
    
    out.zeros(1, p_n_cols);
    
    for(uword col = 0; col < p_n_cols; ++col)
      {
      if(SpProxy<T1>::use_iterator)
        {
        // We must use an iterator; we can't access memory directly.
        typename SpProxy<T1>::const_iterator_type it  = p.begin_col(col);
        typename SpProxy<T1>::const_iterator_type end = p.begin_col(col + 1);
        
        const uword n_zero = p_n_rows - (end.pos() - it.pos());
        
        // in_eT is used just to get the specialization right (complex / noncomplex)
        out.at(0, col) = std::sqrt( op_sp_var::iterator_var(it, end, n_zero, norm_type, in_eT(0)) );
        }
      else
        {
        // We can use direct memory access to calculate the variance.
        out.at(0, col) = std::sqrt(
          op_sp_var::direct_var
            (
            &p.get_values()[p.get_col_ptrs()[col]],
            p.get_col_ptrs()[col + 1] - p.get_col_ptrs()[col],
            p_n_rows,
            norm_type
            )
          );
        }
      }
    }
  else
  if(dim == 1)  // find variance in each row
    {
    arma_debug_print("op_sp_stddev::apply_slow(): dim = 1");
    
    out.zeros(p_n_rows, 1);
    
    for(uword row = 0; row < p_n_rows; ++row)
      {
      // We have to use an iterator here regardless of whether or not we can
      // directly access memory.
      typename SpProxy<T1>::const_row_iterator_type it  = p.begin_row(row);
      typename SpProxy<T1>::const_row_iterator_type end = p.end_row(row);
      
      const uword n_zero = p_n_cols - (end.pos() - it.pos());
      
      out.at(row, 0) = std::sqrt( op_sp_var::iterator_var(it, end, n_zero, norm_type, in_eT(0)) );
      }
    }
  }



template<typename T1>
inline
typename T1::pod_type
op_sp_stddev::stddev_vec
  (
  const T1& X,
  const uword norm_type
  )
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  arma_conform_check( (norm_type > 1), "stddev(): parameter 'norm_type' must be 0 or 1" );
  
  // conditionally unwrap it into a temporary and then directly operate.
  
  const unwrap_spmat<T1> tmp(X);
  
  if(tmp.M.n_elem == 0)
    {
    arma_conform_check(true, "stddev(): object has no elements");
    
    return Datum<T>::nan;
    }
  
  return std::sqrt( op_sp_var::direct_var(tmp.M.values, tmp.M.n_nonzero, tmp.M.n_elem, norm_type) );
  }



//! @}
