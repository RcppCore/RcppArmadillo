// Copyright (C) 2012 Ryan Curtin
// Copyright (C) 2012-2015 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup spop_sum
//! @{



template<typename T1>
arma_hot
inline
void
spop_sum::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_sum>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  arma_debug_check( (dim > 1), "sum(): parameter 'dim' must be 0 or 1" );
  
  const SpProxy<T1> p(in.m);
  
  if(p.is_alias(out) == false)
    {
    spop_sum::apply_noalias(out, p, dim);
    }
  else
    {
    SpMat<eT> tmp;
    
    spop_sum::apply_noalias(tmp, p, dim);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1>
arma_hot
inline
void
spop_sum::apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& p, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  const uword p_n_rows = p.get_n_rows();
  const uword p_n_cols = p.get_n_cols();
  
  if(dim == 0) // find the sum in each column
    {
    out.zeros(1, p_n_cols);
    
    if(p.get_n_nonzero() == 0)  { return; }
    
    if(SpProxy<T1>::must_use_iterator)
      {
      typename SpProxy<T1>::const_iterator_type it     = p.begin();
      typename SpProxy<T1>::const_iterator_type it_end = p.end();
      
      while(it != it_end)
        {
        out.at(0, it.col()) += (*it);
        ++it;
        }
      }
    else
      {
      for(uword col = 0; col < p_n_cols; ++col)
        {
        out.at(0,col) = arrayops::accumulate
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
    out.zeros(p_n_rows, 1);
    
    if(p.get_n_nonzero() == 0)  { return; }
    
    typename SpProxy<T1>::const_iterator_type it     = p.begin();
    typename SpProxy<T1>::const_iterator_type it_end = p.end();
    
    while(it != it_end)
      {
      out.at(it.row(), 0) += (*it);
      ++it;
      }
    }
  }



//! @}
