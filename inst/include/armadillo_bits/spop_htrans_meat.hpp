// Copyright (C) 2012 Ryan Curtin
// Copyright (C) 2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup spop_htrans
//! @{



template<typename T1>
arma_hot
inline
void
spop_htrans::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_htrans>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  spop_strans::apply_proxy(out, in.m);
  }



template<typename T1>
arma_hot
inline
void
spop_htrans::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_htrans>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const SpProxy<T1> p(in.m);
  
  if(p.is_alias(out) == false)
    {
    out.set_size( p.get_n_cols(), p.get_n_rows() );
    
    out.mem_resize(p.get_n_nonzero());
    
    typename SpProxy<T1>::const_row_iterator_type it = p.begin_row();
    
    while(it.pos() < p.get_n_nonzero())
      {
      access::rw(out.values[it.pos()]) = std::conj(*it);
      access::rw(out.row_indices[it.pos()]) = it.col(); // transpose
      ++access::rw(out.col_ptrs[it.row() + 1]);
      
      ++it;
      }
    
    // Fix column pointers.
    const uword out_n_cols = out.n_cols;
    
    for(uword c = 1; c <= out_n_cols; ++c)
      {
      access::rw(out.col_ptrs[c]) += out.col_ptrs[c - 1];
      }
    }
  else
    {
    SpMat<typename T1::elem_type> result( p.get_n_cols(), p.get_n_rows() );
    
    result.mem_resize(p.get_n_nonzero());
    
    typename SpProxy<T1>::const_row_iterator_type it = p.begin_row();
    
    while(it.pos() < p.get_n_nonzero())
      {
      access::rw(result.values[it.pos()]) = std::conj(*it);
      access::rw(result.row_indices[it.pos()]) = it.col(); // transpose
      ++access::rw(result.col_ptrs[it.row() + 1]);
      
      ++it;
      }
    
    // Fix column pointers.
    const uword result_n_cols = result.n_cols;
    
    for(uword c = 1; c <= result_n_cols; ++c)
      {
      access::rw(result.col_ptrs[c]) += result.col_ptrs[c - 1];
      }
    
    out.steal_mem(result);
    }
  }



//! @}
