// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup op_range
//! @{



template<typename T1>
inline
void
op_range::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_range>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  arma_debug_check( (dim > 1), "range(): parameter 'dim' must be 0 or 1");
  
  const quasi_unwrap<T1> U(in.m);
  const Mat<eT>& X = U.M;
  
  if(U.is_alias(out) == false)
    {
    op_range::apply_noalias(out, X, dim);
    }
  else
    {
    Mat<eT> tmp;
    
    op_range::apply_noalias(tmp, X, dim);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT>
inline
void
op_range::apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  // TODO: replace with dedicated implementation which finds min and max at the same time
  out = max(X,dim) - min(X,dim);
  }



template<typename T1>
inline
typename T1::elem_type
op_range::vector_range(const T1& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(expr);
  const Mat<eT>& X = U.M;
  
  const eT*   X_mem = X.memptr();
  const uword N     = X.n_elem;
  
  if(N == 0)
    {
    arma_debug_check(true, "range(): object has no elements");
    
    return Datum<eT>::nan;
    }
  
  // TODO: replace with dedicated implementation which finds min and max at the same time
  return op_max::direct_max(X_mem, N) - op_min::direct_min(X_mem, N);
  }



//! @}
