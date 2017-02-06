// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup glue_polyval
//! @{



template<typename eT>
inline
void
glue_polyval::apply_noalias(Mat<eT>& out, const Mat<eT>& P, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(X.n_rows, X.n_cols);
  
  const eT*   P_mem    = P.memptr();
  const uword P_n_elem = P.n_elem;
  
  out.fill(P_mem[0]);

  for(uword i=1; i < P_n_elem; ++i)
    {
    out = out % X + P_mem[i];
    }
  }



template<typename T1, typename T2>
inline
void
glue_polyval::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_polyval>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UP(expr.A);
  const quasi_unwrap<T2> UX(expr.B);
  
  const Mat<eT>& P = UP.M;
  const Mat<eT>& X = UX.M;
  
  arma_debug_check( ((P.is_vec() == false) && (P.is_empty() == false)), "polyval(): argument P must be a vector" );
  
  if(P.is_empty() || X.is_empty())
    {
    out.zeros(X.n_rows, X.n_cols);
    return;
    }
  
  if(UP.is_alias(out) || UX.is_alias(out))
    {
    Mat<eT> tmp;
    glue_polyval::apply_noalias(tmp, P, X);
    out.steal_mem(tmp);
    }
  else
    {
    glue_polyval::apply_noalias(out, P, X);
    }
  }



//! @}
