// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_trace
//! @{


//! Immediate trace (sum of diagonal elements) of a square dense matrix
template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
trace(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.get_ref());

  arma_debug_check( (A.n_rows != A.n_cols), "trace(): matrix must be square" );
  
  eT val = eT(0);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    val += A.at(i,i);
    }
  
  return val;
  }



template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
trace(const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy<T1> A(X.m);
  
  const u32 N = A.n_elem;
  
  eT val = eT(0);
  
  for(u32 i=0; i<N; ++i)
    {
    val += A[i];
    }
  
  return val;
  }



//! @}
