// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
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

  arma_debug_check( (A.get_n_rows() != A.get_n_cols()), "trace(): matrix must be square sized" );
  
  const uword N   = A.get_n_rows();
        eT  val = eT(0);
  
  for(uword i=0; i<N; ++i)
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
  
  const uword N = A.n_elem;
  
  eT val = eT(0);
  
  for(uword i=0; i<N; ++i)
    {
    val += A[i];
    }
  
  return val;
  }


//! speedup for trace(A*B), where the result of A*B is a square sized matrix
template<typename T1, typename T2>
inline
arma_warn_unused
typename T1::elem_type
trace(const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A, B, "matrix multiply");
  
  arma_debug_check( (A.n_rows != B.n_cols), "trace(): matrix must be square sized" );
  
  const uword N1  = A.n_rows;
  const uword N2  = A.n_cols;
        eT  val = eT(0);
  
  for(uword i=0; i<N1; ++i)
    {
    const eT* B_colmem = B.colptr(i);
          eT  acc      = eT(0);
    
    for(uword j=0; j<N2; ++j)
      {
      acc += A.at(i,j) * B_colmem[j];
      }
    
    val += acc;
    }
  
  return val;
  }



//! @}
