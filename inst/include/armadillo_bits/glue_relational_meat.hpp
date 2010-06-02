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


//! \addtogroup glue_relational
//! @{



template<typename T1, typename T2>
inline
void
glue_rel_lt::apply
  (
        Mat   <u32>& out,
  const mtGlue<u32, T1, T2, glue_rel_lt>& X
  )
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "operator<");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
        u32* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] < B[i]) ? u32(1) : u32(0);
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_rel_gt::apply
  (
        Mat   <u32>& out,
  const mtGlue<u32, T1, T2, glue_rel_gt>& X
  )
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "operator>");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
        u32* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] > B[i]) ? u32(1) : u32(0);
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_rel_lteq::apply
  (
        Mat   <u32>& out,
  const mtGlue<u32, T1, T2, glue_rel_lteq>& X
  )
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "operator<=");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
        u32* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] <= B[i]) ? u32(1) : u32(0);
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_rel_gteq::apply
  (
        Mat   <u32>& out,
  const mtGlue<u32, T1, T2, glue_rel_gteq>& X
  )
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "operator>=");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
        u32* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] >= B[i]) ? u32(1) : u32(0);
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_rel_eq::apply
  (
        Mat   <u32>& out,
  const mtGlue<u32, T1, T2, glue_rel_eq>& X
  )
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "operator==");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
        u32* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] == B[i]) ? u32(1) : u32(0);
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_rel_noteq::apply
  (
        Mat   <u32>& out,
  const mtGlue<u32, T1, T2, glue_rel_noteq>& X
  )
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "operator!=");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
        u32* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] != B[i]) ? u32(1) : u32(0);
    }
  
  }



//! @}
