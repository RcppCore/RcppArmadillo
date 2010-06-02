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


//! \addtogroup op_relational
//! @{



template<typename T1>
inline
void
op_rel_lt_pre::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_lt_pre>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (val < A[i]) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_lt_post::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_lt_post>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] < val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_gt_pre::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_gt_pre>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (val > A[i]) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_gt_post::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_gt_post>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] > val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_lteq_pre::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_lteq_pre>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (val <= A[i]) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_lteq_post::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_lteq_post>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] <= val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_gteq_pre::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_gteq_pre>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (val >= A[i]) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_gteq_post::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_gteq_post>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] >= val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_eq::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_eq>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] == val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_noteq::apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_noteq>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.n_rows, A.n_cols);
  
  const u32  n_elem  = A.n_elem;
  const eT   val     = X.aux;
        u32* out_mem = out.memptr();
   
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A[i] != val) ? u32(1) : u32(0);
    }
  }



//! @}
