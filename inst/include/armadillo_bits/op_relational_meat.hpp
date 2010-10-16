// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (val < A[i]) ? u32(1) : u32(0);
    out_mem[j] = (val < A[j]) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] < val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] < val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (val > A[i]) ? u32(1) : u32(0);
    out_mem[j] = (val > A[j]) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] > val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] > val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (val <= A[i]) ? u32(1) : u32(0);
    out_mem[j] = (val <= A[j]) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] <= val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] <= val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (val >= A[i]) ? u32(1) : u32(0);
    out_mem[j] = (val >= A[j]) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] >= val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] >= val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] == val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] == val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const Proxy<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] != val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] != val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (A[i] != val) ? u32(1) : u32(0);
    }
  }



//
//
//



template<typename T1>
inline
void
op_rel_lt_pre::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_lt_pre>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (val < A[i]) ? u32(1) : u32(0);
    out_mem[j] = (val < A[j]) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (val < A[i]) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_lt_post::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_lt_post>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] < val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] < val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (A[i] < val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_gt_pre::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_gt_pre>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (val > A[i]) ? u32(1) : u32(0);
    out_mem[j] = (val > A[j]) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (val > A[i]) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_gt_post::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_gt_post>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] > val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] > val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (A[i] > val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_lteq_pre::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_lteq_pre>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (val <= A[i]) ? u32(1) : u32(0);
    out_mem[j] = (val <= A[j]) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (val <= A[i]) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_lteq_post::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_lteq_post>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] <= val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] <= val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (A[i] <= val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_gteq_pre::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_gteq_pre>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (val >= A[i]) ? u32(1) : u32(0);
    out_mem[j] = (val >= A[j]) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (val >= A[i]) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_gteq_post::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_gteq_post>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] >= val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] >= val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (A[i] >= val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_eq::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_eq>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] == val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] == val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (A[i] == val) ? u32(1) : u32(0);
    }
  }



template<typename T1>
inline
void
op_rel_noteq::apply(Cube<u32>& out, const mtOpCube<u32, T1, op_rel_noteq>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type;
  
  const ProxyCube<T1> Y(X.m);
  
  out.set_size(Y.get_n_rows(), Y.get_n_cols(), Y.get_n_slices());
  
  const eT      val     = X.aux;
        ea_type A       = Y.get_ea();
        u32*    out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    out_mem[i] = (A[i] != val) ? u32(1) : u32(0);
    out_mem[j] = (A[j] != val) ? u32(1) : u32(0);
    }
  
  if(i < n_elem)
    {
    out_mem[i] = (A[i] != val) ? u32(1) : u32(0);
    }
  }



//! @}
