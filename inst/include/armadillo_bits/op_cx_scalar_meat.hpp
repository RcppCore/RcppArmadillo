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


//! \addtogroup op_cx_scalar
//! @{



template<typename T1>
inline
void
op_cx_scalar_times::apply
  (
        Mat< typename std::complex<typename T1::pod_type> >& out,
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] * k;
    }
  }



template<typename T1>
inline
void
op_cx_scalar_plus::apply
  (
        Mat< typename std::complex<typename T1::pod_type> >& out,
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_plus>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] + k;
    }
  }



template<typename T1>
inline
void
op_cx_scalar_minus_pre::apply
  (
        Mat< typename std::complex<typename T1::pod_type> >& out,
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_pre>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = k - A[i];
    }
  }



template<typename T1>
inline
void
op_cx_scalar_minus_post::apply
  (
        Mat< typename std::complex<typename T1::pod_type> >& out,
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_post>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] - k;
    }
  }



template<typename T1>
inline
void
op_cx_scalar_div_pre::apply
  (
        Mat< typename std::complex<typename T1::pod_type> >& out,
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_pre>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = k / A[i];
    }
  }



template<typename T1>
inline
void
op_cx_scalar_div_post::apply
  (
        Mat< typename std::complex<typename T1::pod_type> >& out,
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_post>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] / k;
    }
  }



//
//
//



template<typename T1>
inline
void
op_cx_scalar_times::apply
  (
           Cube< typename std::complex<typename T1::pod_type> >& out,
  const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] * k;
    }
  }



template<typename T1>
inline
void
op_cx_scalar_plus::apply
  (
           Cube< typename std::complex<typename T1::pod_type> >& out,
  const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_plus>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] + k;
    }
  }



template<typename T1>
inline
void
op_cx_scalar_minus_pre::apply
  (
           Cube< typename std::complex<typename T1::pod_type> >& out,
  const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_pre>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = k - A[i];
    }
  }



template<typename T1>
inline
void
op_cx_scalar_minus_post::apply
  (
           Cube< typename std::complex<typename T1::pod_type> >& out,
  const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_post>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] - k;
    }
  }



template<typename T1>
inline
void
op_cx_scalar_div_pre::apply
  (
           Cube< typename std::complex<typename T1::pod_type> >& out,
  const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_pre>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = k / A[i];
    }
  }



template<typename T1>
inline
void
op_cx_scalar_div_post::apply
  (
           Cube< typename std::complex<typename T1::pod_type> >& out,
  const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_post>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const eT  k       = X.aux_out_eT;
  const uword n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] / k;
    }
  }



//! @}
