// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_misc
//! @{



template<typename T1>
inline
void
op_real::apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_real>& X )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const uword n_elem  = out.n_elem;
        T*  out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::real(A[i]);
    }
  }



template<typename T1>
inline
void
op_real::apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_real>& X )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const uword n_elem  = out.n_elem;
        T*  out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::real(A[i]);
    }
  }



template<typename T1>
inline
void
op_imag::apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_imag>& X )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const uword n_elem  = out.n_elem;
        T*  out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::imag(A[i]);
    }
  }



template<typename T1>
inline
void
op_imag::apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_imag>& X )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const uword n_elem  = out.n_elem;
        T*  out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::imag(A[i]);
    }
  }



template<typename T1>
inline
void
op_abs::apply( Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_abs>& X )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  const uword n_elem  = out.n_elem;
        T*  out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::abs(A[i]);
    }
  }



template<typename T1>
inline
void
op_abs::apply( Cube<typename T1::pod_type>& out, const mtOpCube<typename T1::pod_type, T1, op_abs>& X )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const ProxyCube<T1> A(X.m);
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
  const uword n_elem  = out.n_elem;
        T*  out_mem = out.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::abs(A[i]);
    }
  }



template<typename T1>
inline
void
op_sympd::apply( Mat<typename T1::elem_type>& out, const Op<T1, op_sympd>& X )
  {
  arma_extra_debug_sigprint();
  
  out = X.m;
  }



//! @}
