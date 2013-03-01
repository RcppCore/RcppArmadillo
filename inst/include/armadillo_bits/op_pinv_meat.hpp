// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// Copyright (C) 2011 Stanislav Funiak
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_pinv
//! @{



template<typename eT>
inline
void
op_pinv::direct_pinv(Mat<eT>& out, const Mat<eT>& A, const eT in_tol)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  T tol = access::tmp_real(in_tol);
  
  arma_debug_check((tol < T(0)), "pinv(): tolerance must be >= 0");
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  // economical SVD decomposition 
  Mat<eT> U;
  Col< T> s;
  Mat<eT> V;
  
  const bool status = (n_cols > n_rows) ? auxlib::svd_econ(U,s,V,trans(A),'b') : auxlib::svd_econ(U,s,V,A,'b');
  
  if(status == false)
    {
    out.reset();
    arma_bad("pinv(): svd failed");
    return;
    }
  
  const uword s_n_elem = s.n_elem;
  const T*    s_mem    = s.memptr();
  
  // set tolerance to default if it hasn't been specified as an argument 
  if( (tol == T(0)) && (s_n_elem > 0) )
    {
    tol = (std::max)(n_rows, n_cols) * eop_aux::direct_eps( op_max::direct_max(s_mem, s_n_elem) );
    }
  
  
  // count non zero valued elements in s
  
  uword count = 0;
  
  for(uword i = 0; i < s_n_elem; ++i)
    {
    if(s_mem[i] > tol)
      {
      ++count;
      }
    }
  
  if(count != 0)
    {
    Col<T> s2(count);
    
    T* s2_mem = s2.memptr();
    
    uword count2 = 0;
    
    for(uword i=0; i < s_n_elem; ++i)
      {
      const T val = s_mem[i];
      
      if(val > tol)
        {
        s2_mem[count2] = T(1) / val;
        ++count2;
        }
      }
    
    
    if(n_rows >= n_cols)
      {
      out = ( V.n_cols > count ? V.cols(0,count-1) : V ) * diagmat(s2) * trans( U.n_cols > count ? U.cols(0,count-1) : U );
      }
    else
      {
      out = ( U.n_cols > count ? U.cols(0,count-1) : U ) * diagmat(s2) * trans( V.n_cols > count ? V.cols(0,count-1) : V );
      }
    }
  else
    {
    out.zeros(n_cols, n_rows);
    }
  }



template<typename T1>
inline
void
op_pinv::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_pinv>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  op_pinv::direct_pinv(out, A, in.aux);
  }



//! @}
