// Copyright (C) 2009-2013 Conrad Sanderson
// Copyright (C) 2009-2013 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Dimitrios Bouzas
// Copyright (C) 2011 Stanislav Funiak
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_pinv
//! @{



template<typename T1>
inline
void
op_pinv::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_pinv>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type            eT;
  typedef typename get_pod_type<eT>::result  T;
  
  const bool use_divide_and_conquer = (in.aux_uword_a == 1);
  
  T tol = access::tmp_real(in.aux);
  
  arma_debug_check((tol < T(0)), "pinv(): tolerance must be >= 0");
  
  const Proxy<T1> P(in.m);
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  if( (n_rows*n_cols) == 0 )
    {
    out.set_size(n_cols,n_rows);
    return;
    }
  
  
  // economical SVD decomposition 
  Mat<eT> U;
  Col< T> s;
  Mat<eT> V;
  
  bool status = false;
  
  if(use_divide_and_conquer)
    {
    status = (n_cols > n_rows) ? auxlib::svd_dc_econ(U, s, V, trans(P.Q)) : auxlib::svd_dc_econ(U, s, V, P.Q);
    }
  else
    {
    status = (n_cols > n_rows) ? auxlib::svd_econ(U, s, V, trans(P.Q), 'b') : auxlib::svd_econ(U, s, V, P.Q, 'b');
    }
  
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
    if(s_mem[i] > tol)  { ++count; }
    }
  
  
  if(count > 0)
    {
    Col<T> s2(count);
    
    T* s2_mem = s2.memptr();
    
    uword count2 = 0;
    
    for(uword i=0; i < s_n_elem; ++i)
      {
      const T val = s_mem[i];
      
      if(val > tol)  {  s2_mem[count2] = T(1) / val;  ++count2; }
      }
    
    
    if(n_rows >= n_cols)
      {
      out = ( (V.n_cols > count) ? V.cols(0,count-1) : V ) * diagmat(s2) * trans( (U.n_cols > count) ? U.cols(0,count-1) : U );
      }
    else
      {
      out = ( (U.n_cols > count) ? U.cols(0,count-1) : U ) * diagmat(s2) * trans( (V.n_cols > count) ? V.cols(0,count-1) : V );
      }
    }
  else
    {
    out.zeros(n_cols, n_rows);
    }
  }



//! @}
