// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_inv
//! @{


//! immediate inverse of a matrix, storing the result in a dense matrix
template<typename eT>
inline
void
op_inv::apply(Mat<eT>& out, const Mat<eT>& A, const bool slow)
  {
  arma_extra_debug_sigprint();
  
  // no need to check for aliasing, due to:
  // - auxlib::inv() copies A to out before inversion
  // - for 2x2 and 3x3 matrices the code is alias safe
  
  bool status = auxlib::inv(out, A, slow);
  
  if(status == false)
    {
    out.reset();
    arma_bad("inv(): matrix appears to be singular");
    }
  }



//! immediate inverse of T1, storing the result in a dense matrix
template<typename T1>
inline
void
op_inv::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv>& X)
  {
  arma_extra_debug_sigprint();
  
  const strip_diagmat<T1> strip(X.m);
  
  if(strip.do_diagmat == true)
    {
    op_inv::apply_diag(out, strip.M);
    }
  else
    {
    const uword mode = X.aux_uword_a;
    
    const bool status = (mode == 0) ? auxlib::inv(out, X.m) : auxlib::inv(out, X.m, true);
    
    if(status == false)
      {
      out.reset();
      arma_bad("inv(): matrix appears to be singular");
      }
    }
  }



template<typename T1>
inline
void
op_inv::apply_diag(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy_check<T1> A(X.get_ref(), out);
  
  const uword N = A.n_elem;
  
  out.set_size(N,N);
  
  for(uword col=0; col<N; ++col)
    {
    for(uword row=0; row<col; ++row)   { out.at(row,col) = eT(0); }
    
    out.at(col,col) = eT(1) / A[col];
    
    for(uword row=col+1; row<N; ++row) { out.at(row,col) = eT(0); }
    }
  
  }



//! inverse of T1 (triangular matrices)
template<typename T1>
inline
void
op_inv_tr::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_tr>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = auxlib::inv_tr(out, X.m, X.aux_uword_a);
  
  if(status == false)
    {
    out.reset();
    arma_bad("inv(): matrix appears to be singular");
    }
  }



//! inverse of T1 (symmetric positive definite matrices)
template<typename T1>
inline
void
op_inv_sympd::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_sympd>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = auxlib::inv_sympd(out, X.m, X.aux_uword_a);
  
  if(status == false)
    {
    out.reset();
    arma_bad("inv(): matrix appears to be singular");
    }
  }



//! @}
