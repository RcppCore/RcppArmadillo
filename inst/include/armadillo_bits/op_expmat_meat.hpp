// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------



//! \addtogroup op_expmat
//! @{


//! implementation based on:
//! Cleve Moler, Charles Van Loan.
//! Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later.
//! SIAM Review, Vol. 45, No. 1, 2003, pp. 3-49.
//! http://dx.doi.org/10.1137/S00361445024180


template<typename T1>
inline
void
op_expmat::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_expmat>& expr)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_expmat::apply_direct(out, expr.m);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("expmat(): given matrix appears ill-conditioned");
    }
  }



template<typename T1>
inline
bool
op_expmat::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  if(is_op_diagmat<T1>::value)
    {
    out = expr.get_ref();  // force the evaluation of diagmat()
    
    arma_debug_check( (out.is_square() == false), "expmat(): given matrix must be square sized", [&](){ out.soft_reset(); } );
    
    const uword N = (std::min)(out.n_rows, out.n_cols);
    
    for(uword i=0; i<N; ++i)  { out.at(i,i) = std::exp( out.at(i,i) ); }
    
    return true;
    }
  
  Mat<eT> A = expr.get_ref();
  
  arma_debug_check( (A.is_square() == false), "expmat(): given matrix must be square sized" );
  
  if(A.is_diagmat())
    {
    arma_extra_debug_print("op_expmat: detected diagonal matrix");
    
    const uword N = (std::min)(A.n_rows, A.n_cols);
    
    out.zeros(N,N);
    
    for(uword i=0; i<N; ++i)  { out.at(i,i) = std::exp( A.at(i,i) ); }
    
    return true;
    }
  
  bool do_sym = false;
  
  if( (arma_config::optimise_sym) && (auxlib::crippled_lapack(A) == false) )
    {
    bool is_approx_sym   = false;
    bool is_approx_sympd = false;
    
    sym_helper::analyse_matrix(is_approx_sym, is_approx_sympd, A);
    
    do_sym = ((is_cx<eT>::no) ? (is_approx_sym) : (is_approx_sym && is_approx_sympd));
    }
  
  if(do_sym)
    {
    arma_extra_debug_print("op_expmat: symmetric/hermitian optimisation");
    
    Col< T> eigval;
    Mat<eT> eigvec;
    
    const bool eig_status = eig_sym_helper(eigval, eigvec, A, 'd', "expmat()");
    
    if(eig_status == false)  { return false; }
    
    eigval = exp(eigval);
    
    out = eigvec * diagmat(eigval) * eigvec.t();
    
    return true;
    }
  
  const T norm_val = arma::norm(A, "inf");
  
  if(arma_isfinite(norm_val) == false)  { return false; }
  
  const double log2_val = (norm_val > T(0)) ? double(eop_aux::log2(norm_val)) : double(0);
  
  int exponent = int(0);  std::frexp(log2_val, &exponent);
  
  const uword s = uword( (std::max)(int(0), exponent + int(1)) );
  
  A /= eT(eop_aux::pow(double(2), double(s)));
  
  T c = T(0.5);
  
  Mat<eT> E(A.n_rows, A.n_rows, fill::eye);  E += c * A;
  Mat<eT> D(A.n_rows, A.n_rows, fill::eye);  D -= c * A;
  
  Mat<eT> X = A;
  
  bool positive = true;
  
  const uword N = 6;
  
  for(uword i = 2; i <= N; ++i)
    {
    c = c * T(N - i + 1) / T(i * (2*N - i + 1));
    
    X = A * X;
    
    E += c * X;
    
    if(positive)  { D += c * X; }  else  { D -= c * X; }
    
    positive = (positive) ? false : true;
    }
  
  if( (D.internal_has_nonfinite()) || (E.internal_has_nonfinite()) )  { return false; }
  
  const bool status = solve(out, D, E, solve_opts::no_approx);
  
  if(status == false)  { return false; }
  
  for(uword i=0; i < s; ++i)  { out = out * out; }
  
  return true;
  }



template<typename T1>
inline
void
op_expmat_sym::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_expmat_sym>& in)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_expmat_sym::apply_direct(out, in.m);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("expmat_sym(): transformation failed");
    }
  }



template<typename T1>
inline
bool
op_expmat_sym::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::elem_type eT;
    typedef typename T1::pod_type   T;
    
    const unwrap<T1>   U(expr.get_ref());
    const Mat<eT>& X = U.M;
    
    arma_debug_check( (X.is_square() == false), "expmat_sym(): given matrix must be square sized" );
    
    if((arma_config::debug) && (arma_config::warn_level > 0) && (is_cx<eT>::yes) && (sym_helper::check_diag_imag(X) == false))
      {
      arma_debug_warn_level(1, "inv_sympd(): imaginary components on diagonal are non-zero");
      }
    
    if(is_op_diagmat<T1>::value || X.is_diagmat())
      {
      arma_extra_debug_print("op_expmat_sym: detected diagonal matrix");
      
      out = X;
      
      eT* colmem = out.memptr();
      
      const uword N = X.n_rows;
      
      for(uword i=0; i<N; ++i)
        {
        eT& out_ii      = colmem[i];
         T  out_ii_real = access::tmp_real(out_ii);
         
        out_ii = eT( std::exp(out_ii_real) );
        
        colmem += N;
        }
      
      return true;
      }
    
    Col< T> eigval;
    Mat<eT> eigvec;
    
    const bool status = eig_sym_helper(eigval, eigvec, X, 'd', "expmat_sym()");
    
    if(status == false)  { return false; }
    
    eigval = exp(eigval);
    
    out = eigvec * diagmat(eigval) * eigvec.t();
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(expr);
    arma_stop_logic_error("expmat_sym(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! @}
