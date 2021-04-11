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


//! \addtogroup op_log_det
//! @{



template<typename T1>
inline
bool
op_log_det::apply_direct(typename T1::elem_type& out_val, typename T1::pod_type& out_sign, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(strip_diagmat<T1>::do_diagmat)
    {
    const strip_diagmat<T1> strip(expr.get_ref());
    
    return op_log_det::apply_diagmat(out_val, out_sign, strip.M);
    }
  
  if(strip_trimat<T1>::do_trimat)
    {
    const strip_trimat<T1> strip(expr.get_ref());
    
    return op_log_det::apply_trimat(out_val, out_sign, strip.M);
    }
  
  Mat<eT> A(expr.get_ref());
  
  arma_debug_check( (A.is_square() == false), "log_det(): given matrix must be square sized" );
  
  if(A.is_diagmat())  { return op_log_det::apply_diagmat(out_val, out_sign, A); }
  
  const bool is_triu =                   trimat_helper::is_triu(A);
  const bool is_tril = is_triu ? false : trimat_helper::is_tril(A);
  
  if(is_triu || is_tril)  { return op_log_det::apply_trimat(out_val, out_sign, A); }
  
  return auxlib::log_det(out_val, out_sign, A);
  }



template<typename T1>
inline
bool
op_log_det::apply_diagmat(typename T1::elem_type& out_val, typename T1::pod_type& out_sign, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const diagmat_proxy<T1> A(expr.get_ref());
  
  arma_debug_check( (A.n_rows != A.n_cols), "log_det(): given matrix must be square sized" );
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  if(N == 0)
    {
    out_val  = eT(0);
    out_sign =  T(1);
    
    return true;
    }
  
  eT x = A[0];
  
  T  sign = (is_cx<eT>::no) ?         ( (access::tmp_real(x) < T(0)) ?   T(-1) : T(1) ) : T(1);
  eT val  = (is_cx<eT>::no) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x    ) : std::log(x);
  
  for(uword i=1; i<N; ++i)
    {
    x = A[i];
    
    sign *= (is_cx<eT>::no) ?         ( (access::tmp_real(x) < T(0)) ?   T(-1) : T(1) ) : T(1);
    val  += (is_cx<eT>::no) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x    ) : std::log(x);
    }
  
  out_val  = val;
  out_sign = sign;
  
  return (arma_isnan(out_val) == false);
  }



template<typename T1>
inline
bool
op_log_det::apply_trimat(typename T1::elem_type& out_val, typename T1::pod_type& out_sign, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> P(expr.get_ref());
  
  const uword N = P.get_n_rows();
  
  arma_debug_check( (N != P.get_n_cols()), "log_det(): given matrix must be square sized" );
  
  if(N == 0)
    {
    out_val  = eT(0);
    out_sign =  T(1);
    
    return true;
    }
  
  eT x = P.at(0,0);
  
  T  sign = (is_cx<eT>::no) ?         ( (access::tmp_real(x) < T(0)) ?   T(-1) : T(1) ) : T(1);
  eT val  = (is_cx<eT>::no) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x    ) : std::log(x);
  
  for(uword i=1; i<N; ++i)
    {
    x = P.at(i,i);
    
    sign *= (is_cx<eT>::no) ?         ( (access::tmp_real(x) < T(0)) ?   T(-1) : T(1) ) : T(1);
    val  += (is_cx<eT>::no) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x    ) : std::log(x);
    }
  
  out_val  = val;
  out_sign = sign;
  
  return (arma_isnan(out_val) == false);
  }



//



template<typename T1>
inline
bool
op_log_det_sympd::apply_direct(typename T1::pod_type& out_val, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> A(expr.get_ref());
  
  arma_debug_check( (A.is_square() == false), "log_det_sympd(): given matrix must be square sized" );
  
  if((arma_config::debug) && (auxlib::rudimentary_sym_check(A) == false))
    {
    if(is_cx<eT>::no )  { arma_debug_warn_level(1, "log_det_sympd(): given matrix is not symmetric"); }
    if(is_cx<eT>::yes)  { arma_debug_warn_level(1, "log_det_sympd(): given matrix is not hermitian"); }
    }
  
  return auxlib::log_det_sympd(out_val, A);
  }



//! @}
