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


//! \addtogroup op_inv_gen
//! @{



template<typename T1>
inline
void
op_inv_gen_default::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_gen_default>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_inv_gen_default::apply_direct(out, X.m, "inv()");
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("inv(): matrix is singular");
    }
  }



template<typename T1>
inline
bool
op_inv_gen_default::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, const char* caller_sig)
  {
  arma_extra_debug_sigprint();
  
  return op_inv_gen_full::apply_direct<T1,false>(out, expr, caller_sig, uword(0));
  }



//



template<typename T1>
inline
void
op_inv_gen_full::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_gen_full>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword flags = X.in_aux_uword_a;
  
  const bool status = op_inv_gen_full::apply_direct(out, X.m, "inv()", flags);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("inv(): matrix is singular");
    }
  }



template<typename T1, const bool has_user_flags>
inline
bool
op_inv_gen_full::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, const char* caller_sig, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(has_user_flags == true )  { arma_extra_debug_print("op_inv_gen_full: has_user_flags = true");  }
  if(has_user_flags == false)  { arma_extra_debug_print("op_inv_gen_full: has_user_flags = false"); }
  
  const bool tiny         = has_user_flags && bool(flags & inv_opts::flag_tiny        );
  const bool likely_sympd = has_user_flags && bool(flags & inv_opts::flag_likely_sympd);
  const bool no_sympd     = has_user_flags && bool(flags & inv_opts::flag_no_sympd    );
  
  arma_extra_debug_print("op_inv_gen_full: enabled flags:");
  
  if(tiny        )  { arma_extra_debug_print("tiny");         }
  if(likely_sympd)  { arma_extra_debug_print("likely_sympd"); }
  if(no_sympd    )  { arma_extra_debug_print("no_sympd");     }
  
  arma_debug_check( (no_sympd && likely_sympd), "inv(): options 'no_sympd' and 'likely_sympd' are mutually exclusive" );
  
  out = expr.get_ref();
  
  arma_debug_check( (out.is_square() == false), caller_sig, ": given matrix must be square sized" );
  
  if(tiny && (out.n_rows <= 4) && is_cx<eT>::no)
    {
    arma_extra_debug_print("op_inv_gen_full: attempting tinymatrix optimisation");
    
    const bool status = op_inv_gen_full::apply_tiny(out);
    
    if(status)  { return true; }
    
    arma_extra_debug_print("op_inv_gen_full: tinymatrix optimisation failed");
    
    // fallthrough if optimisation failed
    }
  
  if(is_op_diagmat<T1>::value || out.is_diagmat())
    {
    arma_extra_debug_print("op_inv_gen_full: detected diagonal matrix");
    
    const uword N = out.n_rows;
    
    eT* colmem = out.memptr();
    
    for(uword i=0; i<N; ++i)
      {
      eT& out_ii = colmem[i];
      
      const eT src_val = out_ii;
      const eT inv_val = eT(1) / src_val;
      
      if(src_val == eT(0))  { return false; }
      
      out_ii = inv_val;
      
      colmem += N;
      }
    
    return true;
    }
  
  const strip_trimat<T1> strip(expr.get_ref());
  
  const bool is_triu_expr = strip.do_triu;
  const bool is_tril_expr = strip.do_tril;
  
  const bool is_triu_mat = (is_triu_expr || is_tril_expr) ? false : (                        trimat_helper::is_triu(out));
  const bool is_tril_mat = (is_triu_expr || is_tril_expr) ? false : ((is_triu_mat) ? false : trimat_helper::is_tril(out));
  
  if(is_triu_expr || is_tril_expr || is_triu_mat || is_tril_mat)
    {
    return auxlib::inv_tr(out, ((is_triu_expr || is_triu_mat) ? uword(0) : uword(1)));
    }
  
  const bool try_sympd = arma_config::optimise_sympd && ((no_sympd) ? false : (likely_sympd ? true : sympd_helper::guess_sympd(out)));
  
  if(try_sympd)
    {
    arma_extra_debug_print("op_inv_gen_full: attempting sympd optimisation");
    
    Mat<eT> tmp = out;
    
    bool sympd_state = false;
    
    const bool status = auxlib::inv_sympd(tmp, sympd_state);
    
    if(status)  { out.steal_mem(tmp); return true; }
    
    if((status == false) && (sympd_state == true))  { return false; }
    
    arma_extra_debug_print("op_inv_gen_full: sympd optimisation failed");
    
    // fallthrough if optimisation failed
    }
  
  return auxlib::inv(out);
  }



template<typename eT>
arma_cold
inline
bool
op_inv_gen_full::apply_tiny(Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  // NOTE: assuming matrix X is square sized
  
  const uword N = X.n_rows;
  
  Mat<eT> out(N, N, arma_nozeros_indicator());
  
  constexpr T det_min =        std::numeric_limits<T>::epsilon();
  constexpr T det_max = T(1) / std::numeric_limits<T>::epsilon();
  
  const eT* Xm   =   X.memptr();
        eT* outm = out.memptr();
  
       if(N == 0)  { return true; }
  else if(N == 1)  { outm[0] = eT(1) / Xm[0]; }
  else if(N == 2)
    {
    const eT a = Xm[pos<0,0>::n2];
    const eT b = Xm[pos<0,1>::n2];
    const eT c = Xm[pos<1,0>::n2];
    const eT d = Xm[pos<1,1>::n2];
    
    const eT     det_val = (a*d - b*c);
    const  T abs_det_val = std::abs(det_val);
    
    if((abs_det_val < det_min) || (abs_det_val > det_max))  { return false; }
    
    outm[pos<0,0>::n2] =  d / det_val;
    outm[pos<0,1>::n2] = -b / det_val;
    outm[pos<1,0>::n2] = -c / det_val;
    outm[pos<1,1>::n2] =  a / det_val;
    }
  else if(N == 3)
    {
    const eT     det_val = op_det::apply_tiny(X);
    const  T abs_det_val = std::abs(det_val);
    
    if((abs_det_val < det_min) || (abs_det_val > det_max))  { return false; }
    
    outm[pos<0,0>::n3] =  (Xm[pos<2,2>::n3]*Xm[pos<1,1>::n3] - Xm[pos<2,1>::n3]*Xm[pos<1,2>::n3]) / det_val;
    outm[pos<1,0>::n3] = -(Xm[pos<2,2>::n3]*Xm[pos<1,0>::n3] - Xm[pos<2,0>::n3]*Xm[pos<1,2>::n3]) / det_val;
    outm[pos<2,0>::n3] =  (Xm[pos<2,1>::n3]*Xm[pos<1,0>::n3] - Xm[pos<2,0>::n3]*Xm[pos<1,1>::n3]) / det_val;
    
    outm[pos<0,1>::n3] = -(Xm[pos<2,2>::n3]*Xm[pos<0,1>::n3] - Xm[pos<2,1>::n3]*Xm[pos<0,2>::n3]) / det_val;
    outm[pos<1,1>::n3] =  (Xm[pos<2,2>::n3]*Xm[pos<0,0>::n3] - Xm[pos<2,0>::n3]*Xm[pos<0,2>::n3]) / det_val;
    outm[pos<2,1>::n3] = -(Xm[pos<2,1>::n3]*Xm[pos<0,0>::n3] - Xm[pos<2,0>::n3]*Xm[pos<0,1>::n3]) / det_val;
    
    outm[pos<0,2>::n3] =  (Xm[pos<1,2>::n3]*Xm[pos<0,1>::n3] - Xm[pos<1,1>::n3]*Xm[pos<0,2>::n3]) / det_val;
    outm[pos<1,2>::n3] = -(Xm[pos<1,2>::n3]*Xm[pos<0,0>::n3] - Xm[pos<1,0>::n3]*Xm[pos<0,2>::n3]) / det_val;
    outm[pos<2,2>::n3] =  (Xm[pos<1,1>::n3]*Xm[pos<0,0>::n3] - Xm[pos<1,0>::n3]*Xm[pos<0,1>::n3]) / det_val;
    
    const eT check_val = Xm[pos<0,0>::n3]*outm[pos<0,0>::n3] + Xm[pos<0,1>::n3]*outm[pos<1,0>::n3] + Xm[pos<0,2>::n3]*outm[pos<2,0>::n3];
    
    const  T max_diff  = (is_float<T>::value) ? T(1e-4) : T(1e-10);  // empirically determined; may need tuning
    
    if(std::abs(T(1) - check_val) >= max_diff)  { return false; }
    }
  else if(N == 4)
    {
    const eT     det_val = op_det::apply_tiny(X);
    const  T abs_det_val = std::abs(det_val);
    
    if((abs_det_val < det_min) || (abs_det_val > det_max))  { return false; }
    
    outm[pos<0,0>::n4] = ( Xm[pos<1,2>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<1,3>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,1>::n4] + Xm[pos<1,3>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,2>::n4] - Xm[pos<1,1>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,2>::n4] - Xm[pos<1,2>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,3>::n4] + Xm[pos<1,1>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<1,0>::n4] = ( Xm[pos<1,3>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,2>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,3>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,2>::n4] + Xm[pos<1,0>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,2>::n4] + Xm[pos<1,2>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,3>::n4] - Xm[pos<1,0>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<2,0>::n4] = ( Xm[pos<1,1>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,3>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,0>::n4] + Xm[pos<1,3>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,1>::n4] - Xm[pos<1,0>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<1,1>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,3>::n4] + Xm[pos<1,0>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<3,0>::n4] = ( Xm[pos<1,2>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,1>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,2>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,1>::n4] + Xm[pos<1,0>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,1>::n4] + Xm[pos<1,1>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,2>::n4] - Xm[pos<1,0>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,2>::n4] ) / det_val;
    
    outm[pos<0,1>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,2>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,3>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,1>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,2>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,3>::n4] - Xm[pos<0,1>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<1,1>::n4] = ( Xm[pos<0,2>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,0>::n4] + Xm[pos<0,3>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,0>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,2>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,3>::n4] + Xm[pos<0,0>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<2,1>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,1>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,0>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,1>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,3>::n4] - Xm[pos<0,0>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<3,1>::n4] = ( Xm[pos<0,1>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,2>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,0>::n4] + Xm[pos<0,2>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,0>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,1>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,0>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,2>::n4] ) / det_val;
    
    outm[pos<0,2>::n4] = ( Xm[pos<0,2>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,3>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,3>::n4] + Xm[pos<0,1>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<1,2>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,2>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,3>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<2,2>::n4] = ( Xm[pos<0,1>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,0>::n4] + Xm[pos<0,3>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,3>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,3>::n4] ) / det_val;
    outm[pos<3,2>::n4] = ( Xm[pos<0,2>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,1>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,2>::n4] ) / det_val;
    
    outm[pos<0,3>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,1>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,1>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,2>::n4] + Xm[pos<0,1>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,2>::n4] + Xm[pos<0,2>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,3>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,3>::n4] ) / det_val;
    outm[pos<1,3>::n4] = ( Xm[pos<0,2>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,0>::n4] + Xm[pos<0,3>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,2>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,2>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,3>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,3>::n4] ) / det_val;
    outm[pos<2,3>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,0>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,1>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,1>::n4] + Xm[pos<0,1>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,3>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,3>::n4] ) / det_val;
    outm[pos<3,3>::n4] = ( Xm[pos<0,1>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,0>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,0>::n4] + Xm[pos<0,2>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,1>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,1>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,2>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,2>::n4] ) / det_val;
    
    const eT check_val = Xm[pos<0,0>::n4]*outm[pos<0,0>::n4] + Xm[pos<0,1>::n4]*outm[pos<1,0>::n4] + Xm[pos<0,2>::n4]*outm[pos<2,0>::n4] + Xm[pos<0,3>::n4]*outm[pos<3,0>::n4];
    
    const  T max_diff  = (is_float<T>::value) ? T(1e-4) : T(1e-10);  // empirically determined; may need tuning
    
    if(std::abs(T(1) - check_val) >= max_diff)  { return false; }
    }
  else
    {
    return false;
    }
  
  arrayops::copy(X.memptr(), out.memptr(), out.n_elem);
  
  return true;
  }



template<typename T1>
inline
bool
op_inv_gen_rcond::apply_direct(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  out       = expr.get_ref();
  out_rcond = T(0);
  
  arma_debug_check( (out.is_square() == false), "inv(): given matrix must be square sized" );
  
  const uword N = out.n_rows;
  
  if(is_op_diagmat<T1>::value || out.is_diagmat())
    {
    arma_extra_debug_print("op_inv_gen_rcond: detected diagonal matrix");
    
    eT* colmem = out.memptr();
    
    T max_abs_src_val = T(0);
    T max_abs_inv_val = T(0);
    
    for(uword i=0; i<N; ++i)
      {
      eT& out_ii = colmem[i];
      
      const eT src_val = out_ii;
      const eT inv_val = eT(1) / src_val;
      
      if(src_val == eT(0))  { return false; }
      
      out_ii = inv_val;
      
      const T abs_src_val = std::abs(src_val);
      const T abs_inv_val = std::abs(inv_val);
      
      max_abs_src_val = (abs_src_val > max_abs_src_val) ? abs_src_val : max_abs_src_val;
      max_abs_inv_val = (abs_inv_val > max_abs_inv_val) ? abs_inv_val : max_abs_inv_val;
      
      colmem += N;
      }
    
    out_rcond = T(1) / (max_abs_src_val * max_abs_inv_val);
    
    return true;
    }
  
  const strip_trimat<T1> strip(expr.get_ref());
  
  const bool is_triu_expr = strip.do_triu;
  const bool is_tril_expr = strip.do_tril;
  
  const bool is_triu_mat = (is_triu_expr || is_tril_expr) ? false : (                        trimat_helper::is_triu(out));
  const bool is_tril_mat = (is_triu_expr || is_tril_expr) ? false : ((is_triu_mat) ? false : trimat_helper::is_tril(out));
  
  if(is_triu_expr || is_tril_expr || is_triu_mat || is_tril_mat)
    {
    return auxlib::inv_tr_rcond(out, out_rcond, ((is_triu_expr || is_triu_mat) ? uword(0) : uword(1)));
    }
  
  const bool try_sympd = arma_config::optimise_sympd && ((auxlib::crippled_lapack(out)) ? false : sympd_helper::guess_sympd(out));
  
  if(try_sympd)
    {
    arma_extra_debug_print("op_inv_gen_rcond: attempting sympd optimisation");
    
    Mat<eT> tmp = out;
    
    bool sympd_state = false;
    
    const bool status = auxlib::inv_sympd_rcond(tmp, sympd_state, out_rcond, T(-1));
    
    if(status)  { out.steal_mem(tmp); return true; }
    
    if((status == false) && (sympd_state == true))  { return false; }
    
    arma_extra_debug_print("op_inv_gen_rcond: sympd optimisation failed");
    
    // fallthrough if optimisation failed
    }
  
  return auxlib::inv_rcond(out, out_rcond);
  }



//! @}
