// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup solve_factoriser
//! @{



template<typename eT>
inline
bool
solve_factoriser_trimat_helper<eT>::train(Mat<eT>& AA, typename get_pod_type<eT>::result& out_rcond, const Mat<eT>& A, const uword layout)
  {
  arma_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename get_pod_type<eT>::result T;
    
    out_rcond = T(0);
    
    AA = A;
    
    if(AA.is_empty())  { return true; }
    
    out_rcond = auxlib::rcond_trimat(AA, layout);
    
    return true;
    }
  #else
    {
    arma_ignore(AA);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(layout);
    arma_stop_logic_error("use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
solve_factoriser_trimat_helper<eT>::solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA, const uword layout)
  {
  arma_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(AA.is_empty() || B.is_empty())  { out.zeros(AA.n_rows, B.n_cols); return true; }
    
    out = B;
    
    arma_conform_assert_blas_size(out);
    
    char     uplo  = (layout == 0) ? 'U' : 'L';
    char     trans = 'N';
    char     diag  = 'N';
    blas_int n     = blas_int(AA.n_rows);
    blas_int nrhs  = blas_int( B.n_cols);
    blas_int info  = blas_int(0);
    
    arma_debug_print("lapack::trtrs()");
    lapack::trtrs(&uplo, &trans, &diag, &n, &nrhs, AA.memptr(), &n, out.memptr(), &n, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(B);
    arma_ignore(AA);
    arma_ignore(layout);
    arma_stop_logic_error("use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//



template<typename eT>
inline
bool
solve_factoriser_sympd_helper<eT>::train(Mat<eT>& AA, typename get_pod_type<eT>::result& out_rcond, const Mat<typename get_pod_type<eT>::result>& A)
  {
  arma_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename get_pod_type<eT>::result T;
    
    out_rcond = T(0);
    
    AA = A;
    
    if(AA.is_empty())  { return true; }
    
    arma_conform_assert_blas_size(AA);
    
    char     norm_id  = '1';
    char     uplo     = 'L';
    blas_int n        = blas_int(AA.n_rows);  // assuming AA is square
    blas_int info     = blas_int(0);
    T        norm_val = T(0);
    
    podarray<T> work(AA.n_rows);
    
    arma_debug_print("lapack::lansy()");
    norm_val = (has_blas_float_bug<eT>::value) ? auxlib::norm1_sym(AA) : lapack::lansy(&norm_id, &uplo, &n, AA.memptr(), &n, work.memptr());
    
    arma_debug_print("lapack::potrf()");
    lapack::potrf<eT>(&uplo, &n, AA.memptr(), &n, &info);
    
    if(info != 0)  { AA.reset(); return false; }
    
    out_rcond = auxlib::lu_rcond_sympd<T>(AA, norm_val);
    
    return true;
    }
  #else
    {
    arma_ignore(AA);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_stop_logic_error("use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
solve_factoriser_sympd_helper<eT>::train(Mat<eT>& AA, typename get_pod_type<eT>::result& out_rcond, const Mat< std::complex<typename get_pod_type<eT>::result> >& A)
  {
  arma_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename get_pod_type<eT>::result T;
    
    out_rcond = T(0);
    
    AA = A;
    
    if(AA.is_empty())  { return true; }
    
    arma_conform_assert_blas_size(AA);
    
    char     norm_id  = '1';
    char     uplo     = 'L';
    blas_int n        = blas_int(AA.n_rows);  // assuming AA is square
    blas_int info     = blas_int(0);
    T        norm_val = T(0);
    
    podarray<T> work(AA.n_rows);
    
    arma_debug_print("lapack::lanhe()");
    norm_val = (has_blas_float_bug<eT>::value) ? auxlib::norm1_sym(AA) : lapack::lanhe(&norm_id, &uplo, &n, AA.memptr(), &n, work.memptr());
    
    arma_debug_print("lapack::potrf()");
    lapack::potrf<eT>(&uplo, &n, AA.memptr(), &n, &info);
    
    if(info != 0)  { AA.reset(); return false; }
    
    out_rcond = auxlib::lu_rcond_sympd<T>(AA, norm_val);
    
    return true;
    }
  #else
    {
    arma_ignore(AA);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_stop_logic_error("use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
solve_factoriser_sympd_helper<eT>::solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA)
  {
  arma_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(AA.is_empty() || B.is_empty())  { out.zeros(AA.n_rows, B.n_cols); return true; }
    
    out = B;
    
    arma_conform_assert_blas_size(out);
    
    char     uplo = 'L';
    blas_int n    = blas_int(AA.n_rows);  // assuming AA is square
    blas_int nrhs = blas_int( B.n_cols);
    blas_int info = blas_int(0);
    
    arma_debug_print("lapack::potrs()");
    lapack::potrs<eT>(&uplo, &n, &nrhs, AA.memptr(), &n, out.memptr(), &n, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(B);
    arma_ignore(AA);
    arma_stop_logic_error("use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//



template<typename eT>
inline
bool
solve_factoriser_square_helper<eT>::train(Mat<eT>& AA, podarray<blas_int>& ipiv, typename get_pod_type<eT>::result& out_rcond, const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename get_pod_type<eT>::result T;
    
    out_rcond = T(0);
    
    AA = A;
    
    if(AA.is_empty())  { return true; }
    
    arma_conform_assert_blas_size(AA);
    
    char     norm_id  = '1';
    blas_int n        = blas_int(AA.n_rows);  // assuming AA is square
    blas_int lda      = blas_int(AA.n_rows);
    blas_int info     = blas_int(0);
    T        norm_val = T(0);
    
    podarray<T> junk(1);
    
    ipiv.set_size(AA.n_rows + 2);  // +2 for paranoia
    
    arma_debug_print("lapack::lange()");
    norm_val = (has_blas_float_bug<eT>::value) ? auxlib::norm1_gen(A) : lapack::lange<eT>(&norm_id, &n, &n, AA.memptr(), &lda, junk.memptr());
    
    arma_debug_print("lapack::getrf()");
    lapack::getrf<eT>(&n, &n, AA.memptr(), &n, ipiv.memptr(), &info);
    
    if(info != blas_int(0))  { AA.reset(); ipiv.reset(); return false; }
    
    out_rcond = auxlib::lu_rcond<T>(AA, norm_val);
    
    return true;
    }
  #else
    {
    arma_ignore(AA);
    arma_ignore(ipiv);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_stop_logic_error("use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
solve_factoriser_square_helper<eT>::solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA, const podarray<blas_int>& ipiv)
  {
  arma_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(AA.is_empty() || B.is_empty())  { out.zeros(AA.n_rows, B.n_cols); return true; }
    
    out = B;
    
    arma_conform_assert_blas_size(out);
    
    char     trans = 'N';
    blas_int n     = blas_int(AA.n_rows);  // assuming AA is square
    blas_int lda   = blas_int(AA.n_rows);
    blas_int ldb   = blas_int( B.n_rows);
    blas_int nrhs  = blas_int( B.n_cols);
    blas_int info  = blas_int(0);
    
    arma_debug_print("lapack::getrs()");
    lapack::getrs<eT>(&trans, &n, &nrhs, AA.memptr(), &lda, ipiv.memptr(), out.memptr(), &ldb, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(B);
    arma_ignore(AA);
    arma_ignore(ipiv);
    arma_stop_logic_error("use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//



template<typename eT>
inline
bool
solve_factoriser_rect_helper<eT>::train(Mat<eT>& AA, typename get_pod_type<eT>::result& out_rcond, const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  Mat<eT> A_copy = A;
  
  // obtain the pseudo-inverse of A, under the assumption of full rank
  
  const bool status = auxlib::solve_rect_rcond(AA, out_rcond, A_copy, eye<Mat<eT>>(A.n_rows, A.n_rows));  // A_copy is overwritten
  
  if(status == false)  { AA.reset(); return false; }
  
  return true;
  }



template<typename eT>
inline
bool
solve_factoriser_rect_helper<eT>::solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA)
  {
  arma_debug_sigprint();
  
  if(AA.n_cols != B.n_rows)  { return false; }
  
  out = AA * B;
  
  return true;
  }



//



template<typename eT>
inline
bool
solve_factoriser_approx_helper<eT>::train(Mat<eT>& AA, const Mat<eT>& A)
  {
  arma_debug_sigprint();
  
  Mat<eT> A_copy = A;
  
  // obtain the pseudo-inverse of A; not using pinv() to ensure the processing is the same as standard arma::solve()
  
  const bool status = auxlib::solve_approx_svd(AA, A_copy, eye<Mat<eT>>(A.n_rows, A.n_rows));  // A_copy is overwritten
  
  if(status == false)  { AA.reset(); return false; }
  
  return true;
  }



template<typename eT>
inline
bool
solve_factoriser_approx_helper<eT>::solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA)
  {
  arma_debug_sigprint();
  
  if(AA.n_cols != B.n_rows)  { return false; }
  
  out = AA * B;
  
  return true;
  }




//
//
//



inline
solve_factoriser::~solve_factoriser()
  {
  arma_debug_sigprint();
  
  reset();
  }



inline
solve_factoriser::solve_factoriser()
  {
  arma_debug_sigprint_this(this);
  }



template<typename eT>
inline
bool
solve_factoriser::acquire_AA_ptr()
  {
  arma_debug_sigprint();
  
  AA_ptr_raw = new(std::nothrow) Mat<eT>;
  
  return (AA_ptr_raw != nullptr);
  }



template<typename eT>
inline
void
solve_factoriser::release_AA_ptr()
  {
  arma_debug_sigprint();
  
  if(AA_ptr_raw != nullptr)  { delete static_cast< Mat<eT>* >(AA_ptr_raw); AA_ptr_raw = nullptr; }
  }



inline
void
solve_factoriser::reset()
  {
  arma_debug_sigprint();
  
       if(elem_id == elem_id_enum::elem_float    )  { (*this).release_AA_ptr<    float>(); }
  else if(elem_id == elem_id_enum::elem_double   )  { (*this).release_AA_ptr<   double>(); }
  else if(elem_id == elem_id_enum::elem_cx_float )  { (*this).release_AA_ptr< cx_float>(); }
  else if(elem_id == elem_id_enum::elem_cx_double)  { (*this).release_AA_ptr<cx_double>(); }
  
  stored_rcond_val = double(0);
  layout           = 0;
  n_rows           = 0;
  
  elem_id = elem_id_enum::elem_none;
  fact_id = fact_id_enum::fact_none;
  
  ipiv.reset();
  }



template<typename eT>
inline
void
solve_factoriser::set_elem_id()
  {
  arma_debug_sigprint();
  
       if(    is_float<eT>::value)  { elem_id = elem_id_enum::elem_float;     }
  else if(   is_double<eT>::value)  { elem_id = elem_id_enum::elem_double;    }
  else if( is_cx_float<eT>::value)  { elem_id = elem_id_enum::elem_cx_float;  }
  else if(is_cx_double<eT>::value)  { elem_id = elem_id_enum::elem_cx_double; }
  else                              { elem_id = elem_id_enum::elem_none;      }
  }



template<typename eT>
inline
bool
solve_factoriser::check_elem_id() const
  {
  arma_debug_sigprint();
  
  bool ok = false;
  
       if(    (is_float<eT>::value) && (elem_id == elem_id_enum::elem_float    ))  { ok = true; }
  else if(   (is_double<eT>::value) && (elem_id == elem_id_enum::elem_double   ))  { ok = true; }
  else if( (is_cx_float<eT>::value) && (elem_id == elem_id_enum::elem_cx_float ))  { ok = true; }
  else if((is_cx_double<eT>::value) && (elem_id == elem_id_enum::elem_cx_double))  { ok = true; }
  
  return ok;
  }



inline
double
solve_factoriser::rcond() const
  {
  arma_debug_sigprint();
  
  return stored_rcond_val;
  }



template<typename T1>
inline
bool
solve_factoriser::factorise
  (
  const Base<typename T1::elem_type,T1>& A_expr,
  const solve_opts::opts&                opts,
  const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk
  )
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type            eT;
  typedef typename get_pod_type<eT>::result  T;
  
  if(arma_config::lapack == false)  { arma_stop_logic_error("solve_factoriser::factorise(): use of LAPACK must be enabled"); return false; }
  
  reset();
  
  set_elem_id<eT>();
  
  const bool AA_ptr_ok = acquire_AA_ptr<eT>();
  
  if(AA_ptr_ok == false)  { arma_warn(3, "solve_factoriser::factorise(): initialisation failed" ); return false; }
  
  Mat<eT>* AA_ptr_cooked = static_cast< Mat<eT>* >(AA_ptr_raw);
  Mat<eT>& AA            = (*AA_ptr_cooked);
  
  const uword flags = opts.flags;
  
  const bool fast         = bool(flags & solve_opts::flag_fast        );
  const bool equilibrate  = bool(flags & solve_opts::flag_equilibrate );
  const bool no_approx    = bool(flags & solve_opts::flag_no_approx   );
  const bool no_band      = bool(flags & solve_opts::flag_no_band     );
  const bool no_sympd     = bool(flags & solve_opts::flag_no_sympd    );
  const bool allow_ugly   = bool(flags & solve_opts::flag_allow_ugly  );
  const bool likely_sympd = bool(flags & solve_opts::flag_likely_sympd);
  const bool refine       = bool(flags & solve_opts::flag_refine      );
  const bool no_trimat    = bool(flags & solve_opts::flag_no_trimat   );
  const bool force_approx = bool(flags & solve_opts::flag_force_approx);
  const bool force_sym    = bool(flags & solve_opts::flag_force_sym   );
  const bool scale_thresh = bool(flags & solve_opts::flag_scale_thresh);
  
  arma_debug_print("solve_factoriser::factorise(): enabled flags:");
  
  if(fast        )  { arma_debug_print("fast");         }
  if(equilibrate )  { arma_debug_print("equilibrate");  }
  if(no_approx   )  { arma_debug_print("no_approx");    }
  if(no_band     )  { arma_debug_print("no_band");      }
  if(no_sympd    )  { arma_debug_print("no_sympd");     }
  if(allow_ugly  )  { arma_debug_print("allow_ugly");   }
  if(likely_sympd)  { arma_debug_print("likely_sympd"); }
  if(refine      )  { arma_debug_print("refine");       }
  if(no_trimat   )  { arma_debug_print("no_trimat");    }
  if(force_approx)  { arma_debug_print("force_approx"); }
  if(force_sym   )  { arma_debug_print("force_sym");    }
  if(scale_thresh)  { arma_debug_print("scale_thresh"); }
  
  arma_conform_check( fast,                        "solve_factoriser::factorise(): option 'fast' is not supported"                                );
  arma_conform_check( equilibrate,                 "solve_factoriser::factorise(): option 'equilibrate' is not supported"                         );
  arma_conform_check( refine,                      "solve_factoriser::factorise(): option 'refine' is not supported"                              );
  arma_conform_check( force_sym,                   "solve_factoriser::factorise(): option 'force_sym' is not supported"                           );
  arma_conform_check( (no_sympd  && likely_sympd), "solve_factoriser::factorise(): options 'no_sympd' and 'likely_sympd' are mutually exclusive"  );
  arma_conform_check( (no_approx && force_approx), "solve_factoriser::factorise(): options 'no_approx' and 'force_approx' are mutually exclusive" );
  
  const bool is_triu_expr = is_trimatu_expr<T1>::eval(A_expr.get_ref());
  const bool is_tril_expr = is_trimatl_expr<T1>::eval(A_expr.get_ref());
  
  arma_debug_print("is_triu_expr: ", is_triu_expr);
  arma_debug_print("is_tril_expr: ", is_tril_expr);
  
  const plain_unwrap<T1> U(A_expr.get_ref());
  const Mat<eT>& A     = U.M;
  
  if(arma_config::check_nonfinite && A.internal_has_nonfinite())  { arma_warn(3, "solve_factoriser::factorise(): given matrix has non-finite elements" ); return false; }
  
  n_rows = A.n_rows;
  
  if(A.is_empty())  { AA.copy_size(A); return true; }
  
  T    rcond_val = T(0);
  bool status    = false;
  
  if(force_approx)
    {
    arma_debug_print("solve_factoriser::factorise(): forced approximate factorisation");
    
    if(likely_sympd)  { arma_warn(2, "solve_factoriser::factorise(): option 'likely_sympd' ignored for forced approximate factorisation" ); }
    if(scale_thresh)  { arma_warn(2, "solve_factoriser::factorise(): option 'scale_thresh' ignored for forced approximate factorisation" ); }
    
    status = solve_factoriser_approx_helper<eT>::train(AA, A);
    
    if(status)  { fact_id = fact_id_enum::fact_approx; }
    
    return status;
    }
  
  const T rcond_threshold = (scale_thresh) ? T(std::numeric_limits<T>::epsilon() * (std::max)(uword(1), (std::max)(A.n_rows, A.n_cols))) : T(std::numeric_limits<T>::epsilon());
  
  arma_debug_print("solve_factoriser::factorise(): rcond_threshold: ", rcond_threshold);
  
  if(A.is_square())
    {
    arma_debug_print("solve_factoriser::factorise(): detected square system");
    
    const bool is_triu = (no_trimat || likely_sympd           ) ? false : ((is_triu_expr) ? true : trimat_helper::is_triu(A));
    const bool is_tril = (no_trimat || likely_sympd || is_triu) ? false : ((is_tril_expr) ? true : trimat_helper::is_tril(A));
    
    const bool try_sympd = arma_config::optimise_sym && ( (no_sympd || is_triu || is_tril) ? false : (likely_sympd ? true : sym_helper::guess_sympd(A, uword(16))) );
    
    arma_debug_print("solve_factoriser::factorise(): internal flags:");
    arma_debug_print("is_triu:   ", is_triu  );
    arma_debug_print("is_tril:   ", is_tril  );
    arma_debug_print("try_sympd: ", try_sympd);
    
    if(is_triu || is_tril)
      {
      if(is_triu)  { arma_debug_print("solve_factoriser::factorise(): upper triangular"); }
      if(is_tril)  { arma_debug_print("solve_factoriser::factorise(): lower triangular"); }
      
      layout = (is_triu) ? uword(0) : uword(1);
      
      status = solve_factoriser_trimat_helper<eT>::train(AA, rcond_val, A, layout);
      
      if(status)  { fact_id = fact_id_enum::fact_trimat; }
      }
    else
    if(try_sympd)
      {
      arma_debug_print("solve_factoriser::factorise(): sympd");
      
      status = solve_factoriser_sympd_helper<eT>::train(AA, rcond_val, A);
      
      if(status)
        {
        fact_id = fact_id_enum::fact_sympd;
        }
      else
        {
        arma_debug_print("solve_factoriser::factorise(): sympd failed; retrying with standard square solver");
        
        status = solve_factoriser_square_helper<eT>::train(AA, ipiv, rcond_val, A);
        
        if(status)  { fact_id = fact_id_enum::fact_square; }
        }
      }
    else
      {
      arma_debug_print("solve_factoriser::factorise(): square");
      
      status = solve_factoriser_square_helper<eT>::train(AA, ipiv, rcond_val, A);
      
      if(status)  { fact_id = fact_id_enum::fact_square; }
      }
    }
  else
    {
    arma_debug_print("solve_factoriser::factorise(): detected non-square system");
    
    if(likely_sympd)  { arma_warn(2,  "solve_factoriser::factorise(): option 'likely_sympd' ignored for non-square matrix" ); }
    
    status = solve_factoriser_rect_helper<eT>::train(AA, rcond_val, A);
    
    if(status)  { fact_id = fact_id_enum::fact_rect; }
    }
  
  
  if( (status == true) && (allow_ugly == false) && ((rcond_val < rcond_threshold) || arma_isnan(rcond_val)) )
    {
    status = false;
    }
  
  
  stored_rcond_val = double(rcond_val);
  
  
  if( (status == false) && (no_approx == false) )
    {
    arma_debug_print("solve_factoriser::factorise(): approx");
    
    if(rcond_val == T(0))
      {
      arma_warn(2, "solve_factoriser::factorise(): matrix is singular; attempting approx factorisation");
      }
    else
      {
      arma_warn(2, "solve_factoriser::factorise(): matrix is close to singular; rcond: ", rcond_val, "; attempting approx factorisation");
      }
    
    status = solve_factoriser_approx_helper<eT>::train(AA, A);
    
    if(status)  { fact_id = fact_id_enum::fact_approx; }
    }
  
  
  if(status == false)
    {
      AA.reset();
    ipiv.reset();
    
    fact_id = fact_id_enum::fact_none;
    
    arma_warn(3, "solve_factoriser::factorise(): factorisation failed" );
    }
  
  return status;
  }



template<typename T1>
inline
bool
solve_factoriser::solve
  (
         Mat<typename T1::elem_type>&    X,
  const Base<typename T1::elem_type,T1>& B_expr,
  const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk
  )
  const
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  if(arma_config::lapack == false)  { arma_stop_logic_error("solve_factoriser::solve(): use of LAPACK must be enabled"); return false; }
  if(check_elem_id<eT>() == false)  { arma_stop_logic_error("solve_factoriser::solve(): element type mismatch");         return false; }
  
  if(fact_id == fact_id_enum::fact_none)  { arma_warn(2, "solve_factoriser::solve(): no factorisation available"); X.soft_reset(); return false; }
  
  const Mat<eT>* AA_ptr_cooked = static_cast< const Mat<eT>* >(AA_ptr_raw);
  const Mat<eT>& AA            = (*AA_ptr_cooked);
  
  const quasi_unwrap<T1> U(B_expr.get_ref());
  const Mat<eT>& B     = U.M;
  
  if(n_rows != B.n_rows)  { arma_warn(1, "solve_factoriser::solve(): size mismatch between factorisation and given matrix" ); X.soft_reset(); return false; }

  if(AA.is_empty() || B.is_empty())  { X.zeros(AA.n_rows, B.n_cols); return true; }
  
  bool status = false;
  
  if(U.is_alias(X))
    {
    Mat<eT> tmp;
    
         if(fact_id == fact_id_enum::fact_trimat)  { status = solve_factoriser_trimat_helper<eT>::solve(tmp, B, AA, layout); }
    else if(fact_id == fact_id_enum::fact_sympd )  { status =  solve_factoriser_sympd_helper<eT>::solve(tmp, B, AA        ); }
    else if(fact_id == fact_id_enum::fact_square)  { status = solve_factoriser_square_helper<eT>::solve(tmp, B, AA, ipiv  ); }
    else if(fact_id == fact_id_enum::fact_rect  )  { status =   solve_factoriser_rect_helper<eT>::solve(tmp, B, AA        ); }
    else if(fact_id == fact_id_enum::fact_approx)  { status = solve_factoriser_approx_helper<eT>::solve(tmp, B, AA        ); }
    
    X.steal_mem(tmp);
    }
  else
    {
         if(fact_id == fact_id_enum::fact_trimat)  { status = solve_factoriser_trimat_helper<eT>::solve(X, B, AA, layout); }
    else if(fact_id == fact_id_enum::fact_sympd )  { status =  solve_factoriser_sympd_helper<eT>::solve(X, B, AA        ); }
    else if(fact_id == fact_id_enum::fact_square)  { status = solve_factoriser_square_helper<eT>::solve(X, B, AA, ipiv  ); }
    else if(fact_id == fact_id_enum::fact_rect  )  { status =   solve_factoriser_rect_helper<eT>::solve(X, B, AA        ); }
    else if(fact_id == fact_id_enum::fact_approx)  { status = solve_factoriser_approx_helper<eT>::solve(X, B, AA        ); }
    }
  
  if(status == false)  { arma_warn(3, "solve_factoriser::solve(): solution not found"); X.soft_reset(); }
  
  return status;
  }



//! @}
