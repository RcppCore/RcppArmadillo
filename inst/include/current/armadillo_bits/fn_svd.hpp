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


//! \addtogroup fn_svd
//! @{



template<typename T1>
inline
bool
svd
  (
         Col<typename T1::pod_type>&     S,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> A(X.get_ref());
  
  const bool status = auxlib::svd_dc(S, A);
  
  if(status == false)
    {
    S.soft_reset();
    arma_warn(3, "svd(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
arma_warn_unused
inline
Col<typename T1::pod_type>
svd
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  Col<T> out;
  
  Mat<eT> A(X.get_ref());
  
  const bool status = auxlib::svd_dc(out, A);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("svd(): decomposition failed");
    }
  
  return out;
  }



template<typename T1>
inline
bool
svd
  (
         Mat<typename T1::elem_type>&    U,
         Col<typename T1::pod_type >&    S,
         Mat<typename T1::elem_type>&    V,
  const Base<typename T1::elem_type,T1>& X,
  const char*                            method = "dc",
  const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  arma_conform_check
    (
    ( (void_ptr(&U) == void_ptr(&S)) || (&U == &V) || (void_ptr(&S) == void_ptr(&V)) ),
    "svd(): two or more output objects are the same object"
    );
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  arma_conform_check( ((sig != 's') && (sig != 'd')), "svd(): unknown method specified" );
  
  Mat<eT> A(X.get_ref());
  
  bool status = false;
  
  if(sig == 'd')
    {
    const uword N = (std::min)(A.n_rows, A.n_cols);
    
    const uword N_limit = (is_cx<eT>::yes) ? uword(20000) : uword(23000);
    
    const bool allow_dc = (sizeof(blas_int) >= std::size_t(8)) ? true : (N <= N_limit);
    
    if(allow_dc == false)
      {
      arma_warn(3, "svd(): matrix size too large for divide-and-conquer algorithm; using standard algorithm instead");
      }
    
    status = (allow_dc) ? auxlib::svd_dc(U, S, V, A) : auxlib::svd(U, S, V, A);
    }
  else
    {
    status = auxlib::svd(U, S, V, A);
    }
  
  if(status == false)
    {
    U.soft_reset();
    S.soft_reset();
    V.soft_reset();
    arma_warn(3, "svd(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
inline
bool
svd_econ
  (
         Mat<typename T1::elem_type>&    U,
         Col<typename T1::pod_type >&    S,
         Mat<typename T1::elem_type>&    V,
  const Base<typename T1::elem_type,T1>& X,
  const char                             mode,
  const char*                            method = "dc",
  const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  arma_conform_check
    (
    ( (void_ptr(&U) == void_ptr(&S)) || (&U == &V) || (void_ptr(&S) == void_ptr(&V)) ),
    "svd_econ(): two or more output objects are the same object"
    );
  
  arma_conform_check
    (
    ( (mode != 'l') && (mode != 'r') && (mode != 'b') ),
    "svd_econ(): parameter 'mode' is incorrect"
    );
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  arma_conform_check( ((sig != 's') && (sig != 'd')), "svd_econ(): unknown method specified" );
  
  Mat<eT> A(X.get_ref());
  
  bool status = false;
  
  if( (mode == 'b') && (sig == 'd') )
    {
    const uword N = (std::min)(A.n_rows, A.n_cols);
    
    const uword N_limit = (is_cx<eT>::yes) ? uword(20000) : uword(23000);
    
    const bool allow_dc = (sizeof(blas_int) >= std::size_t(8)) ? true : (N <= N_limit);
    
    if(allow_dc == false)
      {
      arma_warn(3, "svd_econ(): matrix size too large for divide-and-conquer algorithm; using standard algorithm instead");
      }
    
    status = (allow_dc) ? auxlib::svd_dc_econ(U, S, V, A) : auxlib::svd_econ(U, S, V, A, mode);
    }
  else
    {
    status = auxlib::svd_econ(U, S, V, A, mode);
    }
  
  if(status == false)
    {
    U.soft_reset();
    S.soft_reset();
    V.soft_reset();
    arma_warn(3, "svd_econ(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
inline
bool
svd_econ
  (
         Mat<typename T1::elem_type>&    U,
         Col<typename T1::pod_type >&    S,
         Mat<typename T1::elem_type>&    V,
  const Base<typename T1::elem_type,T1>& X,
  const char*                            mode   = "both",
  const char*                            method = "dc",
  const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  return svd_econ(U, S, V, X, ((mode != nullptr) ? mode[0] : char(0)), method);
  }



//! @}
