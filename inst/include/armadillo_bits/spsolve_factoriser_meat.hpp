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


//! \addtogroup spsolve_factoriser
//! @{



template<typename worker_type>
inline
void
spsolve_factoriser::delete_worker()
  {
  arma_extra_debug_sigprint();
  
  if(worker_ptr != nullptr)
    {
    worker_type* ptr = reinterpret_cast<worker_type*>(worker_ptr);
    
    delete ptr;
    
    worker_ptr = nullptr;
    }
  }



inline
void
spsolve_factoriser::cleanup()
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_SUPERLU)
    {
         if(elem_type_indicator == 1)  { delete_worker< superlu_worker<    float> >(); }
    else if(elem_type_indicator == 2)  { delete_worker< superlu_worker<   double> >(); }
    else if(elem_type_indicator == 3)  { delete_worker< superlu_worker< cx_float> >(); }
    else if(elem_type_indicator == 4)  { delete_worker< superlu_worker<cx_double> >(); }
    }
  #endif
  
  worker_ptr          = nullptr;
  elem_type_indicator = 0;
  n_rows              = 0;
  rcond_value         = double(0);
  }



inline
spsolve_factoriser::~spsolve_factoriser()
  {
  arma_extra_debug_sigprint_this(this);
  
  cleanup();
  }



inline
spsolve_factoriser::spsolve_factoriser()
  {
  arma_extra_debug_sigprint_this(this);
  }



inline
void
spsolve_factoriser::reset()
  {
  arma_extra_debug_sigprint();
  
  cleanup();
  }



inline
double
spsolve_factoriser::rcond() const
  {
  arma_extra_debug_sigprint();
  
  return rcond_value;
  }



template<typename T1>
inline
bool
spsolve_factoriser::factorise
  (
  const SpBase<typename T1::elem_type,T1>& A_expr,
  const spsolve_opts_base&                 settings,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  #if defined(ARMA_USE_SUPERLU)
    {
    typedef typename T1::elem_type            eT;
    typedef typename get_pod_type<eT>::result  T;
    
    typedef superlu_worker<eT> worker_type;
    
    //
    
    cleanup();
    
    //
    
    const unwrap_spmat<T1> U(A_expr.get_ref());
    const SpMat<eT>& A =   U.M;
    
    if(A.is_square() == false)
      {
      arma_debug_warn_level(1, "spsolve_factoriser::factorise(): solving under-determined / over-determined systems is currently not supported");
      return false;
      }
    
    n_rows = A.n_rows;
    
    //
    
    superlu_opts superlu_opts_default;
    
    const superlu_opts& opts = (settings.id == 1) ? static_cast<const superlu_opts&>(settings) : superlu_opts_default;
    
    if( (opts.pivot_thresh < double(0)) || (opts.pivot_thresh > double(1)) )
      {
      arma_debug_warn_level(1, "spsolve_factoriser::factorise(): pivot_thresh must be in the [0,1] interval" );
      return false;
      }
    
    //
    
    worker_ptr = new(std::nothrow) worker_type;
    
    if(worker_ptr == nullptr)
      {
      arma_debug_warn_level(3, "spsolve_factoriser::factorise(): could not construct worker object");
      return false;
      }
    
    //
    
         if(    is_float<eT>::value)  { elem_type_indicator = 1; }
    else if(   is_double<eT>::value)  { elem_type_indicator = 2; }
    else if( is_cx_float<eT>::value)  { elem_type_indicator = 3; }
    else if(is_cx_double<eT>::value)  { elem_type_indicator = 4; }
    
    //
    
    worker_type* local_worker_ptr = reinterpret_cast<worker_type*>(worker_ptr);
    worker_type& local_worker_ref = (*local_worker_ptr);
    
    //
    
    T local_rcond_value = T(0);
    
    const bool status = local_worker_ref.factorise(local_rcond_value, A, opts);
    
    rcond_value = double(local_rcond_value);
    
    if( (status == false) || arma_isnan(local_rcond_value) || ((opts.allow_ugly == false) && (local_rcond_value < std::numeric_limits<T>::epsilon())) )
      {
      arma_debug_warn_level(3, "spsolve_factoriser::factorise(): factorisation failed; rcond: ", local_rcond_value);
      delete_worker<worker_type>();
      return false;
      }
    
    return true;
    }
  #else
    {
    arma_ignore(A_expr);
    arma_ignore(settings);
    arma_stop_logic_error("spsolve_factoriser::factorise(): use of SuperLU must be enabled");
    return false;
    }
  #endif
  }



template<typename T1>
inline
bool
spsolve_factoriser::solve
  (
         Mat<typename T1::elem_type>&    X,
  const Base<typename T1::elem_type,T1>& B_expr,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  #if defined(ARMA_USE_SUPERLU)
    {
    typedef typename T1::elem_type eT;
    
    typedef superlu_worker<eT> worker_type;
    
    if(worker_ptr == nullptr)
      {
      arma_debug_warn_level(2, "spsolve_factoriser::solve(): no factorisation available");
      X.soft_reset();
      return false;
      }
    
    bool type_mismatch = false;
    
         if(    (is_float<eT>::value) && (elem_type_indicator != 1) )  { type_mismatch = true; }
    else if(   (is_double<eT>::value) && (elem_type_indicator != 2) )  { type_mismatch = true; }
    else if( (is_cx_float<eT>::value) && (elem_type_indicator != 3) )  { type_mismatch = true; }
    else if((is_cx_double<eT>::value) && (elem_type_indicator != 4) )  { type_mismatch = true; }
    
    if(type_mismatch)
      {
      arma_debug_warn_level(1, "spsolve_factoriser::solve(): matrix type mismatch");
      X.soft_reset();
      return false;
      }
    
    const quasi_unwrap<T1> U(B_expr.get_ref());
    const Mat<eT>& B     = U.M;
    
    if(n_rows != B.n_rows)
      {
      arma_debug_warn_level(1, "spsolve_factoriser::solve(): matrix size mismatch");
      X.soft_reset();
      return false;
      }

    const bool is_alias = U.is_alias(X);
    
    Mat<eT>  tmp;
    Mat<eT>& out = is_alias ? tmp : X;
    
    worker_type* local_worker_ptr = reinterpret_cast<worker_type*>(worker_ptr);
    worker_type& local_worker_ref = (*local_worker_ptr);
    
    const bool status = local_worker_ref.solve(out,B);
    
    if(is_alias)  { X.steal_mem(tmp); }
    
    if(status == false)
      {
      arma_debug_warn_level(3, "spsolve_factoriser::solve(): solution not found");
      X.soft_reset();
      return false;
      }
    
    return true;
    }
  #else
    {
    arma_ignore(X);
    arma_ignore(B_expr);
    arma_stop_logic_error("spsolve_factoriser::solve(): use of SuperLU must be enabled");
    return false;
    }
  #endif
  }



//! @}
