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


//! \addtogroup op_norm2est
//! @{



template<typename eT>
inline
norm2est_randu_filler<eT>::norm2est_randu_filler()
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::mt19937_64::result_type local_seed_type;
  
  local_engine.seed(local_seed_type(123));
  
  typedef typename std::uniform_real_distribution<eT>::param_type local_param_type;
  
  local_u_distr.param(local_param_type(-1.0, +1.0));
  }


template<typename eT>
inline
void
norm2est_randu_filler<eT>::fill(eT* mem, const uword N)
  {
  arma_extra_debug_sigprint();
  
  for(uword i=0; i<N; ++i)  { mem[i] = eT( local_u_distr(local_engine) ); }
  }


//


template<typename T>
inline
norm2est_randu_filler< std::complex<T> >::norm2est_randu_filler()
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::mt19937_64::result_type local_seed_type;
  
  local_engine.seed(local_seed_type(123));
  
  typedef typename std::uniform_real_distribution<T>::param_type local_param_type;
  
  local_u_distr.param(local_param_type(-1.0, +1.0));
  }


template<typename T>
inline
void
norm2est_randu_filler< std::complex<T> >::fill(std::complex<T>* mem, const uword N)
  {
  arma_extra_debug_sigprint();
  
  for(uword i=0; i<N; ++i)
    {
    std::complex<T>& mem_i = mem[i];
    
    mem_i.real( T(local_u_distr(local_engine)) );
    mem_i.imag( T(local_u_distr(local_engine)) );
    }
  }



//
//
//



template<typename T1>
inline
typename T1::pod_type
op_norm2est::norm2est
  (
  const Base<typename T1::elem_type, T1>& X,
  const typename T1::pod_type             tolerance,
  const uword                             max_iter
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type   T;
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (tolerance <      T(0)), "norm2est(): parameter 'tolerance' must be > 0" );
  arma_debug_check( (max_iter  == uword(0)), "norm2est(): parameter 'max_iter' must be > 0"  );
  
  const T tol = (tolerance == T(0)) ? T(1e-6) : T(tolerance);
  
  const quasi_unwrap<T1> U(X.get_ref());
  const Mat<eT>& A     = U.M;
  
  if(A.n_elem == 0)  { return T(0); }
  
  if(A.internal_has_nonfinite())  { arma_debug_warn_level(1, "norm2est(): given matrix has non-finite elements"); }
  
  if((A.n_rows == 1) || (A.n_cols == 1))  { return op_norm::vec_norm_2( Proxy< Mat<eT> >(A) ); }
  
  norm2est_randu_filler<eT> randu_filler;
  
  Col<eT> x(A.n_rows, fill::none);
  Col<eT> y(A.n_cols, fill::none);
  
  randu_filler.fill(y.memptr(), y.n_elem);
  
  T est_old = 0;
  T est_cur = 0;
  
  for(uword i=0; i<max_iter; ++i)
    {
    arma_extra_debug_print(arma_str::format("norm2est(): iteration: %u") % i);
    
    x = A * y;
    
    T x_norm = op_norm::vec_norm_2( Proxy< Col<eT> >(x) );
    
    if(x_norm == T(0) || (arma_isfinite(x_norm) == false) || (x.internal_has_nonfinite()))
      {
      randu_filler.fill(x.memptr(), x.n_elem);
      
      x_norm = op_norm::vec_norm_2( Proxy< Col<eT> >(x) );
      }
    
    if(x_norm != T(0))  { x /= x_norm; }
    
    y = A.t() * x;
    
    est_old = est_cur;
    est_cur = op_norm::vec_norm_2( Proxy< Col<eT> >(y) );
    
    arma_extra_debug_print(arma_str::format("norm2est(): est_old: %e") % est_old);
    arma_extra_debug_print(arma_str::format("norm2est(): est_cur: %e") % est_cur);
    
    if(arma_isfinite(est_cur) == false)  { return est_old; }
    
    if( ((std::abs)(est_cur - est_old)) <= (tol * (std::max)(est_cur,est_old)) )  { break; }
    }
  
  return est_cur;
  }



//
//
//



template<typename T1>
inline
typename T1::pod_type
op_norm2est::norm2est
  (
  const SpBase<typename T1::elem_type, T1>& X,
  const typename T1::pod_type               tolerance,
  const uword                               max_iter
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type   T;
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (tolerance <      T(0)), "norm2est(): parameter 'tolerance' must be > 0" );
  arma_debug_check( (max_iter  == uword(0)), "norm2est(): parameter 'max_iter' must be > 0"  );
  
  const T tol = (tolerance == T(0)) ? T(1e-6) : T(tolerance);
  
  const unwrap_spmat<T1> U(X.get_ref());
  const SpMat<eT>& A   = U.M;
  
  if(A.n_nonzero == 0)  { return T(0); }
  
  if(A.internal_has_nonfinite())  { arma_debug_warn_level(1, "norm2est(): given matrix has non-finite elements"); }
  
  if((A.n_rows == 1) || (A.n_cols == 1))  { return spop_norm::vec_norm_k(A.values, A.n_nonzero, 2); }
  
  norm2est_randu_filler<eT> randu_filler;
  
  Mat<eT> x(A.n_rows, 1, fill::none);
  Mat<eT> y(A.n_cols, 1, fill::none);
  
  randu_filler.fill(y.memptr(), y.n_elem);
  
  T est_old = 0;
  T est_cur = 0;
  
  for(uword i=0; i<max_iter; ++i)
    {
    arma_extra_debug_print(arma_str::format("norm2est(): iteration: %u") % i);
    
    x = A * y;
    
    T x_norm = op_norm::vec_norm_2( Proxy< Mat<eT> >(x) );
    
    if(x_norm == T(0) || (arma_isfinite(x_norm) == false) || (x.internal_has_nonfinite()))
      {
      randu_filler.fill(x.memptr(), x.n_elem);
      
      x_norm = op_norm::vec_norm_2( Proxy< Mat<eT> >(x) );
      }
    
    if(x_norm != T(0))  { x /= x_norm; }
    
    y = A.t() * x;
    
    est_old = est_cur;
    est_cur = op_norm::vec_norm_2( Proxy< Mat<eT> >(y) );
    
    arma_extra_debug_print(arma_str::format("norm2est(): est_old: %e") % est_old);
    arma_extra_debug_print(arma_str::format("norm2est(): est_cur: %e") % est_cur);
    
    if(arma_isfinite(est_cur) == false)  { return est_old; }
    
    if( ((std::abs)(est_cur - est_old)) <= (tol * (std::max)(est_cur,est_old)) )  { break; }
    }
  
  return est_cur;
  }



//! @}
