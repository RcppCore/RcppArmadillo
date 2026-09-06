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



//! \addtogroup glue_powext
//! @{


template<typename T1, typename T2>
inline
void
glue_powext::apply(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_powext>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UA(X.A);
  const quasi_unwrap<T2> UB(X.B);
  
  const Mat<eT>& A = UA.M;
  const Mat<eT>& B = UB.M;
  
  arma_conform_assert_same_size(A, B, "element-wise pow()");
  
  const bool UA_bad_alias = UA.is_alias(out) && (UA.has_subview);  // allow inplace operation
  const bool UB_bad_alias = UB.is_alias(out);
  
  if(UA_bad_alias || UB_bad_alias)
    {
    Mat<eT> tmp;
    
    glue_powext::apply(tmp, A, B);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_powext::apply(out, A, B);
    }
  }



template<typename eT>
inline
void
glue_powext::apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_debug_sigprint();
  
  out.set_size(A.n_rows, A.n_cols);
  
  const uword N = out.n_elem;
  
        eT* out_mem = out.memptr();
  const eT*   A_mem =   A.memptr();
  const eT*   B_mem =   B.memptr();
  
  if( arma_config::openmp && mp_gate<eT>::eval(N) )
    {
    #if defined(ARMA_USE_OPENMP)
      {
      const int n_threads = mp_thread_limit::get();
      
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i<N; ++i)
        {
        out_mem[i] = eop_aux::pow(A_mem[i], B_mem[i]);
        }
      }
    #endif
    }
  else
    {
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = eop_aux::pow(A_mem[i], B_mem[i]);
      }
    }
  }



template<typename T1, typename T2>
inline
void
glue_powext::apply(Cube<typename T1::elem_type>& out, const GlueCube<T1, T2, glue_powext>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> UA(X.A);
  const unwrap_cube<T2> UB(X.B);
  
  const Cube<eT>& A = UA.M;
  const Cube<eT>& B = UB.M;
  
  arma_conform_assert_same_size(A, B, "element-wise pow()");
  
  if(UB.is_alias(out))
    {
    Cube<eT> tmp;
    
    glue_powext::apply(tmp, A, B);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_powext::apply(out, A, B);
    }
  }



template<typename eT>
inline
void
glue_powext::apply(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B)
  {
  arma_debug_sigprint();
  
  out.set_size(A.n_rows, A.n_cols, A.n_slices);
  
  const uword N = out.n_elem;
  
        eT* out_mem = out.memptr();
  const eT*   A_mem =   A.memptr();
  const eT*   B_mem =   B.memptr();
  
  if( arma_config::openmp && mp_gate<eT>::eval(N) )
    {
    #if defined(ARMA_USE_OPENMP)
      {
      const int n_threads = mp_thread_limit::get();
      
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i<N; ++i)
        {
        out_mem[i] = eop_aux::pow(A_mem[i], B_mem[i]);
        }
      }
    #endif
    }
  else
    {
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = eop_aux::pow(A_mem[i], B_mem[i]);
      }
    }
  }



//



template<typename T1, typename T2>
inline
void
glue_powext_cx::apply(Mat<typename T1::elem_type>& out, const mtGlue<typename T1::elem_type, T1, T2, glue_powext_cx>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const quasi_unwrap<T1> UA(X.A);
  const quasi_unwrap<T2> UB(X.B);
  
  const Mat<eT>& A = UA.M;
  const Mat< T>& B = UB.M;
  
  arma_conform_assert_same_size(A, B, "element-wise pow()");
  
  if(UA.is_alias(out) && (UA.has_subview))
    {
    Mat<eT> tmp;
    
    glue_powext_cx::apply(tmp, A, B);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_powext_cx::apply(out, A, B);
    }
  }



template<typename T>
inline
void
glue_powext_cx::apply(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A, const Mat<T>& B)
  {
  arma_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  out.set_size(A.n_rows, A.n_cols);
  
  const uword N = out.n_elem;
  
        eT* out_mem = out.memptr();
  const eT*   A_mem =   A.memptr();
  const  T*   B_mem =   B.memptr();
  
  if( arma_config::openmp && mp_gate<eT>::eval(N) )
    {
    #if defined(ARMA_USE_OPENMP)
      {
      const int n_threads = mp_thread_limit::get();
      
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i<N; ++i)
        {
        out_mem[i] = std::pow(A_mem[i], B_mem[i]);
        }
      }
    #endif
    }
  else
    {
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = std::pow(A_mem[i], B_mem[i]);
      }
    }
  }



template<typename T1, typename T2>
inline
void
glue_powext_cx::apply(Cube<typename T1::elem_type>& out, const mtGlueCube<typename T1::elem_type,T1,T2,glue_powext_cx>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef typename get_pod_type<eT>::result T;
  
  const unwrap_cube<T1> UA(X.A);
  const unwrap_cube<T2> UB(X.B);
  
  const Cube<eT>& A = UA.M;
  const Cube< T>& B = UB.M;
  
  arma_conform_assert_same_size(A, B, "element-wise pow()");
  
  glue_powext_cx::apply(out, A, B);
  }



template<typename T>
inline
void
glue_powext_cx::apply(Cube< std::complex<T> >& out, const Cube< std::complex<T> >& A, const Cube<T>& B)
  {
  arma_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  out.set_size(A.n_rows, A.n_cols, A.n_slices);
  
  const uword N = out.n_elem;
  
        eT* out_mem = out.memptr();
  const eT*   A_mem =   A.memptr();
  const  T*   B_mem =   B.memptr();
  
  if( arma_config::openmp && mp_gate<eT>::eval(N) )
    {
    #if defined(ARMA_USE_OPENMP)
      {
      const int n_threads = mp_thread_limit::get();
      
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i<N; ++i)
        {
        out_mem[i] = std::pow(A_mem[i], B_mem[i]);
        }
      }
    #endif
    }
  else
    {
    for(uword i=0; i<N; ++i)
      {
      out_mem[i] = std::pow(A_mem[i], B_mem[i]);
      }
    }
  }



//! @}
