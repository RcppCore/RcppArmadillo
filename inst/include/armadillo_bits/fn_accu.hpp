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


//! \addtogroup fn_accu
//! @{



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_proxy_linear(const Proxy<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  eT val = eT(0);
  
  typename Proxy<T1>::ea_type Pea = P.get_ea();
  
  const uword n_elem = P.get_n_elem();
  
  if( arma_config::openmp && Proxy<T1>::use_mp && mp_gate<eT>::eval(n_elem) )
    {
    #if defined(ARMA_USE_OPENMP)
      {
      // NOTE: using parallelisation with manual reduction workaround to take into account complex numbers;
      // NOTE: OpenMP versions lower than 4.0 do not support user-defined reduction
      
      const int   n_threads_max = mp_thread_limit::get();
      const uword n_threads_use = (std::min)(uword(podarray_prealloc_n_elem::val), uword(n_threads_max));
      const uword chunk_size    = n_elem / n_threads_use;
      
      podarray<eT> partial_accs(n_threads_use);
      
      #pragma omp parallel for schedule(static) num_threads(int(n_threads_use))
      for(uword thread_id=0; thread_id < n_threads_use; ++thread_id)
        {
        const uword start = (thread_id+0) * chunk_size;
        const uword endp1 = (thread_id+1) * chunk_size;
        
        eT acc = eT(0);
        for(uword i=start; i < endp1; ++i)  { acc += Pea[i]; }
        
        partial_accs[thread_id] = acc;
        }
      
      for(uword thread_id=0; thread_id < n_threads_use; ++thread_id)  { val += partial_accs[thread_id]; }
      
      for(uword i=(n_threads_use*chunk_size); i < n_elem; ++i)  { val += Pea[i]; }
      }
    #endif
    }
  else
    {
    #if defined(__FAST_MATH__)
      {
      if(P.is_aligned())
        {
        typename Proxy<T1>::aligned_ea_type Pea_aligned = P.get_aligned_ea();
        
        for(uword i=0; i<n_elem; ++i)  { val += Pea_aligned.at_alt(i); }
        }
      else
        {
        for(uword i=0; i<n_elem; ++i)  { val += Pea[i]; }
        }
      }
    #else
      {
      eT val1 = eT(0);
      eT val2 = eT(0);
      
      uword i,j;
      for(i=0, j=1; j < n_elem; i+=2, j+=2)  { val1 += Pea[i]; val2 += Pea[j]; }
      
      if(i < n_elem)  { val1 += Pea[i]; }
      
      val = val1 + val2;
      }
    #endif
    }
  
  return val;
  }



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_proxy_at_mp(const Proxy<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  eT val = eT(0);
  
  #if defined(ARMA_USE_OPENMP)
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_cols == 1)
      {
      const int   n_threads_max = mp_thread_limit::get();
      const uword n_threads_use = (std::min)(uword(podarray_prealloc_n_elem::val), uword(n_threads_max));
      const uword chunk_size    = n_rows / n_threads_use;
      
      podarray<eT> partial_accs(n_threads_use);
      
      #pragma omp parallel for schedule(static) num_threads(int(n_threads_use))
      for(uword thread_id=0; thread_id < n_threads_use; ++thread_id)
        {
        const uword start = (thread_id+0) * chunk_size;
        const uword endp1 = (thread_id+1) * chunk_size;
        
        eT acc = eT(0);
        for(uword i=start; i < endp1; ++i)  { acc += P.at(i,0); }
        
        partial_accs[thread_id] = acc;
        }
      
      for(uword thread_id=0; thread_id < n_threads_use; ++thread_id)  { val += partial_accs[thread_id]; }
      
      for(uword i=(n_threads_use*chunk_size); i < n_rows; ++i)  { val += P.at(i,0); }
      }
    else
    if(n_rows == 1)
      {
      const int   n_threads_max = mp_thread_limit::get();
      const uword n_threads_use = (std::min)(uword(podarray_prealloc_n_elem::val), uword(n_threads_max));
      const uword chunk_size    = n_cols / n_threads_use;
      
      podarray<eT> partial_accs(n_threads_use);
      
      #pragma omp parallel for schedule(static) num_threads(int(n_threads_use))
      for(uword thread_id=0; thread_id < n_threads_use; ++thread_id)
        {
        const uword start = (thread_id+0) * chunk_size;
        const uword endp1 = (thread_id+1) * chunk_size;
        
        eT acc = eT(0);
        for(uword i=start; i < endp1; ++i)  { acc += P.at(0,i); }
        
        partial_accs[thread_id] = acc;
        }
      
      for(uword thread_id=0; thread_id < n_threads_use; ++thread_id)  { val += partial_accs[thread_id]; }
      
      for(uword i=(n_threads_use*chunk_size); i < n_cols; ++i)  { val += P.at(0,i); }
      }
    else
      {
      podarray<eT> col_accs(n_cols);
      
      const int n_threads = mp_thread_limit::get();
      
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword col=0; col < n_cols; ++col)
        {
        eT val1 = eT(0);
        eT val2 = eT(0);
        
        uword i,j;
        for(i=0, j=1; j < n_rows; i+=2, j+=2)  { val1 += P.at(i,col); val2 += P.at(j,col); }
        
        if(i < n_rows)  { val1 += P.at(i,col); }
        
        col_accs[col] = val1 + val2;
        }
      
      val = arrayops::accumulate(col_accs.memptr(), n_cols);
      }
    }
  #else
    {
    arma_ignore(P);
    }
  #endif
  
  return val;
  }



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_proxy_at(const Proxy<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(arma_config::openmp && Proxy<T1>::use_mp && mp_gate<eT>::eval(P.get_n_elem()))
    {
    return accu_proxy_at_mp(P);
    }
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  eT val = eT(0);
  
  if(n_rows != 1)
    {
    eT val1 = eT(0);
    eT val2 = eT(0);
    
    for(uword col=0; col < n_cols; ++col)
      {
      uword i,j;
      for(i=0, j=1; j < n_rows; i+=2, j+=2)  { val1 += P.at(i,col); val2 += P.at(j,col); }
      
      if(i < n_rows)  { val1 += P.at(i,col); }
      }
    
    val = val1 + val2;
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)  { val += P.at(0,col); }
    }
  
  return val;
  }



//! accumulate the elements of a matrix
template<typename T1>
arma_warn_unused
arma_hot
inline
typename enable_if2< is_arma_type<T1>::value, typename T1::elem_type >::result
accu(const T1& X)
  {
  arma_debug_sigprint();
  
  if((is_Mat<T1>::value) || (is_subview_col<T1>::value) || (is_Mat<typename Proxy<T1>::stored_type>::value))
    {
    const quasi_unwrap<T1> U(X);
    
    return arrayops::accumulate(U.M.memptr(), U.M.n_elem);
    }
  
  const Proxy<T1> P(X);
  
  return (Proxy<T1>::use_at) ? accu_proxy_at(P) : accu_proxy_linear(P);
  }



template<typename T1, typename functor>
inline
typename T1::elem_type
accu_op_omit_helper(const Proxy<T1>& P, functor is_omitted)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  constexpr eT eT_zero = eT(0);
  
  eT acc = eT(0);
  
  if(Proxy<T1>::use_at)
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    for(uword c=0; c < n_cols; ++c)
    for(uword r=0; r < n_rows; ++r)
      {
      const eT val = P.at(r,c);
      
      acc += is_omitted(val) ? eT_zero : val;
      }
    }
  else
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
    
    const uword n_elem = P.get_n_elem();
    
    eT val1 = eT(0);
    eT val2 = eT(0);
    
    uword i,j;
    for(i=0, j=1; j < n_elem; i+=2, j+=2)
      {
      const eT tmp_i = Pea[i];
      const eT tmp_j = Pea[j];
      
      val1 += is_omitted(tmp_i) ? eT_zero : tmp_i;
      val2 += is_omitted(tmp_j) ? eT_zero : tmp_j;
      }
    
    if(i < n_elem)
      {
      const eT tmp_i = Pea[i];
      
      val1 += is_omitted(tmp_i) ? eT_zero : tmp_i;
      }
    
    acc = val1 + val2;
    }
  
  return acc;
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const Op<T1, op_omit>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword omit_mode = in.aux_uword_a;
  
  if(arma_config::fast_math_warn)
    {
    if(omit_mode == 1)  { arma_warn(1, "omit_nan(): detection of NaN is not reliable in fast math mode"); }
    if(omit_mode == 2)  { arma_warn(1, "omit_nonfinite(): detection of non-finite values is not reliable in fast math mode"); }
    }
  
  auto is_omitted_1 = [](const eT& x) -> bool  { return arma_isnan(x);       };
  auto is_omitted_2 = [](const eT& x) -> bool  { return arma_isnonfinite(x); };
  
  const Proxy<T1> P(in.m);
  
  eT acc = eT(0);
  
  if(omit_mode == 1)  { acc = accu_op_omit_helper(P, is_omitted_1); }
  if(omit_mode == 2)  { acc = accu_op_omit_helper(P, is_omitted_2); }
  
  return acc;
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const eOp<T1,eop_square>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef eOp<T1,eop_square> expr_type;
  
  typedef typename expr_type::proxy_type::stored_type expr_P_stored_type;
  
  if((is_Mat<expr_P_stored_type>::value) || (is_subview_col<expr_P_stored_type>::value))
    {
    const quasi_unwrap<expr_P_stored_type> U(expr.P.Q);
    
    const eT* X_mem = U.M.memptr();
    
    return op_dot::direct_dot(U.M.n_elem, X_mem, X_mem);
    }
  
  const Proxy<expr_type> P(expr);
  
  return (Proxy<expr_type>::use_at) ? accu_proxy_at(P) : accu_proxy_linear(P);
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const eOp<T1,eop_pow>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef eOp<T1,eop_pow> expr_type;
  
  if(arma_config::optimise_powexpr && (expr.aux == eT(2)))
    {
    typedef eOp<T1,eop_square> modified_expr_type;
    
    return accu( reinterpret_cast< const modified_expr_type& >(expr) );
    }
  
  if(arma_config::optimise_powexpr && (expr.aux == eT(0.5)) && is_non_integral<eT>::value)
    {
    typedef eOp<T1,eop_sqrt> modified_expr_type;
    
    return accu( reinterpret_cast< const modified_expr_type& >(expr) );
    }
  
  const Proxy<expr_type> P(expr);
  
  return (Proxy<expr_type>::use_at) ? accu_proxy_at(P) : accu_proxy_linear(P);
  }



//! explicit handling of multiply-and-accumulate
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
accu(const eGlue<T1,T2,eglue_schur>& expr)
  {
  arma_debug_sigprint();
  
  typedef eGlue<T1,T2,eglue_schur> expr_type;
  
  typedef typename expr_type::proxy1_type::stored_type P1_stored_type;
  typedef typename expr_type::proxy2_type::stored_type P2_stored_type;
  
  constexpr bool is_sv = (is_subview<P1_stored_type>::value) || (is_subview<P2_stored_type>::value);
  
  if( (is_sv) && (expr.get_n_rows() >= 4) )
    {
    arma_debug_print("accu(): eglue_schur subview optimisation");
    
    typedef typename T1::elem_type eT;
    
    const sv_keep_unwrap<P1_stored_type>& UA(expr.P1.Q);
    const sv_keep_unwrap<P2_stored_type>& UB(expr.P2.Q);
    
    typedef typename sv_keep_unwrap<T1>::stored_type UA_M_type;
    typedef typename sv_keep_unwrap<T2>::stored_type UB_M_type;
    
    const UA_M_type& A = UA.M;
    const UB_M_type& B = UB.M;
    
    // A and B have the same size (checked by the eGlue constructor)
    
    const uword A_n_rows = A.n_rows;
    const uword A_n_cols = A.n_cols;
    
    eT acc = eT(0);
    
    for(uword c=0; c < A_n_cols; ++c)  { acc += op_dot::direct_dot(A_n_rows, A.colptr(c), B.colptr(c)); }
    
    return acc;
    }
  
  constexpr bool have_direct_mem_1 = (is_Mat<P1_stored_type>::value) || (is_subview_col<P1_stored_type>::value);
  constexpr bool have_direct_mem_2 = (is_Mat<P2_stored_type>::value) || (is_subview_col<P2_stored_type>::value);
  
  if(have_direct_mem_1 && have_direct_mem_2)
    {
    arma_debug_print("accu(): eglue_schur direct_mem optimisation");
    
    const quasi_unwrap<P1_stored_type> tmp1(expr.P1.Q);
    const quasi_unwrap<P2_stored_type> tmp2(expr.P2.Q);
    
    return op_dot::direct_dot(tmp1.M.n_elem, tmp1.M.memptr(), tmp2.M.memptr());
    }
  
  const Proxy<expr_type> P(expr);
  
  return (Proxy<expr_type>::use_at) ? accu_proxy_at(P) : accu_proxy_linear(P);
  }



template<typename T1, typename op_type>
arma_warn_unused
inline
uword
accu(const mtOp<uword,T1,op_type>& X, const typename arma_op_rel_only<op_type>::result* junk1 = nullptr, const typename arma_not_cx<typename T1::elem_type>::result* junk2 = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const Proxy<T1> P(X.m);
  
  uword count = 0;
  
  if(Proxy<T1>::use_at == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
          ea_type A      = P.get_ea();
    const uword   n_elem = P.get_n_elem();
    
    for(uword i=0; i<n_elem; ++i)
      {
      const eT val = A[i];
      
      bool condition;
      
           if(is_same_type<op_type, op_rel_eq       >::yes)  { condition = (val == k  ); }
      else if(is_same_type<op_type, op_rel_noteq    >::yes)  { condition = (val != k  ); }
      else if(is_same_type<op_type, op_rel_lt_pre   >::yes)  { condition = (k   <  val); }
      else if(is_same_type<op_type, op_rel_lt_post  >::yes)  { condition = (val <  k  ); }
      else if(is_same_type<op_type, op_rel_gt_pre   >::yes)  { condition = (k   >  val); }
      else if(is_same_type<op_type, op_rel_gt_post  >::yes)  { condition = (val >  k  ); }
      else if(is_same_type<op_type, op_rel_lteq_pre >::yes)  { condition = (k   <= val); }
      else if(is_same_type<op_type, op_rel_lteq_post>::yes)  { condition = (val <= k  ); }
      else if(is_same_type<op_type, op_rel_gteq_pre >::yes)  { condition = (k   >= val); }
      else if(is_same_type<op_type, op_rel_gteq_post>::yes)  { condition = (val >= k  ); }
      else { condition = false; }
      
      count += (condition) ? uword(1) : uword(0);
      }
    }
  else
    {
    const uword P_n_cols = P.get_n_cols();
    const uword P_n_rows = P.get_n_rows();
    
    for(uword col=0; col < P_n_cols; ++col)
    for(uword row=0; row < P_n_rows; ++row)
      {
      const eT val = P.at(row,col);
      
      bool condition;
      
           if(is_same_type<op_type, op_rel_eq       >::yes)  { condition = (val == k  ); }
      else if(is_same_type<op_type, op_rel_noteq    >::yes)  { condition = (val != k  ); }
      else if(is_same_type<op_type, op_rel_lt_pre   >::yes)  { condition = (k   <  val); }
      else if(is_same_type<op_type, op_rel_lt_post  >::yes)  { condition = (val <  k  ); }
      else if(is_same_type<op_type, op_rel_gt_pre   >::yes)  { condition = (k   >  val); }
      else if(is_same_type<op_type, op_rel_gt_post  >::yes)  { condition = (val >  k  ); }
      else if(is_same_type<op_type, op_rel_lteq_pre >::yes)  { condition = (k   <= val); }
      else if(is_same_type<op_type, op_rel_lteq_post>::yes)  { condition = (val <= k  ); }
      else if(is_same_type<op_type, op_rel_gteq_pre >::yes)  { condition = (k   >= val); }
      else if(is_same_type<op_type, op_rel_gteq_post>::yes)  { condition = (val >= k  ); }
      else { condition = false; }
      
      count += (condition) ? uword(1) : uword(0);
      }
    }
  
  return count;
  }



template<typename T1, typename op_type>
arma_warn_unused
inline
uword
accu(const mtOp<uword,T1,op_type>& X, const typename arma_op_rel_only<op_type>::result* junk1 = nullptr, const typename arma_cx_only<typename T1::elem_type>::result* junk2 = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const Proxy<T1> P(X.m);
  
  uword count = 0;
  
  if(Proxy<T1>::use_at == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
          ea_type A      = P.get_ea();
    const uword   n_elem = P.get_n_elem();
    
    for(uword i=0; i<n_elem; ++i)
      {
      const eT val = A[i];
      
      bool condition;
      
           if(is_same_type<op_type, op_rel_eq   >::yes)  { condition = (val == k); }
      else if(is_same_type<op_type, op_rel_noteq>::yes)  { condition = (val != k); }
      else { condition = false; }
      
      count += (condition) ? uword(1) : uword(0);
      }
    }
  else
    {
    const uword P_n_cols = P.get_n_cols();
    const uword P_n_rows = P.get_n_rows();
    
    for(uword col=0; col < P_n_cols; ++col)
    for(uword row=0; row < P_n_rows; ++row)
      {
      const eT val = P.at(row,col);
      
      bool condition;
      
           if(is_same_type<op_type, op_rel_eq   >::yes)  { condition = (val == k); }
      else if(is_same_type<op_type, op_rel_noteq>::yes)  { condition = (val != k); }
      else { condition = false; }
      
      count += (condition) ? uword(1) : uword(0);
      }
    }
  
  return count;
  }



template<typename T1, typename T2>
arma_warn_unused
inline
uword
accu(const mtGlue<uword,T1,T2,glue_rel_noteq>& X)
  {
  arma_debug_sigprint();
  
  const Proxy<T1> PA(X.A);
  const Proxy<T2> PB(X.B);
  
  arma_conform_assert_same_size(PA, PB, "operator!=");
  
  uword n_nonzero = 0;
  
  if( (Proxy<T1>::use_at == false) && (Proxy<T2>::use_at == false) )
    {
    typedef typename Proxy<T1>::ea_type PA_ea_type;
    typedef typename Proxy<T2>::ea_type PB_ea_type;
    
          PA_ea_type A      = PA.get_ea();
          PB_ea_type B      = PB.get_ea();
    const uword      n_elem = PA.get_n_elem();
    
    for(uword i=0; i < n_elem; ++i)
      {
      n_nonzero += (A[i] != B[i]) ? uword(1) : uword(0);
      }
    }
  else
    {
    const uword PA_n_cols = PA.get_n_cols();
    const uword PA_n_rows = PA.get_n_rows();
    
    if(PA_n_rows == 1)
      {
      for(uword col=0; col < PA_n_cols; ++col)
        {
        n_nonzero += (PA.at(0,col) != PB.at(0,col)) ? uword(1) : uword(0);
        }
      }
    else
      {
      for(uword col=0; col < PA_n_cols; ++col)
      for(uword row=0; row < PA_n_rows; ++row)
        {
        n_nonzero += (PA.at(row,col) != PB.at(row,col)) ? uword(1) : uword(0);
        }
      }
    }
  
  return n_nonzero;
  }



template<typename T1, typename T2>
arma_warn_unused
inline
uword
accu(const mtGlue<uword,T1,T2,glue_rel_eq>& X)
  {
  arma_debug_sigprint();
  
  const Proxy<T1> PA(X.A);
  const Proxy<T2> PB(X.B);
  
  arma_conform_assert_same_size(PA, PB, "operator==");
  
  uword n_nonzero = 0;
  
  if( (Proxy<T1>::use_at == false) && (Proxy<T2>::use_at == false) )
    {
    typedef typename Proxy<T1>::ea_type PA_ea_type;
    typedef typename Proxy<T2>::ea_type PB_ea_type;
    
          PA_ea_type A      = PA.get_ea();
          PB_ea_type B      = PB.get_ea();
    const uword      n_elem = PA.get_n_elem();
    
    for(uword i=0; i < n_elem; ++i)
      {
      n_nonzero += (A[i] == B[i]) ? uword(1) : uword(0);
      }
    }
  else
    {
    const uword PA_n_cols = PA.get_n_cols();
    const uword PA_n_rows = PA.get_n_rows();
    
    if(PA_n_rows == 1)
      {
      for(uword col=0; col < PA_n_cols; ++col)
        {
        n_nonzero += (PA.at(0,col) == PB.at(0,col)) ? uword(1) : uword(0);
        }
      }
    else
      {
      for(uword col=0; col < PA_n_cols; ++col)
      for(uword row=0; row < PA_n_rows; ++row)
        {
        n_nonzero += (PA.at(row,col) == PB.at(row,col)) ? uword(1) : uword(0);
        }
      }
    }
  
  return n_nonzero;
  }



//! accumulate the elements of a subview (submatrix)
template<typename eT>
arma_warn_unused
arma_hot
inline
eT
accu(const subview<eT>& X)
  {
  arma_debug_sigprint();  
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(X_n_rows == 1)
    {
    const Mat<eT>& m = X.m;
    
    const uword col_offset = X.aux_col1;
    const uword row_offset = X.aux_row1;
    
    eT val1 = eT(0);
    eT val2 = eT(0);
    
    uword i,j;
    for(i=0, j=1; j < X_n_cols; i+=2, j+=2)
      {
      val1 += m.at(row_offset, col_offset + i);
      val2 += m.at(row_offset, col_offset + j);
      }
    
    if(i < X_n_cols)  { val1 += m.at(row_offset, col_offset + i); }
    
    return val1 + val2;
    }
  
  if(X_n_cols == 1)  { return arrayops::accumulate( X.colptr(0), X_n_rows ); }
  
  eT val = eT(0);
  
  for(uword col=0; col < X_n_cols; ++col)
    {
    val += arrayops::accumulate( X.colptr(col), X_n_rows );
    }
  
  return val;
  }



template<typename eT>
arma_warn_unused
arma_hot
inline
eT
accu(const subview_col<eT>& X)
  {
  arma_debug_sigprint();  
  
  return arrayops::accumulate( X.colmem, X.n_rows );
  }



//



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_cube_proxy_linear(const ProxyCube<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  eT val = eT(0);
    
  typename ProxyCube<T1>::ea_type Pea = P.get_ea();
  
  const uword n_elem = P.get_n_elem();
  
  if( arma_config::openmp && ProxyCube<T1>::use_mp && mp_gate<eT>::eval(n_elem) )
    {
    #if defined(ARMA_USE_OPENMP)
      {
      // NOTE: using parallelisation with manual reduction workaround to take into account complex numbers;
      // NOTE: OpenMP versions lower than 4.0 do not support user-defined reduction
      
      const int   n_threads_max = mp_thread_limit::get();
      const uword n_threads_use = (std::min)(uword(podarray_prealloc_n_elem::val), uword(n_threads_max));
      const uword chunk_size    = n_elem / n_threads_use;
      
      podarray<eT> partial_accs(n_threads_use);
      
      #pragma omp parallel for schedule(static) num_threads(int(n_threads_use))
      for(uword thread_id=0; thread_id < n_threads_use; ++thread_id)
        {
        const uword start = (thread_id+0) * chunk_size;
        const uword endp1 = (thread_id+1) * chunk_size;
        
        eT acc = eT(0);
        for(uword i=start; i < endp1; ++i)  { acc += Pea[i]; }
        
        partial_accs[thread_id] = acc;
        }
      
      for(uword thread_id=0; thread_id < n_threads_use; ++thread_id)  { val += partial_accs[thread_id]; }
      
      for(uword i=(n_threads_use*chunk_size); i < n_elem; ++i)  { val += Pea[i]; }
      }
    #endif
    }
  else
    {
    #if defined(__FAST_MATH__)
      {
      if(P.is_aligned())
        {
        typename ProxyCube<T1>::aligned_ea_type Pea_aligned = P.get_aligned_ea();
        
        for(uword i=0; i<n_elem; ++i)  { val += Pea_aligned.at_alt(i); }
        }
      else
        {
        for(uword i=0; i<n_elem; ++i)  { val += Pea[i]; }
        }
      }
    #else
      {
      eT val1 = eT(0);
      eT val2 = eT(0);
      
      uword i,j;
      for(i=0, j=1; j<n_elem; i+=2, j+=2)  { val1 += Pea[i]; val2 += Pea[j]; }
      
      if(i < n_elem) { val1 += Pea[i]; }
      
      val = val1 + val2;
      }
    #endif
    }
  
  return val;
  }



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_cube_proxy_at_mp(const ProxyCube<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  eT val = eT(0);
  
  #if defined(ARMA_USE_OPENMP)
    {
    const uword n_rows   = P.get_n_rows();
    const uword n_cols   = P.get_n_cols();
    const uword n_slices = P.get_n_slices();
    
    podarray<eT> slice_accs(n_slices);
    
    const int n_threads = mp_thread_limit::get();
    
    #pragma omp parallel for schedule(static) num_threads(n_threads)
    for(uword slice = 0; slice < n_slices; ++slice)
      {
      eT val1 = eT(0);
      eT val2 = eT(0);
      
      for(uword col = 0; col < n_cols; ++col)
        {
        uword i,j;
        for(i=0, j=1; j<n_rows; i+=2, j+=2)  { val1 += P.at(i,col,slice);  val2 += P.at(j,col,slice); }
        
        if(i < n_rows)  { val1 += P.at(i,col,slice); }
        }
      
      slice_accs[slice] = val1 + val2;
      }
    
    val = arrayops::accumulate(slice_accs.memptr(), slice_accs.n_elem);
    }
  #else
    {
    arma_ignore(P);
    }
  #endif
  
  return val;
  }



template<typename T1>
arma_hot
inline
typename T1::elem_type
accu_cube_proxy_at(const ProxyCube<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(arma_config::openmp && ProxyCube<T1>::use_mp && mp_gate<eT>::eval(P.get_n_elem()))
    {
    return accu_cube_proxy_at_mp(P);
    }
  
  const uword n_rows   = P.get_n_rows();
  const uword n_cols   = P.get_n_cols();
  const uword n_slices = P.get_n_slices();
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  for(uword slice = 0; slice < n_slices; ++slice)
  for(uword col   = 0; col   < n_cols;   ++col  )
    {
    uword i,j;
    for(i=0, j=1; j<n_rows; i+=2, j+=2)  { val1 += P.at(i,col,slice); val2 += P.at(j,col,slice); }
    
    if(i < n_rows)  { val1 += P.at(i,col,slice); }
    }
  
  return (val1 + val2);
  }



//! accumulate the elements of a cube
template<typename T1>
arma_warn_unused
arma_hot
inline
typename T1::elem_type
accu(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_debug_sigprint();
  
  if((is_Cube<T1>::value) || (is_Cube<typename ProxyCube<T1>::stored_type>::value))
    {
    const unwrap_cube<T1> U(X.get_ref());
    
    return arrayops::accumulate(U.M.memptr(), U.M.n_elem);
    }
  
  const ProxyCube<T1> P(X.get_ref());
  
  return (ProxyCube<T1>::use_at) ? accu_cube_proxy_at(P) : accu_cube_proxy_linear(P);
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const eOpCube<T1,eop_square>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef eOpCube<T1,eop_square> expr_type;
  
  typedef typename expr_type::proxy_type::stored_type expr_P_stored_type;
  
  if(is_Cube<expr_P_stored_type>::value)
    {
    const unwrap_cube<expr_P_stored_type> U(expr.P.Q);
    
    const eT* X_mem = U.M.memptr();
    
    return op_dot::direct_dot(U.M.n_elem, X_mem, X_mem);
    }
  
  const ProxyCube<expr_type> P(expr);
  
  return (ProxyCube<expr_type>::use_at) ? accu_cube_proxy_at(P) : accu_cube_proxy_linear(P);
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const eOpCube<T1,eop_pow>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef eOpCube<T1,eop_pow> expr_type;
  
  if(arma_config::optimise_powexpr && (expr.aux == eT(2)))
    {
    typedef eOpCube<T1,eop_square> modified_expr_type;
    
    return accu( reinterpret_cast< const modified_expr_type& >(expr) );
    }
  
  if(arma_config::optimise_powexpr && (expr.aux == eT(0.5)) && is_non_integral<eT>::value)
    {
    typedef eOpCube<T1,eop_sqrt> modified_expr_type;
    
    return accu( reinterpret_cast< const modified_expr_type& >(expr) );
    }
  
  const ProxyCube<expr_type> P(expr);
  
  return (ProxyCube<expr_type>::use_at) ? accu_cube_proxy_at(P) : accu_cube_proxy_linear(P);
  }



//! explicit handling of multiply-and-accumulate (cube version)
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
accu(const eGlueCube<T1,T2,eglue_schur>& expr)
  {
  arma_debug_sigprint();
  
  typedef eGlueCube<T1,T2,eglue_schur> expr_type;
  
  typedef typename expr_type::proxy1_type::stored_type P1_stored_type;
  typedef typename expr_type::proxy2_type::stored_type P2_stored_type;
  
  if(is_Cube<P1_stored_type>::value && is_Cube<P2_stored_type>::value)
    {
    const unwrap_cube<P1_stored_type> tmp1(expr.P1.Q);
    const unwrap_cube<P2_stored_type> tmp2(expr.P2.Q);
    
    return op_dot::direct_dot(tmp1.M.n_elem, tmp1.M.memptr(), tmp2.M.memptr());
    }
  
  const ProxyCube<expr_type> P(expr);
  
  return (ProxyCube<expr_type>::use_at) ? accu_cube_proxy_at(P) : accu_cube_proxy_linear(P);
  }



template<typename T1, typename functor>
inline
typename T1::elem_type
accu_cube_omit_helper(const ProxyCube<T1>& P, functor is_omitted)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  constexpr eT eT_zero = eT(0);
  
  eT acc = eT(0);
  
  if(ProxyCube<T1>::use_at)
    {
    const uword n_r = P.get_n_rows();
    const uword n_c = P.get_n_cols();
    const uword n_s = P.get_n_slices();
    
    for(uword s=0; s < n_s; ++s)
    for(uword c=0; c < n_c; ++c)
    for(uword r=0; r < n_r; ++r)
      {
      const eT val = P.at(r,c,s);
      
      acc += is_omitted(val) ? eT_zero : val;
      }
    }
  else
    {
    typename ProxyCube<T1>::ea_type Pea = P.get_ea();
    
    const uword n_elem = P.get_n_elem();
    
    eT val1 = eT(0);
    eT val2 = eT(0);
    
    uword i,j;
    for(i=0, j=1; j < n_elem; i+=2, j+=2)
      {
      const eT tmp_i = Pea[i];
      const eT tmp_j = Pea[j];
      
      val1 += is_omitted(tmp_i) ? eT_zero : tmp_i;
      val2 += is_omitted(tmp_j) ? eT_zero : tmp_j;
      }
    
    if(i < n_elem)
      {
      const eT tmp_i = Pea[i];
      
      val1 += is_omitted(tmp_i) ? eT_zero : tmp_i;
      }
    
    acc = val1 + val2;
    }
  
  return acc;
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const CubeToMatOp<T1, op_omit_cube>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const ProxyCube<T1> P(in.m);
  
  const uword omit_mode = in.aux_uword;
  
  if(arma_config::fast_math_warn)
    {
    if(omit_mode == 1)  { arma_warn(1, "omit_nan(): detection of NaN is not reliable in fast math mode"); }
    if(omit_mode == 2)  { arma_warn(1, "omit_nonfinite(): detection of non-finite values is not reliable in fast math mode"); }
    }
  
  auto is_omitted_1 = [](const eT& x) -> bool  { return arma_isnan(x);       };
  auto is_omitted_2 = [](const eT& x) -> bool  { return arma_isnonfinite(x); };
  
  eT acc = eT(0);
  
  if(omit_mode == 1)  { acc = accu_cube_omit_helper(P, is_omitted_1); }
  if(omit_mode == 2)  { acc = accu_cube_omit_helper(P, is_omitted_2); }
  
  return acc;
  }



//



template<typename T>
arma_warn_unused
inline
typename arma_scalar_only<T>::result
accu(const T& x)
  {
  return x;
  }



//! accumulate values in a sparse object
template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const SpBase<typename T1::elem_type,T1>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> P(expr.get_ref());
  
  const uword N = P.get_n_nonzero();
  
  if(N == 0)  { return eT(0); }
  
  if(SpProxy<T1>::use_iterator == false)
    {
    // direct counting
    return arrayops::accumulate(P.get_values(), N);
    }
  
  if(is_SpSubview<typename SpProxy<T1>::stored_type>::value)
    {
    const SpSubview<eT>& sv = reinterpret_cast< const SpSubview<eT>& >(P.Q);
    
    if(sv.n_rows == sv.m.n_rows)
      {
      const SpMat<eT>& m   = sv.m;
      const uword      col = sv.aux_col1;
      
      return arrayops::accumulate(&(m.values[ m.col_ptrs[col] ]), N);
      }
    }
  
  typename SpProxy<T1>::const_iterator_type it = P.begin();
  
  eT val = eT(0);
  
  for(uword i=0; i < N; ++i)  { val += (*it); ++it; }
  
  return val;
  }



//! explicit handling of accu(A + B), where A and B are sparse matrices
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
accu(const SpGlue<T1,T2,spglue_plus>& expr)
  {
  arma_debug_sigprint();
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  arma_conform_assert_same_size(UA.M.n_rows, UA.M.n_cols, UB.M.n_rows, UB.M.n_cols, "addition");
  
  return (accu(UA.M) + accu(UB.M));
  }



//! explicit handling of accu(A - B), where A and B are sparse matrices
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
accu(const SpGlue<T1,T2,spglue_minus>& expr)
  {
  arma_debug_sigprint();
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  arma_conform_assert_same_size(UA.M.n_rows, UA.M.n_cols, UB.M.n_rows, UB.M.n_cols, "subtraction");
  
  return (accu(UA.M) - accu(UB.M));
  }



//! explicit handling of accu(A % B), where A and B are sparse matrices
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
accu(const SpGlue<T1,T2,spglue_schur>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> px(expr.A);
  const SpProxy<T2> py(expr.B);
  
  arma_conform_assert_same_size(px.get_n_rows(), px.get_n_cols(), py.get_n_rows(), py.get_n_cols(), "element-wise multiplication");
  
  if( (px.get_n_nonzero() == 0) && (py.get_n_nonzero() == 0) )  { return eT(0); }
  
  typedef typename SpProxy<T1>::stored_type px_Q_type;
  typedef typename SpProxy<T2>::stored_type py_Q_type;
  
  if(is_SpMat<px_Q_type>::value && is_SpMat<py_Q_type>::value)
    {
    const unwrap_spmat<px_Q_type> UX(px.Q);
    const unwrap_spmat<py_Q_type> UY(py.Q);
    
    const SpMat<eT>& X = UX.M;
    const SpMat<eT>& Y = UY.M;
    
    if(&X == &Y)  { return op_dot::direct_dot(X.n_nonzero, X.values, X.values); }
    }
  
  typename SpProxy<T1>::const_iterator_type x_it     = px.begin();
  typename SpProxy<T1>::const_iterator_type x_it_end = px.end();
  
  typename SpProxy<T2>::const_iterator_type y_it     = py.begin();
  typename SpProxy<T2>::const_iterator_type y_it_end = py.end();
  
  eT acc = eT(0);
  
  while( (x_it != x_it_end) || (y_it != y_it_end) )
    {
    if(x_it == y_it)
      {
      acc += ((*x_it) * (*y_it));
      
      ++x_it;
      ++y_it;
      }
    else
      {
      const uword x_it_col = x_it.col();
      const uword x_it_row = x_it.row();
      
      const uword y_it_col = y_it.col();
      const uword y_it_row = y_it.row();
      
      if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
        {
        acc += (*x_it) * eT(0);  // in case (*x_it) is inf or nan
        
        ++x_it;
        }
      else // x is closer to the end
        {
        acc += eT(0) * (*y_it);  // in case (*y_it) is inf or nan
        
        ++y_it;
        }
      }
    }
  
  return acc;
  }



template<typename T1, typename spop_type>
arma_warn_unused
inline
typename T1::elem_type
accu(const SpOp<T1, spop_type>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  constexpr bool is_vectorise = \
       (is_same_type<spop_type, spop_vectorise_row>::yes)
    || (is_same_type<spop_type, spop_vectorise_col>::yes)
    || (is_same_type<spop_type, spop_vectorise_all>::yes);
  
  if(is_vectorise)  { return accu(expr.m); }
  
  const SpMat<eT> tmp = expr;
  
  return accu(tmp);
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const SpOp<T1, spop_square>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_SpSubview_col<T1>::value)
    {
    const SpSubview_col<eT>& svcol = reinterpret_cast<const SpSubview_col<eT>&>(expr.m);
    
    if(svcol.n_nonzero == 0)  { return eT(0); }
    
    if(svcol.n_rows == svcol.m.n_rows)
      {
      arma_debug_print("accu(): SpSubview_col spop_square optimisation");
      
      const SpMat<eT>& m   = svcol.m;
      const uword      col = svcol.aux_col1;
      
      const eT* ptr = &(m.values[ m.col_ptrs[col] ]);
      
      return op_dot::direct_dot(svcol.n_nonzero, ptr, ptr);
      }
    }
  
  const SpProxy<T1> P(expr.m);
  
  const uword N = P.get_n_nonzero();
  
  if(N == 0)  { return eT(0); }
  
  if(SpProxy<T1>::use_iterator == false)
    {
    return op_dot::direct_dot(N, P.get_values(), P.get_values());
    }
  else
    {
    typename SpProxy<T1>::const_iterator_type it = P.begin();
    
    eT acc = eT(0);
    
    for(uword i=0; i < N; ++i)  { const eT tmp = (*it); acc += (tmp*tmp); ++it; }
    
    return acc;
    }
  }



template<typename T1, typename functor>
inline
typename T1::elem_type
accu_spop_omit_helper(const T1& expr, functor is_omitted)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  constexpr eT eT_zero = eT(0);
  
  if(is_SpSubview_col<T1>::value)
    {
    const SpSubview_col<eT>& svcol = reinterpret_cast<const SpSubview_col<eT>&>(expr);
    
    if(svcol.n_nonzero == 0)  { return eT(0); }
    
    if(svcol.n_rows == svcol.m.n_rows)
      {
      arma_debug_print("accu_spop_omit_helper(): SpSubview_col optimisation");
      
      const SpMat<eT>& m   = svcol.m;
      const uword      col = svcol.aux_col1;
      
      const eT* vals = &(m.values[ m.col_ptrs[col] ]);
      
      const uword N = svcol.n_nonzero;
      
      eT acc = eT(0);
      
      for(uword i=0; i < N; ++i)
        {
        const eT tmp = vals[i];
        
        acc += is_omitted(tmp) ? eT_zero : tmp;
        }
      
      return acc;
      }
    }
  
  const SpProxy<T1> P(expr);
  
  const uword N = P.get_n_nonzero();
  
  if(N == 0)  { return eT(0); }
  
  eT acc = eT(0);
  
  if(SpProxy<T1>::use_iterator == false)
    {
    const eT* vals = P.get_values();
    
    for(uword i=0; i < N; ++i)
      {
      const eT tmp = vals[i];
      
      acc += is_omitted(tmp) ? eT_zero : tmp;
      }
    }
  else
    {
    typename SpProxy<T1>::const_iterator_type it = P.begin();
    
    for(uword i=0; i < N; ++i)
      {
      const eT tmp = (*it);
      
      acc += is_omitted(tmp) ? eT_zero : tmp;
      
      ++it;
      }
    }
  
  return acc;
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
accu(const SpOp<T1, spop_omit>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword omit_mode = expr.aux_uword_a;
  
  if(arma_config::fast_math_warn)
    {
    if(omit_mode == 1)  { arma_warn(1, "omit_nan(): detection of NaN is not reliable in fast math mode"); }
    if(omit_mode == 2)  { arma_warn(1, "omit_nonfinite(): detection of non-finite values is not reliable in fast math mode"); }
    }
  
  auto is_omitted_1 = [](const eT& x) -> bool { return arma_isnan(x);       };
  auto is_omitted_2 = [](const eT& x) -> bool { return arma_isnonfinite(x); };
  
  eT acc = eT(0);
  
  if(omit_mode == 1)  { acc = accu_spop_omit_helper(expr.m, is_omitted_1); }
  if(omit_mode == 2)  { acc = accu_spop_omit_helper(expr.m, is_omitted_2); }
  
  return acc;
  }



template<typename T1, typename spop_type>
arma_warn_unused
inline
uword
accu(const mtSpOp<uword,T1,spop_type>& X, const typename arma_spop_rel_only<spop_type>::result* junk1 = nullptr, const typename arma_not_cx<typename T1::elem_type>::result* junk2 = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const SpProxy<T1> P(X.m);
  
  const uword n_zeros = P.get_n_elem() - P.get_n_nonzero();
  
  const eT zero = eT(0);
  
  // shortcuts
  
  if( (is_same_type<spop_type, spop_rel_eq   >::yes) && (k == zero) )  { return n_zeros;           }
  if( (is_same_type<spop_type, spop_rel_noteq>::yes) && (k == zero) )  { return P.get_n_nonzero(); }
  
  // take into account all implicit zeros
  
  bool use_n_zeros;
  
       if(is_same_type<spop_type, spop_rel_eq       >::yes)  { use_n_zeros = (zero == k   ); }
  else if(is_same_type<spop_type, spop_rel_noteq    >::yes)  { use_n_zeros = (zero != k   ); }
  else if(is_same_type<spop_type, spop_rel_lt_pre   >::yes)  { use_n_zeros = (k    <  zero); }
  else if(is_same_type<spop_type, spop_rel_lt_post  >::yes)  { use_n_zeros = (zero <  k   ); }
  else if(is_same_type<spop_type, spop_rel_gt_pre   >::yes)  { use_n_zeros = (k    >  zero); }
  else if(is_same_type<spop_type, spop_rel_gt_post  >::yes)  { use_n_zeros = (zero >  k   ); }
  else if(is_same_type<spop_type, spop_rel_lteq_pre >::yes)  { use_n_zeros = (k    <= zero); }
  else if(is_same_type<spop_type, spop_rel_lteq_post>::yes)  { use_n_zeros = (zero <= k   ); }
  else if(is_same_type<spop_type, spop_rel_gteq_pre >::yes)  { use_n_zeros = (k    >= zero); }
  else if(is_same_type<spop_type, spop_rel_gteq_post>::yes)  { use_n_zeros = (zero >= k   ); }
  else { use_n_zeros = false; }
  
  uword count = (use_n_zeros) ? n_zeros : 0;
  
  typename SpProxy<T1>::const_iterator_type it     = P.begin();
  typename SpProxy<T1>::const_iterator_type it_end = P.end();
  
  // take into account all non-zero elements
  
  for(; it != it_end; ++it)
    {
    const eT val = (*it);
    
    bool condition;
    
         if(is_same_type<spop_type, spop_rel_eq       >::yes)  { condition = (val == k  ); }
    else if(is_same_type<spop_type, spop_rel_noteq    >::yes)  { condition = (val != k  ); }
    else if(is_same_type<spop_type, spop_rel_lt_pre   >::yes)  { condition = (k   <  val); }
    else if(is_same_type<spop_type, spop_rel_lt_post  >::yes)  { condition = (val <  k  ); }
    else if(is_same_type<spop_type, spop_rel_gt_pre   >::yes)  { condition = (k   >  val); }
    else if(is_same_type<spop_type, spop_rel_gt_post  >::yes)  { condition = (val >  k  ); }
    else if(is_same_type<spop_type, spop_rel_lteq_pre >::yes)  { condition = (k   <= val); }
    else if(is_same_type<spop_type, spop_rel_lteq_post>::yes)  { condition = (val <= k  ); }
    else if(is_same_type<spop_type, spop_rel_gteq_pre >::yes)  { condition = (k   >= val); }
    else if(is_same_type<spop_type, spop_rel_gteq_post>::yes)  { condition = (val >= k  ); }
    else { condition = false; }
    
    count += (condition) ? uword(1) : uword(0);
    }
  
  return count;
  }



template<typename T1, typename spop_type>
arma_warn_unused
inline
uword
accu(const mtSpOp<uword,T1,spop_type>& X, const typename arma_spop_rel_only<spop_type>::result* junk1 = nullptr, const typename arma_cx_only<typename T1::elem_type>::result* junk2 = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const SpProxy<T1> P(X.m);
  
  const uword n_zeros = P.get_n_elem() - P.get_n_nonzero();
  
  const eT zero = eT(0);
  
  // shortcuts
  
  if( (is_same_type<spop_type, spop_rel_eq   >::yes) && (k == zero) )  { return n_zeros;           }
  if( (is_same_type<spop_type, spop_rel_noteq>::yes) && (k == zero) )  { return P.get_n_nonzero(); }
  
  // take into account all implicit zeros
  
  bool use_n_zeros;
  
       if(is_same_type<spop_type, spop_rel_eq   >::yes)  { use_n_zeros = (zero == k); }
  else if(is_same_type<spop_type, spop_rel_noteq>::yes)  { use_n_zeros = (zero != k); }
  else { use_n_zeros = false; }
  
  uword count = (use_n_zeros) ? n_zeros : 0;
  
  typename SpProxy<T1>::const_iterator_type it     = P.begin();
  typename SpProxy<T1>::const_iterator_type it_end = P.end();
  
  // take into account all non-zero elements
  
  for(; it != it_end; ++it)
    {
    const eT val = (*it);
    
    bool condition;
    
         if(is_same_type<spop_type, spop_rel_eq   >::yes)  { condition = (val == k); }
    else if(is_same_type<spop_type, spop_rel_noteq>::yes)  { condition = (val != k); }
    else { condition = false; }
    
    count += (condition) ? uword(1) : uword(0);
    }
  
  return count;
  }



//! @}
