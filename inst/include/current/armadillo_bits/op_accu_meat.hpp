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


//! \addtogroup op_accu
//! @{



template<typename T1>
inline
typename T1::elem_type
op_accu_mat::apply_proxy_linear(const Proxy<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  eT val = eT(0);
  
  typename Proxy<T1>::ea_type Pea = P.get_ea();
  
  const uword n_elem = P.get_n_elem();
  
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
  
  return val;
  }



template<typename T1>
inline
typename T1::elem_type
op_accu_mat::apply_proxy_at(const Proxy<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
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



template<typename T1>
inline
typename T1::elem_type
op_accu_mat::apply(const T1& X)
  {
  arma_debug_sigprint();
  
  if( (quasi_unwrap<T1>::has_orig_mem) || (is_Mat<typename Proxy<T1>::stored_type>::value) || (arma_config::openmp && Proxy<T1>::use_mp) )
    {
    const quasi_unwrap<T1> U(X);
    
    return arrayops::accumulate(U.M.memptr(), U.M.n_elem);
    }
  
  const Proxy<T1> P(X);
  
  return (Proxy<T1>::use_at) ? op_accu_mat::apply_proxy_at(P) : op_accu_mat::apply_proxy_linear(P);
  }



template<typename T1, typename functor>
inline
typename T1::elem_type
op_accu_mat::apply_omit_helper(const Proxy<T1>& P, functor is_omitted)
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
inline
typename T1::elem_type
op_accu_mat::apply(const Op<T1, op_omit>& in)
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
  
  if(omit_mode == 1)  { acc = op_accu_mat::apply_omit_helper(P, is_omitted_1); }
  if(omit_mode == 2)  { acc = op_accu_mat::apply_omit_helper(P, is_omitted_2); }
  
  return acc;
  }



template<typename T1>
inline
typename T1::elem_type
op_accu_mat::apply(const eOp<T1,eop_square>& expr)
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
  
  return (Proxy<expr_type>::use_at) ? op_accu_mat::apply_proxy_at(P) : op_accu_mat::apply_proxy_linear(P);
  }



template<typename T1>
inline
typename T1::elem_type
op_accu_mat::apply(const eOp<T1,eop_pow>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef eOp<T1,eop_pow> expr_type;
  
  if(arma_config::optimise_powexpr && (expr.aux == eT(2)))
    {
    typedef eOp<T1,eop_square> modified_expr_type;
    
    return op_accu_mat::apply( reinterpret_cast< const modified_expr_type& >(expr) );
    }
  
  if(arma_config::optimise_powexpr && (expr.aux == eT(0.5)) && is_real_or_cx<eT>::value)
    {
    typedef eOp<T1,eop_sqrt> modified_expr_type;
    
    return op_accu_mat::apply( reinterpret_cast< const modified_expr_type& >(expr) );
    }
  
  const Proxy<expr_type> P(expr);
  
  return (Proxy<expr_type>::use_at) ? op_accu_mat::apply_proxy_at(P) : op_accu_mat::apply_proxy_linear(P);
  }



template<typename T1, typename T2>
inline
typename T1::elem_type
op_accu_mat::apply(const eGlue<T1,T2,eglue_schur>& expr)
  {
  arma_debug_sigprint();
  
  typedef eGlue<T1,T2,eglue_schur> expr_type;
  
  typedef typename expr_type::proxy1_type::stored_type P1_stored_type;
  typedef typename expr_type::proxy2_type::stored_type P2_stored_type;
  
  constexpr bool is_sv = (is_subview<P1_stored_type>::value) || (is_subview<P2_stored_type>::value);
  
  if( (is_sv) && (expr.get_n_rows() >= 4) )
    {
    arma_debug_print("op_accu_mat::apply(): eglue_schur subview optimisation");
    
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
    arma_debug_print("op_accu_mat::apply(): eglue_schur direct_mem optimisation");
    
    const quasi_unwrap<P1_stored_type> tmp1(expr.P1.Q);
    const quasi_unwrap<P2_stored_type> tmp2(expr.P2.Q);
    
    return op_dot::direct_dot(tmp1.M.n_elem, tmp1.M.memptr(), tmp2.M.memptr());
    }
  
  const Proxy<expr_type> P(expr);
  
  return (Proxy<expr_type>::use_at) ? op_accu_mat::apply_proxy_at(P) : op_accu_mat::apply_proxy_linear(P);
  }



template<typename T1, typename op_type>
inline
uword
op_accu_mat::apply(const mtOp<uword,T1,op_type>& X, const typename arma_op_rel_only<op_type>::result* junk1, const typename arma_not_cx<typename T1::elem_type>::result* junk2)
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
inline
uword
op_accu_mat::apply(const mtOp<uword,T1,op_type>& X, const typename arma_op_rel_only<op_type>::result* junk1, const typename arma_cx_only<typename T1::elem_type>::result* junk2)
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
inline
uword
op_accu_mat::apply(const mtGlue<uword,T1,T2,glue_rel_noteq>& X)
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
inline
uword
op_accu_mat::apply(const mtGlue<uword,T1,T2,glue_rel_eq>& X)
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



template<typename eT>
inline
eT
op_accu_mat::apply(const subview<eT>& X)
  {
  arma_debug_sigprint();  
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(X_n_rows == 1)
    {
    const uword X_m_n_rows = X.m.n_rows;
    
    const eT* mem_ptr = X.colptr(0);
    
    eT val1 = eT(0);
    eT val2 = eT(0);
    
    uword j;
    
    for(j=1; j < X_n_cols; j+=2)
      {
      val1 += (*mem_ptr); mem_ptr += X_m_n_rows;
      val2 += (*mem_ptr); mem_ptr += X_m_n_rows;
      }
    
    if((j-1) < X_n_cols)
      {
      val1 += (*mem_ptr);
      }
    
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
inline
eT
op_accu_mat::apply(const subview_col<eT>& X)
  {
  arma_debug_sigprint();  
  
  return arrayops::accumulate( X.colmem, X.n_rows );
  }



template<typename eT>
inline
eT
op_accu_mat::apply(const subview_row<eT>& X)
  {
  arma_debug_sigprint();  
  
  const uword X_m_n_rows = X.m.n_rows;
  const uword X_n_cols   = X.n_cols;
  
  const eT* mem_ptr = X.rowmem;
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  uword j;
  
  for(j=1; j < X_n_cols; j+=2)
    {
    val1 += (*mem_ptr); mem_ptr += X_m_n_rows;
    val2 += (*mem_ptr); mem_ptr += X_m_n_rows;
    }
  
  if((j-1) < X_n_cols)
    {
    val1 += (*mem_ptr);
    }
  
  return val1 + val2;
  }



//



template<typename T1>
inline
typename T1::elem_type
op_accu_cube::apply_proxy_linear(const ProxyCube<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  eT val = eT(0);
  
  typename ProxyCube<T1>::ea_type Pea = P.get_ea();
  
  const uword n_elem = P.get_n_elem();
  
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
    
    if(i < n_elem)  { val1 += Pea[i]; }
    
    val = val1 + val2;
    }
  #endif
  
  return val;
  }



template<typename T1>
inline
typename T1::elem_type
op_accu_cube::apply_proxy_at(const ProxyCube<T1>& P)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
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



template<typename T1>
inline
typename T1::elem_type
op_accu_cube::apply(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_debug_sigprint();
  
  if( (is_Cube<T1>::value) || (is_Cube<typename ProxyCube<T1>::stored_type>::value) || (arma_config::openmp && ProxyCube<T1>::use_mp) )
    {
    const unwrap_cube<T1> U(X.get_ref());
    
    return arrayops::accumulate(U.M.memptr(), U.M.n_elem);
    }
  
  const ProxyCube<T1> P(X.get_ref());
  
  return (ProxyCube<T1>::use_at) ? op_accu_cube::apply_proxy_at(P) : op_accu_cube::apply_proxy_linear(P);
  }



template<typename T1>
inline
typename T1::elem_type
op_accu_cube::apply(const eOpCube<T1,eop_square>& expr)
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
  
  return (ProxyCube<expr_type>::use_at) ? op_accu_cube::apply_proxy_at(P) : op_accu_cube::apply_proxy_linear(P);
  }



template<typename T1>
inline
typename T1::elem_type
op_accu_cube::apply(const eOpCube<T1,eop_pow>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef eOpCube<T1,eop_pow> expr_type;
  
  if(arma_config::optimise_powexpr && (expr.aux == eT(2)))
    {
    typedef eOpCube<T1,eop_square> modified_expr_type;
    
    return op_accu_cube::apply( reinterpret_cast< const modified_expr_type& >(expr) );
    }
  
  if(arma_config::optimise_powexpr && (expr.aux == eT(0.5)) && is_real_or_cx<eT>::value)
    {
    typedef eOpCube<T1,eop_sqrt> modified_expr_type;
    
    return op_accu_cube::apply( reinterpret_cast< const modified_expr_type& >(expr) );
    }
  
  const ProxyCube<expr_type> P(expr);
  
  return (ProxyCube<expr_type>::use_at) ? op_accu_cube::apply_proxy_at(P) : op_accu_cube::apply_proxy_linear(P);
  }



template<typename T1, typename T2>
inline
typename T1::elem_type
op_accu_cube::apply(const eGlueCube<T1,T2,eglue_schur>& expr)
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
  
  return (ProxyCube<expr_type>::use_at) ? op_accu_cube::apply_proxy_at(P) : op_accu_cube::apply_proxy_linear(P);
  }



template<typename T1, typename functor>
inline
typename T1::elem_type
op_accu_cube::apply_omit_helper(const ProxyCube<T1>& P, functor is_omitted)
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
inline
typename T1::elem_type
op_accu_cube::apply(const CubeToMatOp<T1, op_omit_cube>& in)
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
  
  if(omit_mode == 1)  { acc = op_accu_cube::apply_omit_helper(P, is_omitted_1); }
  if(omit_mode == 2)  { acc = op_accu_cube::apply_omit_helper(P, is_omitted_2); }
  
  return acc;
  }



template<typename eT>
inline
eT
op_accu_cube::apply(const subview_cube<eT>& sv)
  {
  arma_debug_sigprint();  
  
  if(sv.n_elem == 0)  { return eT(0); }
  
  const uword sv_nr = sv.n_rows;
  const uword sv_nc = sv.n_cols;
  const uword sv_ns = sv.n_slices;
  
  eT acc = eT(0);
  
  if( (sv_nr == 1)  && (sv_nc == 1) && (sv.aux_slice1 == 0) )
    {
    const uword sv_m_n_elem_slice = sv.m.n_elem_slice;
    
    const eT* sv_m_ptr = &( sv.m.at(sv.aux_row1, sv.aux_col1, 0) );
    
    for(uword s=0; s < sv_ns; ++s)
      {
      acc += (*sv_m_ptr);  sv_m_ptr += sv_m_n_elem_slice;
      }
    }
  else
    {
    for(uword s=0; s < sv_ns; ++s)
    for(uword c=0; c < sv_nc; ++c)
      {
      acc += arrayops::accumulate(sv.slice_colptr(s,c), sv_nr);
      }
    }
  
  return acc;
  }



//! @}
