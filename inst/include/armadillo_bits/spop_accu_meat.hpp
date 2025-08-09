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


//! \addtogroup spop_accu
//! @{



template<typename T1>
inline
typename T1::elem_type
spop_accu::apply(const SpBase<typename T1::elem_type,T1>& expr)
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



template<typename T1, typename T2>
inline
typename T1::elem_type
spop_accu::apply(const SpGlue<T1,T2,spglue_plus>& expr)
  {
  arma_debug_sigprint();
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  arma_conform_assert_same_size(UA.M.n_rows, UA.M.n_cols, UB.M.n_rows, UB.M.n_cols, "addition");
  
  return (spop_accu::apply(UA.M) + spop_accu::apply(UB.M));
  }



template<typename T1, typename T2>
inline
typename T1::elem_type
spop_accu::apply(const SpGlue<T1,T2,spglue_minus>& expr)
  {
  arma_debug_sigprint();
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  arma_conform_assert_same_size(UA.M.n_rows, UA.M.n_cols, UB.M.n_rows, UB.M.n_cols, "subtraction");
  
  return (spop_accu::apply(UA.M) - spop_accu::apply(UB.M));
  }



template<typename T1, typename T2>
inline
typename T1::elem_type
spop_accu::apply(const SpGlue<T1,T2,spglue_schur>& expr)
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
inline
typename T1::elem_type
spop_accu::apply(const SpOp<T1, spop_type>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  constexpr bool is_vectorise = \
       (is_same_type<spop_type, spop_vectorise_row>::yes)
    || (is_same_type<spop_type, spop_vectorise_col>::yes)
    || (is_same_type<spop_type, spop_vectorise_all>::yes);
  
  if(is_vectorise)  { return spop_accu::apply(expr.m); }
  
  const SpMat<eT> tmp = expr;
  
  return spop_accu::apply(tmp);
  }



template<typename T1>
inline
typename T1::elem_type
spop_accu::apply(const SpOp<T1, spop_square>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_SpSubview_col<T1>::value)
    {
    const SpSubview_col<eT>& svcol = reinterpret_cast<const SpSubview_col<eT>&>(expr.m);
    
    if(svcol.n_nonzero == 0)  { return eT(0); }
    
    if(svcol.n_rows == svcol.m.n_rows)
      {
      arma_debug_print("spop_accu::apply(): SpSubview_col spop_square optimisation");
      
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
spop_accu::apply_spop_omit_helper(const T1& expr, functor is_omitted)
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
      arma_debug_print("spop_accu::apply_spop_omit_helper(): SpSubview_col optimisation");
      
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
inline
typename T1::elem_type
spop_accu::apply(const SpOp<T1, spop_omit>& expr)
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
  
  if(omit_mode == 1)  { acc = spop_accu::apply_spop_omit_helper(expr.m, is_omitted_1); }
  if(omit_mode == 2)  { acc = spop_accu::apply_spop_omit_helper(expr.m, is_omitted_2); }
  
  return acc;
  }



template<typename T1, typename spop_type>
inline
uword
spop_accu::apply(const mtSpOp<uword,T1,spop_type>& X, const typename arma_spop_rel_only<spop_type>::result* junk1, const typename arma_not_cx<typename T1::elem_type>::result* junk2)
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
inline
uword
spop_accu::apply(const mtSpOp<uword,T1,spop_type>& X, const typename arma_spop_rel_only<spop_type>::result* junk1, const typename arma_cx_only<typename T1::elem_type>::result* junk2)
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
