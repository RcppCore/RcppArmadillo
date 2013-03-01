// Copyright (C) 2010-2013 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2013 Conrad Sanderson
// Copyright (C) 2010 Dimitrios Bouzas
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_find
//! @{



template<typename T1>
inline
uword
op_find::helper
  (
  Mat<uword>& indices,
  const Base<typename T1::elem_type, T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.get_ref());
  
  const uword n_elem = A.get_n_elem();
  
  indices.set_size(n_elem, 1);
  
  uword* indices_mem = indices.memptr();
  uword  n_nz        = 0;
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type PA = A.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      if(PA[i] != eT(0))  { indices_mem[n_nz] = i;  ++n_nz; }
      }
    }
  else
    {
    const uword n_rows = A.get_n_rows();
    const uword n_cols = A.get_n_cols();
    
    uword i = 0;
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      if(A.at(row,col) != eT(0))  { indices_mem[n_nz] = i; ++n_nz; }
      
      ++i;
      }
    }
  
  return n_nz;
  }



template<typename T1, typename op_type>
inline
uword
op_find::helper
  (
  Mat<uword>& indices,
  const mtOp<uword, T1, op_type>& X,
  const typename arma_op_rel_only<op_type>::result           junk1,
  const typename arma_not_cx<typename T1::elem_type>::result junk2
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::elem_type eT;
  
  const eT val = X.aux;
  
  const Proxy<T1> A(X.m);
  
  const uword n_elem = A.get_n_elem();
  
  indices.set_size(n_elem, 1);
  
  uword* indices_mem = indices.memptr();
  uword  n_nz        = 0;
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type PA = A.get_ea();
    
    uword i,j;
    for(i=0, j=1; j < n_elem; i+=2, j+=2)
      {
      const eT tpi = PA[i];
      const eT tpj = PA[j];
      
      bool not_zero_i;
      bool not_zero_j;
      
           if(is_same_type<op_type, op_rel_lt_pre   >::value == true)  { not_zero_i = (val <  tpi); }
      else if(is_same_type<op_type, op_rel_lt_post  >::value == true)  { not_zero_i = (tpi <  val); }
      else if(is_same_type<op_type, op_rel_gt_pre   >::value == true)  { not_zero_i = (val >  tpi); }
      else if(is_same_type<op_type, op_rel_gt_post  >::value == true)  { not_zero_i = (tpi >  val); }
      else if(is_same_type<op_type, op_rel_lteq_pre >::value == true)  { not_zero_i = (val <= tpi); }
      else if(is_same_type<op_type, op_rel_lteq_post>::value == true)  { not_zero_i = (tpi <= val); }
      else if(is_same_type<op_type, op_rel_gteq_pre >::value == true)  { not_zero_i = (val >= tpi); }
      else if(is_same_type<op_type, op_rel_gteq_post>::value == true)  { not_zero_i = (tpi >= val); }
      else if(is_same_type<op_type, op_rel_eq       >::value == true)  { not_zero_i = (tpi == val); }
      else if(is_same_type<op_type, op_rel_noteq    >::value == true)  { not_zero_i = (tpi != val); }
      else not_zero_i = false;
      
           if(is_same_type<op_type, op_rel_lt_pre   >::value == true)  { not_zero_j = (val <  tpj); }
      else if(is_same_type<op_type, op_rel_lt_post  >::value == true)  { not_zero_j = (tpj <  val); }
      else if(is_same_type<op_type, op_rel_gt_pre   >::value == true)  { not_zero_j = (val >  tpj); }
      else if(is_same_type<op_type, op_rel_gt_post  >::value == true)  { not_zero_j = (tpj >  val); }
      else if(is_same_type<op_type, op_rel_lteq_pre >::value == true)  { not_zero_j = (val <= tpj); }
      else if(is_same_type<op_type, op_rel_lteq_post>::value == true)  { not_zero_j = (tpj <= val); }
      else if(is_same_type<op_type, op_rel_gteq_pre >::value == true)  { not_zero_j = (val >= tpj); }
      else if(is_same_type<op_type, op_rel_gteq_post>::value == true)  { not_zero_j = (tpj >= val); }
      else if(is_same_type<op_type, op_rel_eq       >::value == true)  { not_zero_j = (tpj == val); }
      else if(is_same_type<op_type, op_rel_noteq    >::value == true)  { not_zero_j = (tpj != val); }
      else not_zero_j = false;
      
      if(not_zero_i == true)  { indices_mem[n_nz] = i;  ++n_nz; }
      if(not_zero_j == true)  { indices_mem[n_nz] = j;  ++n_nz; }
      }
    
    if(i < n_elem)
      {
      bool not_zero;
      
      const eT tmp = PA[i];
      
           if(is_same_type<op_type, op_rel_lt_pre   >::value == true)  { not_zero = (val <  tmp); }
      else if(is_same_type<op_type, op_rel_lt_post  >::value == true)  { not_zero = (tmp <  val); }
      else if(is_same_type<op_type, op_rel_gt_pre   >::value == true)  { not_zero = (val >  tmp); }
      else if(is_same_type<op_type, op_rel_gt_post  >::value == true)  { not_zero = (tmp >  val); }
      else if(is_same_type<op_type, op_rel_lteq_pre >::value == true)  { not_zero = (val <= tmp); }
      else if(is_same_type<op_type, op_rel_lteq_post>::value == true)  { not_zero = (tmp <= val); }
      else if(is_same_type<op_type, op_rel_gteq_pre >::value == true)  { not_zero = (val >= tmp); }
      else if(is_same_type<op_type, op_rel_gteq_post>::value == true)  { not_zero = (tmp >= val); }
      else if(is_same_type<op_type, op_rel_eq       >::value == true)  { not_zero = (tmp == val); }
      else if(is_same_type<op_type, op_rel_noteq    >::value == true)  { not_zero = (tmp != val); }
      else not_zero = false;
      
      if(not_zero == true)  { indices_mem[n_nz] = i;  ++n_nz; }
      }
    }
  else
    {
    const uword n_rows = A.get_n_rows();
    const uword n_cols = A.get_n_cols();
    
    uword i = 0;
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      const eT tmp = A.at(row,col);
      
      bool not_zero;
      
           if(is_same_type<op_type, op_rel_lt_pre   >::value == true)  { not_zero = (val <  tmp); }
      else if(is_same_type<op_type, op_rel_lt_post  >::value == true)  { not_zero = (tmp <  val); }
      else if(is_same_type<op_type, op_rel_gt_pre   >::value == true)  { not_zero = (val >  tmp); }
      else if(is_same_type<op_type, op_rel_gt_post  >::value == true)  { not_zero = (tmp >  val); }
      else if(is_same_type<op_type, op_rel_lteq_pre >::value == true)  { not_zero = (val <= tmp); }
      else if(is_same_type<op_type, op_rel_lteq_post>::value == true)  { not_zero = (tmp <= val); }
      else if(is_same_type<op_type, op_rel_gteq_pre >::value == true)  { not_zero = (val >= tmp); }
      else if(is_same_type<op_type, op_rel_gteq_post>::value == true)  { not_zero = (tmp >= val); }
      else if(is_same_type<op_type, op_rel_eq       >::value == true)  { not_zero = (tmp == val); }
      else if(is_same_type<op_type, op_rel_noteq    >::value == true)  { not_zero = (tmp != val); }
      else not_zero = false;
      
      if(not_zero == true)  { indices_mem[n_nz] = i;  ++n_nz; }
      
      ++i;
      }
    }
  
  return n_nz;
  }



template<typename T1, typename op_type>
inline
uword
op_find::helper
  (
  Mat<uword>& indices,
  const mtOp<uword, T1, op_type>& X,
  const typename arma_op_rel_only<op_type>::result            junk1,
  const typename arma_cx_only<typename T1::elem_type>::result junk2
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type;
  
  const eT val = X.aux;
  
  const Proxy<T1> A(X.m);
  
  ea_type     PA     = A.get_ea();
  const uword n_elem = A.get_n_elem();
  
  indices.set_size(n_elem, 1);
  
  uword* indices_mem = indices.memptr();
  uword  n_nz        = 0;
  
  for(uword i=0; i<n_elem; ++i)
    {
    const eT tmp = PA[i];
    
    bool not_zero;
    
         if(is_same_type<op_type, op_rel_eq   >::value == true)  { not_zero = (tmp == val); }
    else if(is_same_type<op_type, op_rel_noteq>::value == true)  { not_zero = (tmp != val); }
    else not_zero = false;
    
    if(not_zero == true) { indices_mem[n_nz] = i;  ++n_nz; }
    }
  
  return n_nz;
  }



template<typename T1, typename T2, typename glue_type>
inline
uword
op_find::helper
  (
  Mat<uword>& indices,
  const mtGlue<uword, T1, T2, glue_type>& X,
  const typename arma_glue_rel_only<glue_type>::result       junk1,
  const typename arma_not_cx<typename T1::elem_type>::result junk2,
  const typename arma_not_cx<typename T2::elem_type>::result junk3
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  arma_ignore(junk3);
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "relational operator");
  
  ea_type1 PA = A.get_ea();
  ea_type2 PB = B.get_ea();
  
  const uword n_elem = B.get_n_elem();
  
  indices.set_size(n_elem, 1);
  
  uword* indices_mem = indices.memptr();
  uword  n_nz        = 0;
  
  for(uword i=0; i<n_elem; ++i)
    {
    const eT1 tmp1 = PA[i];
    const eT2 tmp2 = PB[i];
    
    bool not_zero;
    
         if(is_same_type<glue_type, glue_rel_lt    >::value == true)  { not_zero = (tmp1 <  tmp2); }
    else if(is_same_type<glue_type, glue_rel_gt    >::value == true)  { not_zero = (tmp1 >  tmp2); }
    else if(is_same_type<glue_type, glue_rel_lteq  >::value == true)  { not_zero = (tmp1 <= tmp2); }
    else if(is_same_type<glue_type, glue_rel_gteq  >::value == true)  { not_zero = (tmp1 >= tmp2); }
    else if(is_same_type<glue_type, glue_rel_eq    >::value == true)  { not_zero = (tmp1 == tmp2); }
    else if(is_same_type<glue_type, glue_rel_noteq >::value == true)  { not_zero = (tmp1 != tmp2); }
    else not_zero = false;
    
    if(not_zero == true)  { indices_mem[n_nz] = i;  ++n_nz; }
    }
  
  return n_nz;
  }



template<typename T1, typename T2, typename glue_type>
inline
uword
op_find::helper
  (
  Mat<uword>& indices,
  const mtGlue<uword, T1, T2, glue_type>& X,
  const typename arma_glue_rel_only<glue_type>::result        junk1,
  const typename arma_cx_only<typename T1::elem_type>::result junk2,
  const typename arma_cx_only<typename T2::elem_type>::result junk3
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  arma_ignore(junk3);
  
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "relational operator");
  
  ea_type1 PA = A.get_ea();
  ea_type2 PB = B.get_ea();
  
  const uword n_elem = B.get_n_elem();
  
  indices.set_size(n_elem, 1);
  
  uword* indices_mem = indices.memptr();
  uword  n_nz        = 0;
  
  for(uword i=0; i<n_elem; ++i)
    {
    bool not_zero;
    
         if(is_same_type<glue_type, glue_rel_eq    >::value == true)  { not_zero = (PA[i] == PB[i]); }
    else if(is_same_type<glue_type, glue_rel_noteq >::value == true)  { not_zero = (PA[i] != PB[i]); }
    else not_zero = false;
    
    if(not_zero == true)  { indices_mem[n_nz] = i;  ++n_nz; }
    }
  
  return n_nz;
  }



template<typename T1>
inline
void
op_find::apply(Mat<uword>& out, const mtOp<uword, T1, op_find>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword k    = X.aux_uword_a;
  const uword type = X.aux_uword_b;
  
  Mat<uword> indices;
  const uword n_nz = op_find::helper(indices, X.m);
  
  if(n_nz > 0)
    {
    if(type == 0)   // "first"
      {
      out = (k > 0 && k <= n_nz) ? indices.rows(0,      k-1   ) : indices.rows(0, n_nz-1);
      }
    else   // "last"
      {
      out = (k > 0 && k <= n_nz) ? indices.rows(n_nz-k, n_nz-1) : indices.rows(0, n_nz-1);
      }
    }
  else
    {
    out.set_size(0,1);  // empty column vector
    }
  }



//! @}
