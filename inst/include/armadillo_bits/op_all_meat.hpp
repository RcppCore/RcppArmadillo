// Copyright (C) 2013 Conrad Sanderson
// Copyright (C) 2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_all
//! @{



template<typename T1>
inline
bool
op_all::all_vec_helper(const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X.get_ref());
  
  const uword n_elem = P.get_n_elem();
  
  uword count = 0;
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      if(Pea[i] != eT(0))  { ++count; }
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      if(P.at(row,col) != eT(0))  { ++count; }
      }
    }
  
  // NOTE: for empty vectors it makes more sense to return false, but we need to return true for compatability with Octave
  return (n_elem == count);
  }



template<typename T1, typename op_type>
inline
bool
op_all::all_vec_helper
  (
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
  
  const Proxy<T1> P(X.m);
  
  const uword n_elem = P.get_n_elem();
  
  uword count = 0;
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
    
    for(uword i=0; i < n_elem; ++i)
      {
      const eT tmp = Pea[i];
      
           if(is_same_type<op_type, op_rel_lt_pre   >::yes)  { if(val <  tmp) { ++count; } }
      else if(is_same_type<op_type, op_rel_lt_post  >::yes)  { if(tmp <  val) { ++count; } }
      else if(is_same_type<op_type, op_rel_gt_pre   >::yes)  { if(val >  tmp) { ++count; } }
      else if(is_same_type<op_type, op_rel_gt_post  >::yes)  { if(tmp >  val) { ++count; } }
      else if(is_same_type<op_type, op_rel_lteq_pre >::yes)  { if(val <= tmp) { ++count; } }
      else if(is_same_type<op_type, op_rel_lteq_post>::yes)  { if(tmp <= val) { ++count; } }
      else if(is_same_type<op_type, op_rel_gteq_pre >::yes)  { if(val >= tmp) { ++count; } }
      else if(is_same_type<op_type, op_rel_gteq_post>::yes)  { if(tmp >= val) { ++count; } }
      else if(is_same_type<op_type, op_rel_eq       >::yes)  { if(tmp == val) { ++count; } }
      else if(is_same_type<op_type, op_rel_noteq    >::yes)  { if(tmp != val) { ++count; } }
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      const eT tmp = P.at(row,col);
      
           if(is_same_type<op_type, op_rel_lt_pre   >::yes)  { if(val <  tmp) { ++count; } }
      else if(is_same_type<op_type, op_rel_lt_post  >::yes)  { if(tmp <  val) { ++count; } }
      else if(is_same_type<op_type, op_rel_gt_pre   >::yes)  { if(val >  tmp) { ++count; } }
      else if(is_same_type<op_type, op_rel_gt_post  >::yes)  { if(tmp >  val) { ++count; } }
      else if(is_same_type<op_type, op_rel_lteq_pre >::yes)  { if(val <= tmp) { ++count; } }
      else if(is_same_type<op_type, op_rel_lteq_post>::yes)  { if(tmp <= val) { ++count; } }
      else if(is_same_type<op_type, op_rel_gteq_pre >::yes)  { if(val >= tmp) { ++count; } }
      else if(is_same_type<op_type, op_rel_gteq_post>::yes)  { if(tmp >= val) { ++count; } }
      else if(is_same_type<op_type, op_rel_eq       >::yes)  { if(tmp == val) { ++count; } }
      else if(is_same_type<op_type, op_rel_noteq    >::yes)  { if(tmp != val) { ++count; } }
      }
    }
  
  return (n_elem == count);
  }



template<typename T1, typename T2, typename glue_type>
inline
bool
op_all::all_vec_helper
  (
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
  
  const uword n_elem = A.get_n_elem();
  
  uword count = 0;
  
  const bool prefer_at_accessor = Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor;
  
  if(prefer_at_accessor == false)
    {
    ea_type1 PA = A.get_ea();
    ea_type2 PB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      const eT1 tmp1 = PA[i];
      const eT2 tmp2 = PB[i];
      
           if(is_same_type<glue_type, glue_rel_lt    >::yes)  { if(tmp1 <  tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_gt    >::yes)  { if(tmp1 >  tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_lteq  >::yes)  { if(tmp1 <= tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_gteq  >::yes)  { if(tmp1 >= tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_eq    >::yes)  { if(tmp1 == tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_noteq >::yes)  { if(tmp1 != tmp2) { ++count; } }
      }
    }
  else
    {
    const uword n_rows = A.get_n_rows();
    const uword n_cols = A.get_n_cols();
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      const eT1 tmp1 = A.at(row,col);
      const eT2 tmp2 = B.at(row,col);
      
           if(is_same_type<glue_type, glue_rel_lt    >::yes)  { if(tmp1 <  tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_gt    >::yes)  { if(tmp1 >  tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_lteq  >::yes)  { if(tmp1 <= tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_gteq  >::yes)  { if(tmp1 >= tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_eq    >::yes)  { if(tmp1 == tmp2) { ++count; } }
      else if(is_same_type<glue_type, glue_rel_noteq >::yes)  { if(tmp1 != tmp2) { ++count; } }
      }
    }
  
  return (n_elem == count);
  }



template<typename T1>
inline
bool
op_all::all_vec(T1& X)
  {
  arma_extra_debug_sigprint();
  
  return op_all::all_vec_helper(X);
  }



template<typename T1>
inline
void
op_all::apply_helper(Mat<uword>& out, const Proxy<T1>& P, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  typedef typename Proxy<T1>::elem_type eT;
  
  if(dim == 0)  // traverse rows (ie. process each column)
    {
    out.zeros(1, n_cols);
    
    if(out.n_elem == 0)  { return; }
    
    uword* out_mem = out.memptr();
    
    if(is_Mat<typename Proxy<T1>::stored_type>::value == true)
      {
      const unwrap<typename Proxy<T1>::stored_type> U(P.Q);
      
      for(uword col=0; col < n_cols; ++col)
        {
        const eT* colmem = U.M.colptr(col);
        
        uword count = 0;
        
        for(uword row=0; row < n_rows; ++row)
          {
          if(colmem[row] != eT(0))  { ++count; }
          }
        
        out_mem[col] = (n_rows == count) ? uword(1) : uword(0); 
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
        {
        uword count = 0;
        
        for(uword row=0; row < n_rows; ++row)
          {
          if(P.at(row,col) != eT(0))  { ++count; }
          }
        
        out_mem[col] = (n_rows == count) ? uword(1) : uword(0); 
        }
      }
    }
  else
    {
    out.zeros(n_rows, 1);
    
    uword* out_mem = out.memptr();
    
    // internal dual use of 'out': keep the counts for each row
    
    if(is_Mat<typename Proxy<T1>::stored_type>::value == true)
      {
      const unwrap<typename Proxy<T1>::stored_type> U(P.Q);
      
      for(uword col=0; col < n_cols; ++col)
        {
        const eT* colmem = U.M.colptr(col);
        
        for(uword row=0; row < n_rows; ++row)
          {
          if(colmem[row] != eT(0))  { ++out_mem[row]; }
          }
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
        {
        for(uword row=0; row < n_rows; ++row)
          {
          if(P.at(row,col) != eT(0))  { ++out_mem[row]; }
          }
        }
      }
    
    
    // see what the counts tell us
    
    for(uword row=0; row < n_rows; ++row)
      {
      out_mem[row] = (n_cols == out_mem[row]) ? uword(1) : uword(0);
      }
    
    }
  }



template<typename T1>
inline
void
op_all::apply(Mat<uword>& out, const mtOp<uword, T1, op_all>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword dim = X.aux_uword_a;
  
  const Proxy<T1> P(X.m);
  
  if(P.is_alias(out) == false)
    {
    op_all::apply_helper(out, P, dim);
    }
  else
    {
    Mat<uword> out2;
    
    op_all::apply_helper(out2, P, dim);
    
    out.steal_mem(out2);
    }
  }



//! @}
