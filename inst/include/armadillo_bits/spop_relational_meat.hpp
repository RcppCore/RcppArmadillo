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


//! \addtogroup spop_relational
//! @{

// NOTE: in general, relational operations between sparse matrices and scalars
// NOTE: have the risk of producing sparse matrices with most elements as non-zeros.
// NOTE: these operations should only be used as an argument to the accu() function,
// NOTE: which aims to omit the generation of intermediate sparse matrices.


#undef operator_rel

#undef arma_applier_spmat_pre
#undef arma_applier_spmat_post


#define arma_applier_spmat_pre(operator_rel)\
  {\
  const uword zero_comp_val = (k operator_rel eT(0)) ? uword(1) : uword(0);\
  \
  const uword n_cols = A.n_cols;\
  const uword n_rows = A.n_rows;\
  const uword n_elem = A.n_elem;\
  \
  Mat<uword> tmp(n_rows, n_cols);\
  \
  typename SpMat<eT>::const_iterator it     = A.begin();\
  typename SpMat<eT>::const_iterator it_end = A.end();\
  \
  uword last_pos = 0;\
  uword offset   = 0;\
  \
  for(; it != it_end; ++it)\
    {\
    const uword cur_pos = it.row() + n_rows * it.col();\
    \
    if((cur_pos - last_pos) > offset)\
      {\
      arrayops::inplace_set( (tmp.memptr() + last_pos + offset), zero_comp_val, (cur_pos - last_pos - offset) );\
      }\
    \
    tmp.at(cur_pos) = (k operator_rel (*it)) ? uword(1) : uword(0);\
    \
    last_pos = cur_pos;\
    offset   = 1;\
    }\
  \
  if(last_pos < n_elem)\
    {\
    arrayops::inplace_set( (tmp.memptr() + last_pos + offset), zero_comp_val, (n_elem - last_pos - offset) );\
    }\
  \
  out = tmp;\
  }



#define arma_applier_spmat_post(operator_rel)\
  {\
  const uword zero_comp_val = (eT(0) operator_rel k) ? uword(1) : uword(0);\
  \
  const uword n_cols = A.n_cols;\
  const uword n_rows = A.n_rows;\
  const uword n_elem = A.n_elem;\
  \
  Mat<uword> tmp(n_rows, n_cols);\
  \
  typename SpMat<eT>::const_iterator it     = A.begin();\
  typename SpMat<eT>::const_iterator it_end = A.end();\
  \
  uword last_pos = 0;\
  uword offset   = 0;\
  \
  for(; it != it_end; ++it)\
    {\
    const uword cur_pos = it.row() + n_rows * it.col();\
    \
    if((cur_pos - last_pos) > offset)\
      {\
      arrayops::inplace_set( (tmp.memptr() + last_pos + offset), zero_comp_val, (cur_pos - last_pos - offset) );\
      }\
    \
    tmp.at(cur_pos) = ((*it) operator_rel k) ? uword(1) : uword(0);\
    \
    last_pos = cur_pos;\
    offset   = 1;\
    }\
  \
  if(last_pos < n_elem)\
    {\
    arrayops::inplace_set( (tmp.memptr() + last_pos + offset), zero_comp_val, (n_elem - last_pos - offset) );\
    }\
  \
  out = tmp;\
  }



template<typename T1>
inline
void
spop_rel_lt_pre::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_lt_pre>& X)
  {
  arma_debug_sigprint();
  
  // operation: scalar < spmat
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k > eT(0))
    {
    arma_debug_print("optimisation: k > 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( k < (*it) )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_pre( < );
    }
  }



template<typename T1>
inline
void
spop_rel_gt_pre::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_gt_pre>& X)
  {
  arma_debug_sigprint();
  
  // operation: scalar > spmat
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k < eT(0))
    {
    arma_debug_print("optimisation: k < 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( k > (*it) )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_pre( > );
    }
  }



template<typename T1>
inline
void
spop_rel_lteq_pre::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_lteq_pre>& X)
  {
  arma_debug_sigprint();
  
  // operation: scalar <= spmat
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k > eT(0))
    {
    arma_debug_print("optimisation: k > 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( k <= (*it) )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_pre( <= );
    }
  }



template<typename T1>
inline
void
spop_rel_gteq_pre::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_gteq_pre>& X)
  {
  arma_debug_sigprint();
  
  // operation: scalar >= spmat
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k < eT(0))
    {
    arma_debug_print("optimisation: k < 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( k >= (*it) )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_pre( >= );
    }
  }



template<typename T1>
inline
void
spop_rel_lt_post::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_lt_post>& X)
  {
  arma_debug_sigprint();
  
  // operation: spmat < scalar
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k < eT(0))
    {
    arma_debug_print("optimisation: k < 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( (*it) < k )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_post( < );
    }
  }



template<typename T1>
inline
void
spop_rel_gt_post::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_gt_post>& X)
  {
  arma_debug_sigprint();
  
  // operation: spmat > scalar
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k > eT(0))
    {
    arma_debug_print("optimisation: k > 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( (*it) > k )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_post( > );
    }
  }



template<typename T1>
inline
void
spop_rel_lteq_post::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_lteq_post>& X)
  {
  arma_debug_sigprint();
  
  // operation: spmat <= scalar
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k < eT(0))
    {
    arma_debug_print("optimisation: k < 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( (*it) <= k )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_post( <= );
    }
  }



template<typename T1>
inline
void
spop_rel_gteq_post::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_gteq_post>& X)
  {
  arma_debug_sigprint();
  
  // operation: spmat >= scalar
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k > eT(0))
    {
    arma_debug_print("optimisation: k > 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( (*it) >= k )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_post( >= );
    }
  }



template<typename T1>
inline
void
spop_rel_eq::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_eq>& X)
  {
  arma_debug_sigprint();
  
  // operation: spmat == scalar
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k != eT(0))
    {
    arma_debug_print("optimisation: k != 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      if( (*it) == k )  { tmp.at(it.row(), it.col()) = uword(1); }
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_post( == );
    }
  }



template<typename T1>
inline
void
spop_rel_noteq::apply(SpMat<uword>& out, const mtSpOp<uword, T1, spop_rel_noteq>& X)
  {
  arma_debug_sigprint();
  
  // operation: spmat != scalar
  
  typedef typename T1::elem_type eT;
  
  const eT k = X.aux;
  
  const unwrap_spmat<T1> U(X.m);
  const SpMat<eT>& A =   U.M;
  
  if(k == eT(0))
    {
    arma_debug_print("optimisation: k = 0");
    
    SpMat<uword> tmp(A.n_rows, A.n_cols);
    
    typename SpMat<eT>::const_iterator it     = A.begin();
    typename SpMat<eT>::const_iterator it_end = A.end();
    
    for(; it != it_end; ++it)
      {
      tmp.at(it.row(), it.col()) = uword(1);
      }
    
    out.steal_mem(tmp);
    }
  else
    {
    arma_applier_spmat_post( != );
    }
  }



#undef arma_applier_spmat_pre
#undef arma_applier_spmat_post



//! @}
