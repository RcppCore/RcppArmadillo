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


//! \addtogroup spop_vecnorm
//! @{



template<typename T1>
inline
void
spop_vecnorm::apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type,T1,spop_vecnorm>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const uword k   = expr.aux_uword_a;
  const uword dim = expr.aux_uword_b;
  
  arma_debug_check( (k   == 0), "vecnorm(): unsupported vector norm type"   );
  arma_debug_check( (dim >  1), "vecnorm(): parameter 'dim' must be 0 or 1" );
  
  const unwrap_spmat<T1> U(expr.m);
  const SpMat<eT>&   X = U.M;
  
  X.sync();
  
  if(dim == 0)
    {
    Mat<T> tmp;
    
    spop_vecnorm::apply_direct(tmp, X, k);
    
    out = tmp;
    }
  else
  if(dim == 1)
    {
      Mat< T> tmp;
    SpMat<eT> Xt;
    
    spop_strans::apply_noalias(Xt, X);
    
    spop_vecnorm::apply_direct(tmp, Xt, k);
    
    out = tmp.t();
    }
  }



template<typename eT>
inline
void
spop_vecnorm::apply_direct(Mat< typename get_pod_type<eT>::result >& out, const SpMat<eT>& X, const uword k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  out.zeros(1, X.n_cols);
  
  T* out_mem = out.memptr();
  
  for(uword col=0; col < X.n_cols; ++col)
    {
    const uword      col_offset = X.col_ptrs[col    ];
    const uword next_col_offset = X.col_ptrs[col + 1];
    
    const eT* start_ptr = &X.values[     col_offset];
    const eT*   end_ptr = &X.values[next_col_offset];
    
    const uword n_elem = end_ptr - start_ptr;
    
    T out_val = T(0);
    
    if(n_elem > 0)
      {
      const Col<eT> tmp(const_cast<eT*>(start_ptr), n_elem, false, false);
      
      const Proxy< Col<eT> > P(tmp);
      
      if(k == uword(1))  { out_val = op_norm::vec_norm_1(P); }
      if(k == uword(2))  { out_val = op_norm::vec_norm_2(P); }
      }
    
    out_mem[col] = out_val;
    }
  }



//



template<typename T1>
inline
void
spop_vecnorm_ext::apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type,T1,spop_vecnorm_ext>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const uword method_id = expr.aux_uword_a;
  const uword dim       = expr.aux_uword_b;
  
  arma_debug_check( (method_id == 0), "vecnorm(): unsupported vector norm type"   );
  arma_debug_check( (dim       >  1), "vecnorm(): parameter 'dim' must be 0 or 1" );
  
  const unwrap_spmat<T1> U(expr.m);
  const SpMat<eT>&   X = U.M;
  
  X.sync();
  
  if(dim == 0)
    {
    Mat<T> tmp;
    
    spop_vecnorm_ext::apply_direct(tmp, X, method_id);
    
    out = tmp;
    }
  else
  if(dim == 1)
    {
      Mat< T> tmp;
    SpMat<eT> Xt;
    
    spop_strans::apply_noalias(Xt, X);
    
    spop_vecnorm_ext::apply_direct(tmp, Xt, method_id);
    
    out = tmp.t();
    }
  }



template<typename eT>
inline
void
spop_vecnorm_ext::apply_direct(Mat< typename get_pod_type<eT>::result >& out, const SpMat<eT>& X, const uword method_id)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  out.zeros(1, X.n_cols);
  
  T* out_mem = out.memptr();
  
  for(uword col=0; col < X.n_cols; ++col)
    {
    const uword      col_offset = X.col_ptrs[col    ];
    const uword next_col_offset = X.col_ptrs[col + 1];
    
    const eT* start_ptr = &X.values[     col_offset];
    const eT*   end_ptr = &X.values[next_col_offset];
    
    const uword n_elem = end_ptr - start_ptr;
    
    T out_val = T(0);
    
    if(n_elem > 0)
      {
      const Col<eT> tmp(const_cast<eT*>(start_ptr), n_elem, false, false);
      
      const Proxy< Col<eT> > P(tmp);
      
      if(method_id == uword(1))
        {
        out_val = op_norm::vec_norm_max(P);
        }
      else
      if(method_id == uword(2))
        {
        const T tmp_val = op_norm::vec_norm_min(P);
        
        out_val = (n_elem < X.n_rows) ? T((std::min)(T(0), tmp_val)) : T(tmp_val);
        }
      }
    
    out_mem[col] = out_val;
    }
  }



//! @}
