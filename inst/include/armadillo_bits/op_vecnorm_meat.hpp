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


//! \addtogroup op_vecnorm
//! @{



template<typename T1>
inline
void
op_vecnorm::apply(Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_vecnorm>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type  in_eT;
  typedef typename T1::pod_type  out_eT;
  
  const quasi_unwrap<T1> U(in.m);
  const Mat<in_eT>&  X = U.M;
  
  const uword k    = in.aux_uword_a;
  const uword dim  = in.aux_uword_b;
  
  arma_debug_check( (k   == 0), "vecnorm(): unsupported vector norm type"   );
  arma_debug_check( (dim >  1), "vecnorm(): parameter 'dim' must be 0 or 1" );
  
  if(U.is_alias(out))
    {
    Mat<out_eT> tmp;
    
    op_vecnorm::apply_noalias(tmp, X, k, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_vecnorm::apply_noalias(out, X, k, dim);
    }
  }




template<typename in_eT>
inline
void
op_vecnorm::apply_noalias(Mat<typename get_pod_type<in_eT>::result>& out, const Mat<in_eT>& X, const uword k, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<in_eT>::result out_eT;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_vecnorm::apply(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols);
    
    if(X_n_rows > 0)
      {
      out_eT* out_mem = out.memptr();
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        op_vecnorm::apply_rawmem( out_mem[col], X.colptr(col), X_n_rows, k );
        }
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_vecnorm::apply(): dim = 1");
    
    out.set_size(X_n_rows, (X_n_cols > 0) ? 1 : 0);
    
    if(X_n_cols > 0)
      {
      podarray<in_eT> dat(X_n_cols);
      
      in_eT*  dat_mem = dat.memptr();
      out_eT* out_mem = out.memptr();
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        dat.copy_row(X, row);
        
        op_vecnorm::apply_rawmem( out_mem[row], dat_mem, X_n_cols, k );
        }
      }
    }
  }



template<typename in_eT>
inline
void
op_vecnorm::apply_rawmem(typename get_pod_type<in_eT>::result& out_val, const in_eT* mem, const uword N, const uword k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<in_eT>::result out_eT;
  
  const Col<in_eT> tmp(const_cast<in_eT*>(mem), N, false, false);
  
  const Proxy< Col<in_eT> > P(tmp);
  
  if(P.get_n_elem() == 0)  { out_val = out_eT(0); return; }
  
  if(k == uword(1))  { out_val = op_norm::vec_norm_1(P); return; }
  if(k == uword(2))  { out_val = op_norm::vec_norm_2(P); return; }
  
  out_val = op_norm::vec_norm_k(P, int(k));
  }



//



template<typename T1>
inline
void
op_vecnorm_ext::apply(Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_vecnorm_ext>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type  in_eT;
  typedef typename T1::pod_type  out_eT;
  
  const quasi_unwrap<T1> U(in.m);
  const Mat<in_eT>&  X = U.M;
  
  const uword method_id = in.aux_uword_a;
  const uword dim       = in.aux_uword_b;
  
  arma_debug_check( (method_id == 0), "vecnorm(): unsupported vector norm type"   );
  arma_debug_check( (dim       >  1), "vecnorm(): parameter 'dim' must be 0 or 1" );
  
  if(U.is_alias(out))
    {
    Mat<out_eT> tmp;
    
    op_vecnorm_ext::apply_noalias(tmp, X, method_id, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_vecnorm_ext::apply_noalias(out, X, method_id, dim);
    }
  }




template<typename in_eT>
inline
void
op_vecnorm_ext::apply_noalias(Mat<typename get_pod_type<in_eT>::result>& out, const Mat<in_eT>& X, const uword method_id, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<in_eT>::result out_eT;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_vecnorm_ext::apply(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols);
    
    if(X_n_rows > 0)
      {
      out_eT* out_mem = out.memptr();
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        op_vecnorm_ext::apply_rawmem( out_mem[col], X.colptr(col), X_n_rows, method_id );
        }
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_vecnorm_ext::apply(): dim = 1");
    
    out.set_size(X_n_rows, (X_n_cols > 0) ? 1 : 0);
    
    if(X_n_cols > 0)
      {
      podarray<in_eT> dat(X_n_cols);
      
      in_eT*  dat_mem = dat.memptr();
      out_eT* out_mem = out.memptr();
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        dat.copy_row(X, row);
        
        op_vecnorm_ext::apply_rawmem( out_mem[row], dat_mem, X_n_cols, method_id );
        }
      }
    }
  }



template<typename in_eT>
inline
void
op_vecnorm_ext::apply_rawmem(typename get_pod_type<in_eT>::result& out_val, const in_eT* mem, const uword N, const uword method_id)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<in_eT>::result out_eT;
  
  const Col<in_eT> tmp(const_cast<in_eT*>(mem), N, false, false);
  
  const Proxy< Col<in_eT> > P(tmp);
  
  if(P.get_n_elem() == 0)  { out_val = out_eT(0); return; }
  
  if(method_id == uword(1))  { out_val = op_norm::vec_norm_max(P); return; }
  if(method_id == uword(2))  { out_val = op_norm::vec_norm_min(P); return; }
  
  out_val = out_eT(0);
  }



//! @}
