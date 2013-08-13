// Copyright (C) 2008-2012 Conrad Sanderson
// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_diagvec
//! @{



template<typename T1>
inline
void
op_diagvec::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagvec>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword a = X.aux_uword_a;
  const uword b = X.aux_uword_b;
  
  const uword row_offset = (b >  0) ? a : 0;
  const uword col_offset = (b == 0) ? a : 0;
  
  const Proxy<T1> P(X.m);
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  arma_debug_check
    (
    ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
    "diagvec(): requested diagonal is out of bounds"
    );
  
  const uword len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  if(is_Mat<typename Proxy<T1>::stored_type>::value)
    {
    op_diagvec::apply_unwrap(out, P.Q, row_offset, col_offset, len);
    }
  else
    {
    if(P.is_alias(out) == false)
      {
      op_diagvec::apply_proxy(out, P, row_offset, col_offset, len);
      }
    else
      {
      Mat<eT> tmp;
      
      op_diagvec::apply_proxy(tmp, P, row_offset, col_offset, len);
      
      out.steal_mem(tmp);
      }
    }
  }



template<typename T1>
arma_hot
inline
void
op_diagvec::apply_unwrap(Mat<typename T1::elem_type>& out, const T1& X, const uword row_offset, const uword col_offset, const uword len)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp_A(X, out);
  const Mat<eT>& A =     tmp_A.M;
  
  out.set_size(len, 1);
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < len; i+=2, j+=2)
    {
    const eT tmp_i = A.at( i + row_offset, i + col_offset );
    const eT tmp_j = A.at( j + row_offset, j + col_offset );
    
    out_mem[i] = tmp_i;
    out_mem[j] = tmp_j;
    }
  
  if(i < len)
    {
    out_mem[i] = A.at( i + row_offset, i + col_offset );
    }
  }



template<typename T1>
arma_hot
inline
void
op_diagvec::apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword row_offset, const uword col_offset, const uword len)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out.set_size(len, 1);
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < len; i+=2, j+=2)
    {
    const eT tmp_i = P.at( i + row_offset, i + col_offset );
    const eT tmp_j = P.at( j + row_offset, j + col_offset );
    
    out_mem[i] = tmp_i;
    out_mem[j] = tmp_j;
    }
  
  if(i < len)
    {
    out_mem[i] = P.at( i + row_offset, i + col_offset );
    }
  }



//! @}
