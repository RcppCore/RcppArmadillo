// Copyright (C) 2008-2013 Conrad Sanderson
// Copyright (C) 2008-2013 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Ryan Curtin
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_htrans
//! @{



template<typename eT>
arma_hot
arma_inline
void
op_htrans::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_mat_noalias(out, A);
  }



template<typename eT>
arma_hot
inline
void
op_htrans::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (A_n_cols == 1) || (A_n_rows == 1) )
    {
    const uword n_elem = A.n_elem;
    
    const eT* A_mem   = A.memptr();
          eT* out_mem = out.memptr();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = std::conj(A_mem[i]);
      }
    }
  else
    {
    for(uword in_row = 0; in_row < A_n_rows; ++in_row)
      {
      const uword out_col = in_row;
      
      for(uword in_col = 0; in_col < A_n_cols; ++in_col)
        {
        const uword out_row = in_col;
        out.at(out_row, out_col) = std::conj( A.at(in_row, in_col) );
        }
      }
    }
  }



template<typename eT>
arma_hot
arma_inline
void
op_htrans::apply_mat_inplace(Mat<eT>& out, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_mat_inplace(out);
  }



template<typename eT>
arma_hot
inline
void
op_htrans::apply_mat_inplace(Mat<eT>& out, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword n_rows = out.n_rows;
  const uword n_cols = out.n_cols;
    
  if(n_rows == n_cols)
    {
    arma_extra_debug_print("doing in-place hermitian transpose of a square matrix");
    
    for(uword col=0; col < n_cols; ++col)
      {
      eT* coldata = out.colptr(col);
      
      out.at(col,col) = std::conj( out.at(col,col) );
      
      for(uword row=(col+1); row < n_rows; ++row)
        {
        const eT val1 = std::conj(coldata[row]);
        const eT val2 = std::conj(out.at(col,row));
        
        out.at(col,row) = val1;
        coldata[row]    = val2;
        }
      }
    }
  else
    {
    Mat<eT> tmp;
    
    op_htrans::apply_mat_noalias(tmp, out);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT>
arma_hot
arma_inline
void
op_htrans::apply_mat(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_mat(out, A);
  }



template<typename eT>
arma_hot
inline
void
op_htrans::apply_mat(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(&out != &A)
    {
    op_htrans::apply_mat_noalias(out, A);
    }
  else
    {
    op_htrans::apply_mat_inplace(out);
    }
  }



template<typename T1>
arma_hot
inline
void
op_htrans::apply_proxy(Mat<typename T1::elem_type>& out, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X);
  
  // allow detection of in-place transpose
  if( (is_Mat<typename Proxy<T1>::stored_type>::value == true) && (Proxy<T1>::fake_mat == false) )
    {
    const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    
    op_htrans::apply_mat(out, tmp.M);
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    const bool is_alias = P.is_alias(out);
    
    if( (resolves_to_vector<T1>::value == true) && (Proxy<T1>::prefer_at_accessor == false) )
      {
      if(is_alias == false)
        {
        out.set_size(n_cols, n_rows);
        
        eT* out_mem = out.memptr();
        
        const uword n_elem = P.get_n_elem();
        
        typename Proxy<T1>::ea_type Pea = P.get_ea();
        
        for(uword i=0; i < n_elem; ++i)
          {
          out_mem[i] = std::conj(Pea[i]);
          }
        }
      else  // aliasing
        {
        Mat<eT> out2(n_cols, n_rows);
        
        eT* out_mem = out2.memptr();
        
        const uword n_elem = P.get_n_elem();
        
        typename Proxy<T1>::ea_type Pea = P.get_ea();
        
        for(uword i=0; i < n_elem; ++i)
          {
          out_mem[i] = std::conj(Pea[i]);
          }
        
        out.steal_mem(out2);
        }
      }
    else
      {
      if(is_alias == false)
        {
        out.set_size(n_cols, n_rows);
        
        for(uword k=0; k < n_cols; ++k)
        for(uword i=0; i < n_rows; ++i)
          {
          out.at(k,i) = std::conj(P.at(i,k));
          }
        }
      else // aliasing
        {
        Mat<eT> out2(n_cols, n_rows);
        
        for(uword k=0; k < n_cols; ++k)
        for(uword i=0; i < n_rows; ++i)
          {
          out2.at(k,i) = std::conj(P.at(i,k));
          }
        
        out.steal_mem(out2);
        }
      }
    }
  }



template<typename T1>
arma_hot
inline
void
op_htrans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_proxy(out, in.m);
  }



template<typename T1>
arma_hot
inline
void
op_htrans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_htrans::apply_proxy(out, in.m);
  }



template<typename T1>
arma_hot
inline
void
op_htrans::apply(Mat<typename T1::elem_type>& out, const Op< Op<T1, op_trimat>, op_htrans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m.m);
  const Mat<eT>& A = tmp.M;
  
  const bool upper = in.m.aux_uword_a;
  
  op_trimat::apply_htrans(out, A, upper);
  }



//
// op_htrans2



template<typename eT>
arma_hot
arma_inline
void
op_htrans2::apply_noalias(Mat<eT>& out, const Mat<eT>& A, const eT val)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (A_n_cols == 1) || (A_n_rows == 1) )
    {
    const uword n_elem = A.n_elem;
    
    const eT* A_mem   = A.memptr();
          eT* out_mem = out.memptr();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = val * std::conj(A_mem[i]);
      }
    }
  else
    {
    for(uword in_row = 0; in_row < A_n_rows; ++in_row)
      {
      const uword out_col = in_row;
      
      for(uword in_col = 0; in_col < A_n_cols; ++in_col)
        {
        const uword out_row = in_col;
        out.at(out_row, out_col) = val * std::conj( A.at(in_row, in_col) );
        }
      }
    }
  }



template<typename eT>
arma_hot
inline
void
op_htrans2::apply(Mat<eT>& out, const Mat<eT>& A, const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(&out != &A)
    {
    op_htrans2::apply_noalias(out, A, val);
    }
  else
    {
    const uword n_rows = out.n_rows;
    const uword n_cols = out.n_cols;
      
    if(n_rows == n_cols)
      {
      arma_extra_debug_print("doing in-place hermitian transpose of a square matrix");
      
      // TODO: do multiplication while swapping
      
      for(uword col=0; col < n_cols; ++col)
        {
        eT* coldata = out.colptr(col);
        
        out.at(col,col) = std::conj( out.at(col,col) );
        
        for(uword row=(col+1); row < n_rows; ++row)
          {
          const eT val1 = std::conj(coldata[row]);
          const eT val2 = std::conj(out.at(col,row));
          
          out.at(col,row) = val1;
          coldata[row]    = val2;
          }
        }
      
      arrayops::inplace_mul( out.memptr(), val, out.n_elem );
      }
    else
      {
      Mat<eT> tmp;
      op_htrans2::apply_noalias(tmp, A, val);
      
      out.steal_mem(tmp);
      }
    }
  }



template<typename T1>
arma_hot
inline
void
op_htrans2::apply_proxy(Mat<typename T1::elem_type>& out, const T1& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X);
  
  // allow detection of in-place transpose
  if( (is_Mat<typename Proxy<T1>::stored_type>::value == true) && (Proxy<T1>::fake_mat == false) )
    {
    const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    
    op_htrans2::apply(out, tmp.M, val);
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    const bool is_alias = P.is_alias(out);
    
    if( (resolves_to_vector<T1>::value == true) && (Proxy<T1>::prefer_at_accessor == false) )
      {
      if(is_alias == false)
        {
        out.set_size(n_cols, n_rows);
        
        eT* out_mem = out.memptr();
        
        const uword n_elem = P.get_n_elem();
        
        typename Proxy<T1>::ea_type Pea = P.get_ea();
        
        for(uword i=0; i < n_elem; ++i)
          {
          out_mem[i] = val * std::conj(Pea[i]);
          }
        }
      else  // aliasing
        {
        Mat<eT> out2(n_cols, n_rows);
        
        eT* out_mem = out2.memptr();
        
        const uword n_elem = P.get_n_elem();
        
        typename Proxy<T1>::ea_type Pea = P.get_ea();
        
        for(uword i=0; i < n_elem; ++i)
          {
          out_mem[i] = val * std::conj(Pea[i]);
          }
        
        out.steal_mem(out2);
        }
      }
    else
      {
      if(is_alias == false)
        {
        out.set_size(n_cols, n_rows);
        
        for(uword k=0; k < n_cols; ++k)
        for(uword i=0; i < n_rows; ++i)
          {
          out.at(k,i) = val * std::conj(P.at(i,k));
          }
        }
      else // aliasing
        {
        Mat<eT> out2(n_cols, n_rows);
        
        for(uword k=0; k < n_cols; ++k)
        for(uword i=0; i < n_rows; ++i)
          {
          out2.at(k,i) = val * std::conj(P.at(i,k));
          }
        
        out.steal_mem(out2);
        }
      }
    }
  }



template<typename T1>
arma_hot
inline
void
op_htrans2::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans2::apply_proxy(out, in.m, in.aux);
  }



template<typename T1>
arma_hot
inline
void
op_htrans2::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_htrans2::apply_proxy(out, in.m, in.aux);
  }



//! @}
