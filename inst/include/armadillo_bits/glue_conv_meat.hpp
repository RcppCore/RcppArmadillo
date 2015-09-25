// Copyright (C) 2010-2015 Conrad Sanderson
// Copyright (C) 2010-2015 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup glue_conv
//! @{


//! rudimentary implementation of the convolution operation

template<typename eT>
inline
void
glue_conv::apply_noalias(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const bool A_is_col)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& h = (A.n_elem <= B.n_elem) ? A : B;
  const Mat<eT>& x = (A.n_elem <= B.n_elem) ? B : A;
  
  const uword   h_n_elem = h.n_elem;
  const uword   x_n_elem = x.n_elem;
  const uword out_n_elem = h_n_elem + x_n_elem - 1;
  
  if( (h_n_elem == 0) || (x_n_elem == 0) )  { out.reset(); return; }
  
  (A_is_col) ? out.set_size(out_n_elem, 1) : out.set_size(1, out_n_elem);
  
  const eT*   h_mem = h.memptr();
  const eT*   x_mem = x.memptr();
        eT* out_mem = out.memptr();
  
  
  for(uword out_i = 0; out_i < (h_n_elem-1); ++out_i)
    {
    eT acc = eT(0);
    
    uword h_i = out_i;
    
    for(uword x_i = 0; x_i <= out_i; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
    
    out_mem[out_i] = acc;
    }
  
  
  for(uword out_i = h_n_elem-1; out_i < out_n_elem - (h_n_elem-1); ++out_i)
    {
    eT acc = eT(0);
   
    uword h_i = h_n_elem - 1;
    
    for(uword x_i = out_i - h_n_elem + 1; x_i <= out_i; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
      
    out_mem[out_i] = acc;
    }
  
  
  for(uword out_i = out_n_elem - (h_n_elem-1); out_i < out_n_elem; ++out_i)
    {
    eT acc = eT(0);
    
    uword h_i = h_n_elem - 1;
    
    for(uword x_i = out_i - h_n_elem + 1; x_i < x_n_elem; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
    
    out_mem[out_i] = acc;
    }
  }



template<typename T1, typename T2>
inline
void
glue_conv::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_conv>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UA(X.A);
  const quasi_unwrap<T2> UB(X.B);
  
  arma_debug_check
    (
    ( ((UA.M.is_vec() == false) && (UA.M.is_empty() == false)) || ((UB.M.is_vec() == false) && (UB.M.is_empty() == false)) ),
    "conv(): given object is not a vector"
    );
  
  const bool A_is_col = ((T1::is_col) || (UA.M.n_cols == 1));
  
  if(UA.is_alias(out) || UB.is_alias(out))
    {
    Mat<eT> tmp;
    
    glue_conv::apply_noalias(tmp, UA.M, UB.M, A_is_col);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_conv::apply_noalias(out, UA.M, UB.M, A_is_col);
    }
  }



//! @}
