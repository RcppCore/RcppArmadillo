// Copyright (C) 2012-2015 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup spop_misc
//! @{



namespace priv
  {
  template<typename eT>
  struct functor_scalar_times
    {
    const eT k;
    
    functor_scalar_times(const eT in_k) : k(in_k) {}
    
    arma_inline eT operator()(const eT val) const { return val * k; }
    };
  }



template<typename T1>
inline
void
spop_scalar_times::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(in.aux != eT(0))
    {
    out.init_xform(in.m, priv::functor_scalar_times<eT>(in.aux));
    }
  else
    {
    const SpProxy<T1> P(in.m);
    
    out.zeros( P.get_n_rows(), P.get_n_cols() );
    }
  }



namespace priv
  {
  struct functor_square
    {
    template<typename eT>
    arma_inline eT operator()(const eT val) const { return val*val; }
    };
  }



template<typename T1>
inline
void
spop_square::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_square>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform(in.m, priv::functor_square());
  }



namespace priv
  {
  struct functor_sqrt
    {
    template<typename eT>
    arma_inline eT operator()(const eT val) const { return eop_aux::sqrt(val); }
    };
  }



template<typename T1>
inline
void
spop_sqrt::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_sqrt>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform(in.m, priv::functor_sqrt());
  }



namespace priv
  {
  struct functor_abs
    {
    template<typename eT>
    arma_inline eT operator()(const eT val) const { return eop_aux::arma_abs(val); }
    };
  }



template<typename T1>
inline
void
spop_abs::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_abs>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform(in.m, priv::functor_abs());
  }



namespace priv
  {
  struct functor_cx_abs
    {
    template<typename T>
    arma_inline T operator()(const std::complex<T>& val) const { return std::abs(val); }
    };
  }



template<typename T1>
inline
void
spop_cx_abs::apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_cx_abs>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform_mt(in.m, priv::functor_cx_abs());
  }



namespace priv
  {
  struct functor_real
    {
    template<typename T>
    arma_inline T operator()(const std::complex<T>& val) const { return val.real(); }
    };
  }



template<typename T1>
inline
void
spop_real::apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_real>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform_mt(in.m, priv::functor_real());
  }



namespace priv
  {
  struct functor_imag
    {
    template<typename T>
    arma_inline T operator()(const std::complex<T>& val) const { return val.imag(); }
    };
  }



template<typename T1>
inline
void
spop_imag::apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_imag>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform_mt(in.m, priv::functor_imag());
  }



namespace priv
  {
  struct functor_conj
    {
    template<typename eT>
    arma_inline eT operator()(const eT val) const { return eop_aux::conj(val); }
    };
  }



template<typename T1>
inline
void
spop_conj::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_conj>& in)
  {
  arma_extra_debug_sigprint();
  
  out.init_xform(in.m, priv::functor_conj());
  }



template<typename T1>
inline
void
spop_repmat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_repmat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(in.m);
  const SpMat<eT>& X =   U.M;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword copies_per_row = in.aux_uword_a;
  const uword copies_per_col = in.aux_uword_b;
  
  // out.set_size(X_n_rows * copies_per_row, X_n_cols * copies_per_col);
  // 
  // const uword out_n_rows = out.n_rows;
  // const uword out_n_cols = out.n_cols;
  // 
  // if( (out_n_rows > 0) && (out_n_cols > 0) )
  //   {
  //   for(uword col = 0; col < out_n_cols; col += X_n_cols)
  //   for(uword row = 0; row < out_n_rows; row += X_n_rows)
  //     {
  //     out.submat(row, col, row+X_n_rows-1, col+X_n_cols-1) = X;
  //     }
  //   }
  
  SpMat<eT> tmp(X_n_rows * copies_per_row, X_n_cols);
  
  if(tmp.n_elem > 0)
    {
    for(uword row = 0; row < tmp.n_rows; row += X_n_rows)
      {
      tmp.submat(row, 0, row+X_n_rows-1, X_n_cols-1) = X;
      }
    }
  
  // tmp contains copies of the input matrix, so no need to check for aliasing
  
  out.set_size(X_n_rows * copies_per_row, X_n_cols * copies_per_col);
  
  const uword out_n_rows = out.n_rows;
  const uword out_n_cols = out.n_cols;
  
  if( (out_n_rows > 0) && (out_n_cols > 0) )
    {
    for(uword col = 0; col < out_n_cols; col += X_n_cols)
      {
      out.submat(0, col, out_n_rows-1, col+X_n_cols-1) = tmp;
      }
    }
  }



template<typename T1>
inline
void
spop_reshape::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_reshape>& in)
  {
  arma_extra_debug_sigprint();
  
  out = in.m;
  
  out.reshape(in.aux_uword_a, in.aux_uword_b);
  }



template<typename T1>
inline
void
spop_resize::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_resize>& in)
  {
  arma_extra_debug_sigprint();
  
  out = in.m;
  
  out.resize(in.aux_uword_a, in.aux_uword_b);
  }



template<typename T1>
inline
void
spop_diagmat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_diagmat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> p(in.m);
  
  if(p.is_alias(out) == false)
    {
    spop_diagmat::apply_noalias(out, p);
    }
  else
    {
    SpMat<eT> tmp;
    
    spop_diagmat::apply_noalias(tmp, p);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1>
inline
void
spop_diagmat::apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& p)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = p.get_n_rows();
  const uword n_cols = p.get_n_cols();
  
  const bool p_is_vec = (n_rows == 1) || (n_cols == 1);
  
  if(p_is_vec)    // generate a diagonal matrix out of a vector
    {
    const uword N = (n_rows == 1) ? n_cols : n_rows;
    
    out.zeros(N, N);
    
    if(p.get_n_nonzero() == 0)  { return; }
    
    typename SpProxy<T1>::const_iterator_type it     = p.begin();
    typename SpProxy<T1>::const_iterator_type it_end = p.end();
      
    if(n_cols == 1)
      {
      while(it != it_end)
        {
        const uword row = it.row();
        
        out.at(row,row) = (*it);
        
        ++it;
        }
      }
    else
    if(n_rows == 1)
      {
      while(it != it_end)
        {
        const uword col = it.col();
        
        out.at(col,col) = (*it);
        
        ++it;
        }
      }
    }
  else   // generate a diagonal matrix out of a matrix
    {
    arma_debug_check( (n_rows != n_cols), "diagmat(): given matrix is not square" );
    
    out.zeros(n_rows, n_rows);
    
    if(p.get_n_nonzero() == 0)  { return; }
    
    typename SpProxy<T1>::const_iterator_type it     = p.begin();
    typename SpProxy<T1>::const_iterator_type it_end = p.end();
      
    while(it != it_end)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      if(row == col)
        {
        out.at(row,row) = (*it);
        }
      
      ++it;
      }
    }
  }



//! @}
