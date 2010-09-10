// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup eop_core
//! @{



template<typename eop_type>
template<typename T1>
arma_hot
arma_inline
typename T1::elem_type
eop_core<eop_type>::get_elem(const eOp<T1, eop_type>& x, const u32 i)
  {
  typedef typename T1::elem_type eT;
  
       if(is_generator<eop_type>::value                == true) { return eop_aux::generate<eT,eop_type>();                       }
  else if(is_same_type<eop_type, eop_ones_diag>::value == true) { return ((i % x.P.n_rows) == (i / x.P.n_rows)) ? eT(1) : eT(0); }
  else                                                          { return eop_core<eop_type>::process(x, x.P[i]);                 }
  }



template<typename eop_type>
template<typename T1>
arma_hot
arma_inline
typename T1::elem_type
eop_core<eop_type>::get_elem(const eOp<T1, eop_type>& x, const u32 row, const u32 col)
  {
  typedef typename T1::elem_type eT;
  
       if(is_generator<eop_type>::value                == true) { return eop_aux::generate<eT,eop_type>();                 }
  else if(is_same_type<eop_type, eop_ones_diag>::value == true) { return (row == col) ? eT(1) : eT(0);                     }
  else                                                          { return eop_core<eop_type>::process(x, x.P.at(row, col)); }
  }



template<typename eop_type>
template<typename T1>
arma_hot
arma_inline
typename T1::elem_type
eop_core<eop_type>::process(const eOp<T1, eop_type>& x, const typename T1::elem_type val)
  {
  typedef typename T1::elem_type eT;
  
  // the optimiser will keep only one return statement
  
       if(is_same_type<eop_type, eop_scalar_plus      >::value == true) { return val + x.aux;              }
  else if(is_same_type<eop_type, eop_scalar_minus_pre >::value == true) { return x.aux - val;              }
  else if(is_same_type<eop_type, eop_scalar_minus_post>::value == true) { return val - x.aux;              }
  else if(is_same_type<eop_type, eop_scalar_times     >::value == true) { return val * x.aux;              }
  else if(is_same_type<eop_type, eop_scalar_div_pre   >::value == true) { return x.aux / val;              }
  else if(is_same_type<eop_type, eop_scalar_div_post  >::value == true) { return val / x.aux;              }
  else if(is_same_type<eop_type, eop_square           >::value == true) { return val*val;                  }
  else if(is_same_type<eop_type, eop_neg              >::value == true) { return eop_aux::neg(val);        }
  else if(is_same_type<eop_type, eop_sqrt             >::value == true) { return eop_aux::sqrt(val);       }
  else if(is_same_type<eop_type, eop_log10            >::value == true) { return eop_aux::log10(val);      }
  else if(is_same_type<eop_type, eop_log              >::value == true) { return eop_aux::log(val);        }
  else if(is_same_type<eop_type, eop_trunc_log        >::value == true) { return    arma::trunc_log(val);  }
  else if(is_same_type<eop_type, eop_exp              >::value == true) { return eop_aux::exp(val);        }
  else if(is_same_type<eop_type, eop_trunc_exp        >::value == true) { return    arma::trunc_exp(val);  }
  else if(is_same_type<eop_type, eop_cos              >::value == true) { return eop_aux::cos(val);        }
  else if(is_same_type<eop_type, eop_sin              >::value == true) { return eop_aux::sin(val);        }
  else if(is_same_type<eop_type, eop_tan              >::value == true) { return eop_aux::tan(val);        }
  else if(is_same_type<eop_type, eop_acos             >::value == true) { return eop_aux::acos(val);       }
  else if(is_same_type<eop_type, eop_asin             >::value == true) { return eop_aux::asin(val);       }
  else if(is_same_type<eop_type, eop_atan             >::value == true) { return eop_aux::atan(val);       }
  else if(is_same_type<eop_type, eop_cosh             >::value == true) { return eop_aux::cosh(val);       }
  else if(is_same_type<eop_type, eop_sinh             >::value == true) { return eop_aux::sinh(val);       }
  else if(is_same_type<eop_type, eop_tanh             >::value == true) { return eop_aux::tanh(val);       }
  else if(is_same_type<eop_type, eop_acosh            >::value == true) { return eop_aux::acosh(val);      }
  else if(is_same_type<eop_type, eop_asinh            >::value == true) { return eop_aux::asinh(val);      }
  else if(is_same_type<eop_type, eop_atanh            >::value == true) { return eop_aux::atanh(val);      }
  else if(is_same_type<eop_type, eop_eps              >::value == true) { return eop_aux::direct_eps(val); }
  else if(is_same_type<eop_type, eop_abs              >::value == true) { return eop_aux::arma_abs(val);   }
  else if(is_same_type<eop_type, eop_conj             >::value == true) { return eop_aux::conj(val);       }
  else if(is_same_type<eop_type, eop_pow              >::value == true) { return eop_aux::pow(val, x.aux); }
  else if(is_same_type<eop_type, eop_pow_int          >::value == true)
    {
    const int exponent = (x.aux_u32_b == 0) ? int(x.aux_u32_a) : -int(x.aux_u32_a);
    
    return eop_aux::pow_int(val, exponent);
    }
  else
    {
    arma_stop("eop_core::process(): unhandled eop_type");
    return eT(0);
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
arma_inline
void
eop_core<eop_type>::apply(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  if(is_Mat<T1>::value == true)
    {
    eop_core<eop_type>::apply_unwrap(out, x);
    }
  else
    {
    eop_core<eop_type>::apply_proxy(out, x);
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_proxy(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  // eop_type::apply_proxy() function is not allowed to unwrap things
  // (in order to get the input into a common format).
  // the proxy class is already providing objects with element access
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P = x.P;
  
  out.set_size(P.n_rows, P.n_cols);
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;
    
  if(is_generator<eop_type>::value == true)
    {
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      out_mem[i] = eop_aux::generate<eT,eop_type>();
      out_mem[j] = eop_aux::generate<eT,eop_type>();
      }
    
    if(i < n_elem)
      {
      out_mem[i] = eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      for(u32 col=0; col<P.n_rows; ++col)
        {
        for(u32 row=0; row<col; ++row)          { out.at(row,col) = eT(0); }
          
        out.at(col,col) = eT(1);
        
        for(u32 row=col+1; row<P.n_rows; ++row) { out.at(row,col) = eT(0); }
        }
      }
    else
      {
      u32 i,j;
      
      for(i=0, j=1; j<n_elem; i+=2, j+=2)
        {
        out_mem[i] = eop_core<eop_type>::process(x, P[i]);
        out_mem[j] = eop_core<eop_type>::process(x, P[j]);
        }
      
      if(i < n_elem)
        {
        out_mem[i] = eop_core<eop_type>::process(x, P[i]);
        }
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_unwrap(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P = x.P;
  
//   cout << "*** P.n_rows = " << P.n_rows << endl;
//   cout << "*** P.n_cols = " << P.n_cols << endl;
  
  out.set_size(P.n_rows, P.n_cols);
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;
    
  const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
  const Mat<eT>& A = tmp.M;
  
  if(is_generator<eop_type>::value == true)
    {
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      out_mem[i] = eop_aux::generate<eT,eop_type>();
      out_mem[j] = eop_aux::generate<eT,eop_type>();
      }
    
    if(i < n_elem)
      {
      out_mem[i] = eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      for(u32 col=0; col<P.n_rows; ++col)
        {
        for(u32 row=0; row<col; ++row)          { out.at(row,col) = eT(0); }
          
        out.at(col,col) = eT(1);
        
        for(u32 row=col+1; row<P.n_rows; ++row) { out.at(row,col) = eT(0); }
        }
      }
    else
      {
      const eT* A_mem = A.memptr();
      
      u32 i,j;
      
      for(i=0, j=1; j<n_elem; i+=2, j+=2)
        {
        out_mem[i] = eop_core<eop_type>::process(x, A_mem[i]);
        out_mem[j] = eop_core<eop_type>::process(x, A_mem[j]);
        }
      
      if(i < n_elem)
        {
        out_mem[i] = eop_core<eop_type>::process(x, A_mem[i]);
        }
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_inplace_plus(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P = x.P;
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, P.n_rows, P.n_cols, "matrix addition");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;
    
  if(is_generator<eop_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] += eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      for(u32 row=0; row<P.n_rows; ++row)
        {
        out.at(row,row) += eT(1);
        }
      }
    else
      {
      for(u32 i=0; i<n_elem; ++i)
        {
        out_mem[i] += eop_core<eop_type>::process(x, P[i]);
        }
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_inplace_minus(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P = x.P;
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, P.n_rows, P.n_cols, "matrix subtraction");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;
    
  if(is_generator<eop_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] -= eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      for(u32 row=0; row<P.n_rows; ++row)
        {
        out.at(row,row) -= eT(1);
        }
      }
    else
      {
      for(u32 i=0; i<n_elem; ++i)
        {
        out_mem[i] -= eop_core<eop_type>::process(x, P[i]);
        }
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_inplace_schur(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P = x.P;
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, P.n_rows, P.n_cols, "element-wise matrix multiplication");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;
  
  if(is_generator<eop_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] *= eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      for(u32 col=0; col<P.n_rows; ++col)
        {
        for(u32 row=0;     row<col;      ++row) { out.at(row,col) = eT(0); }
        for(u32 row=col+1; row<P.n_rows; ++row) { out.at(row,col) = eT(0); }
        }
      }
    else
      {
      for(u32 i=0; i<n_elem; ++i)
        {
        out_mem[i] *= eop_core<eop_type>::process(x, P[i]);
        }
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_inplace_div(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P = x.P;
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, P.n_rows, P.n_cols, "element-wise matrix division");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;
  
  if(is_generator<eop_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] /= eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      for(u32 col=0; col<P.n_rows; ++col)
        {
        const eT zero = eT(0);
        
        for(u32 row=0;     row<col;      ++row) { out.at(row,col) /= zero; }
        for(u32 row=col+1; row<P.n_rows; ++row) { out.at(row,col) /= zero; }
        }
      }
    else
      {
      for(u32 i=0; i<n_elem; ++i)
        {
        out_mem[i] /= eop_core<eop_type>::process(x, P[i]);
        }
      }
    }
  }



//! @}
