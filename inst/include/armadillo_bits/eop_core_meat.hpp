// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
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



//
// matrices



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 n_rows = x.get_n_rows();
  const u32 n_cols = x.get_n_cols();
  const u32 n_elem = x.get_n_elem();
  
  out.set_size(n_rows, n_cols);
  
  if(is_generator<eop_type>::value == true)
    {
         if(is_same_type<eop_type, eop_ones_diag>::value == true) { out.eye();   }
    else if(is_same_type<eop_type, eop_ones_full>::value == true) { out.ones();  }
    else if(is_same_type<eop_type, eop_zeros    >::value == true) { out.zeros(); }
    else if(is_same_type<eop_type, eop_randu    >::value == true) { out.randu(); }
    else if(is_same_type<eop_type, eop_randn    >::value == true) { out.randn(); }
    }
  else
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] = tmp_i;
      out_mem[j] = tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] = eop_core<eop_type>::process(P[i], k);
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
  
  const u32 n_rows = x.get_n_rows();
  const u32 n_cols = x.get_n_cols();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "addition");
  
  if(is_generator<eop_type>::value == true)
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      const u32 N = (std::min)(n_rows, n_cols);
      
      for(u32 i=0; i<N; ++i)
        {
        out.at(i,i) += eT(1);
        }
      }
    else
      {
            eT* out_mem = out.memptr();
      const u32 n_elem  = out.n_elem;
      
      for(u32 i=0; i<n_elem; ++i)
        {
        out_mem[i] += eop_aux::generate<eT,eop_type>();
        }
      }
    }
  else
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    const u32     n_elem  = out.n_elem;
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] += tmp_i;
      out_mem[j] += tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] += eop_core<eop_type>::process(P[i], k);
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
  
  const u32 n_rows = x.get_n_rows();
  const u32 n_cols = x.get_n_cols();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "subtraction");
  
  if(is_generator<eop_type>::value == true)
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      const u32 N = (std::min)(n_rows, n_cols);
      
      for(u32 i=0; i<N; ++i)
        {
        out.at(i,i) -= eT(1);
        }
      }
    else
      {
            eT* out_mem = out.memptr();
      const u32 n_elem  = out.n_elem;
      
      for(u32 i=0; i<n_elem; ++i)
        {
        out_mem[i] -= eop_aux::generate<eT,eop_type>();
        }
      }
    }
  else
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    const u32     n_elem  = out.n_elem;
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] -= tmp_i;
      out_mem[j] -= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] -= eop_core<eop_type>::process(P[i], k);
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
  
  const u32 n_rows = x.get_n_rows();
  const u32 n_cols = x.get_n_cols();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "element-wise multiplication");
  
  if(is_generator<eop_type>::value == true)
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      const u32 N = (std::min)(n_rows, n_cols);
      
      for(u32 i=0; i<N; ++i)
        {
        for(u32 row=0;   row<i;      ++row) { out.at(row,i) = eT(0); }
        for(u32 row=i+1; row<n_rows; ++row) { out.at(row,i) = eT(0); }
        }
      }
    else
      {
            eT* out_mem = out.memptr();
      const u32 n_elem  = out.n_elem;
      
      for(u32 i=0; i<n_elem; ++i)
        {
        out_mem[i] *= eop_aux::generate<eT,eop_type>();
        }
      }
    }
  else
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    const u32     n_elem  = out.n_elem;
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] *= tmp_i;
      out_mem[j] *= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] *= eop_core<eop_type>::process(P[i], k);
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
  
  const u32 n_rows = x.get_n_rows();
  const u32 n_cols = x.get_n_cols();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "element-wise division");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = out.n_elem;
  
  if(is_generator<eop_type>::value == true)
    {
    if(is_same_type<eop_type, eop_ones_diag>::value == true)
      {
      const u32 N = (std::min)(n_rows, n_cols);
      
      for(u32 i=0; i<N; ++i)
        {
        const eT zero = eT(0);
        
        for(u32 row=0;   row<i;      ++row) { out.at(row,i) /= zero; }
        for(u32 row=i+1; row<n_rows; ++row) { out.at(row,i) /= zero; }
        }
      }
    else
      {
            eT* out_mem = out.memptr();
      const u32 n_elem  = out.n_elem;
      
      for(u32 i=0; i<n_elem; ++i)
        {
        out_mem[i] /= eop_aux::generate<eT,eop_type>();
        }
      }
    }
  else
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    const u32     n_elem  = out.n_elem;
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] /= tmp_i;
      out_mem[j] /= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] /= eop_core<eop_type>::process(P[i], k);
      }
    }
  }



//
// cubes



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 n_rows   = x.get_n_rows();
  const u32 n_cols   = x.get_n_cols();
  const u32 n_slices = x.get_n_slices();
  const u32 n_elem   = x.get_n_elem();
  
  out.set_size(n_rows, n_cols, n_slices);
  
  if(is_generator<eop_type>::value == true)
    {
         if(is_same_type<eop_type, eop_ones_full>::value == true) { out.ones();  }
    else if(is_same_type<eop_type, eop_zeros    >::value == true) { out.zeros(); }
    else if(is_same_type<eop_type, eop_randu    >::value == true) { out.randu(); }
    else if(is_same_type<eop_type, eop_randn    >::value == true) { out.randn(); }
    }
  else
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      out_mem[i] = eop_core<eop_type>::process(P[i], k);
      out_mem[j] = eop_core<eop_type>::process(P[j], k);
      }
    
    if(i < n_elem)
      {
      out_mem[i] = eop_core<eop_type>::process(P[i], k);
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_inplace_plus(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 n_rows   = x.get_n_rows();
  const u32 n_cols   = x.get_n_cols();
  const u32 n_slices = x.get_n_slices();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, n_rows, n_cols, n_slices, "addition");
  
  eT* out_mem = out.memptr();
  
  if(is_generator<eop_type>::value == true)
    {
          eT* out_mem = out.memptr();
    const u32 n_elem  = out.n_elem;
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] += eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    const u32     n_elem  = out.n_elem;
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] += tmp_i;
      out_mem[j] += tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] += eop_core<eop_type>::process(P[i], k);
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_inplace_minus(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 n_rows   = x.get_n_rows();
  const u32 n_cols   = x.get_n_cols();
  const u32 n_slices = x.get_n_slices();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, n_rows, n_cols, n_slices, "subtraction");
  
  eT* out_mem = out.memptr();
  
  if(is_generator<eop_type>::value == true)
    {
          eT* out_mem = out.memptr();
    const u32 n_elem  = out.n_elem;
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] -= eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    const u32     n_elem  = out.n_elem;
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] -= tmp_i;
      out_mem[j] -= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] -= eop_core<eop_type>::process(P[i], k);
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_inplace_schur(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 n_rows   = x.get_n_rows();
  const u32 n_cols   = x.get_n_cols();
  const u32 n_slices = x.get_n_slices();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, n_rows, n_cols, n_slices, "element-wise multiplication");
  
  eT* out_mem = out.memptr();
  
  if(is_generator<eop_type>::value == true)
    {
          eT* out_mem = out.memptr();
    const u32 n_elem  = out.n_elem;
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] *= eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    const u32     n_elem  = out.n_elem;
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] *= tmp_i;
      out_mem[j] *= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] *= eop_core<eop_type>::process(P[i], k);
      }
    }
  }



template<typename eop_type>
template<typename T1>
arma_hot
inline
void
eop_core<eop_type>::apply_inplace_div(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 n_rows   = x.get_n_rows();
  const u32 n_cols   = x.get_n_cols();
  const u32 n_slices = x.get_n_slices();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, n_rows, n_cols, n_slices, "element-wise division");
  
  eT* out_mem = out.memptr();
  
  if(is_generator<eop_type>::value == true)
    {
          eT* out_mem = out.memptr();
    const u32 n_elem  = out.n_elem;
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] /= eop_aux::generate<eT,eop_type>();
      }
    }
  else
    {
    typedef typename ProxyCube<T1>::ea_type ea_type;
    
    const eT      k       = x.aux;
          ea_type P       = x.P.get_ea();
          eT*     out_mem = out.memptr();
    const u32     n_elem  = out.n_elem;
    
    u32 i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = eop_core<eop_type>::process(P[i], k);
      const eT tmp_j = eop_core<eop_type>::process(P[j], k);
      
      out_mem[i] /= tmp_i;
      out_mem[j] /= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] /= eop_core<eop_type>::process(P[i], k);
      }
    }
  }



//
// common



template<typename eop_type>
template<typename eT>
arma_hot
arma_pure
arma_inline
eT
eop_core<eop_type>::process(const eT val, const eT k)
  {
  arma_ignore(val);
  arma_ignore(k);
  
  arma_stop("eop_core::process(): unhandled eop_type");
  return eT(0);
  }



template<> template<typename eT> arma_hot arma_const arma_inline eT
eop_core<eop_scalar_plus      >::process(const eT val, const eT k) { return val + k;                  }

template<> template<typename eT> arma_hot arma_const arma_inline eT
eop_core<eop_scalar_minus_pre >::process(const eT val, const eT k) { return k - val;                  }

template<> template<typename eT> arma_hot arma_const arma_inline eT
eop_core<eop_scalar_minus_post>::process(const eT val, const eT k) { return val - k;                  }

template<> template<typename eT> arma_hot arma_const arma_inline eT
eop_core<eop_scalar_times     >::process(const eT val, const eT k) { return val * k;                  }

template<> template<typename eT> arma_hot arma_const arma_inline eT
eop_core<eop_scalar_div_pre   >::process(const eT val, const eT k) { return k / val;                  }

template<> template<typename eT> arma_hot arma_const arma_inline eT
eop_core<eop_scalar_div_post  >::process(const eT val, const eT k) { return val / k;                  }

template<> template<typename eT> arma_hot arma_const arma_inline eT
eop_core<eop_square           >::process(const eT val, const eT  ) { return val*val;                  }

template<> template<typename eT> arma_hot arma_const arma_inline eT
eop_core<eop_neg              >::process(const eT val, const eT  ) { return eop_aux::neg(val);        }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_sqrt             >::process(const eT val, const eT  ) { return eop_aux::sqrt(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_log              >::process(const eT val, const eT  ) { return eop_aux::log(val);        }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_log2             >::process(const eT val, const eT  ) { return eop_aux::log2(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_log10            >::process(const eT val, const eT  ) { return eop_aux::log10(val);      }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_trunc_log        >::process(const eT val, const eT  ) { return    arma::trunc_log(val);  }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_exp              >::process(const eT val, const eT  ) { return eop_aux::exp(val);        }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_exp2             >::process(const eT val, const eT  ) { return eop_aux::exp2(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_exp10            >::process(const eT val, const eT  ) { return eop_aux::exp10(val);      }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_trunc_exp        >::process(const eT val, const eT  ) { return    arma::trunc_exp(val);  }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_cos              >::process(const eT val, const eT  ) { return eop_aux::cos(val);        }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_sin              >::process(const eT val, const eT  ) { return eop_aux::sin(val);        }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_tan              >::process(const eT val, const eT  ) { return eop_aux::tan(val);        }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_acos             >::process(const eT val, const eT  ) { return eop_aux::acos(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_asin             >::process(const eT val, const eT  ) { return eop_aux::asin(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_atan             >::process(const eT val, const eT  ) { return eop_aux::atan(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_cosh             >::process(const eT val, const eT  ) { return eop_aux::cosh(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_sinh             >::process(const eT val, const eT  ) { return eop_aux::sinh(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_tanh             >::process(const eT val, const eT  ) { return eop_aux::tanh(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_acosh            >::process(const eT val, const eT  ) { return eop_aux::acosh(val);      }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_asinh            >::process(const eT val, const eT  ) { return eop_aux::asinh(val);      }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_atanh            >::process(const eT val, const eT  ) { return eop_aux::atanh(val);      }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_eps              >::process(const eT val, const eT  ) { return eop_aux::direct_eps(val); }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_abs              >::process(const eT val, const eT  ) { return eop_aux::arma_abs(val);   }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_conj             >::process(const eT val, const eT  ) { return eop_aux::conj(val);       }

template<> template<typename eT> arma_hot arma_pure arma_inline eT
eop_core<eop_pow              >::process(const eT val, const eT k) { return eop_aux::pow(val, k);     }



//! @}
