// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup subview_elem1
//! @{


template<typename eT, typename T1>
inline
subview_elem1<eT,T1>::~subview_elem1()
  {
  arma_debug_sigprint();
  }


template<typename eT, typename T1>
arma_inline
subview_elem1<eT,T1>::subview_elem1(const Mat<eT>& in_m, const Base<uword,T1>& in_a)
  : m(in_m)
  , a(in_a)
  {
  arma_debug_sigprint();
  }



template<typename eT, typename T1>
arma_inline
subview_elem1<eT,T1>::subview_elem1(const Cube<eT>& in_q, const Base<uword,T1>& in_a)
  : fake_m( const_cast< eT* >(in_q.memptr()), in_q.n_elem, 1, false )
  ,      m( fake_m )
  ,      a( in_a   )
  {
  arma_debug_sigprint();
  }



template<typename eT, typename T1>
template<typename op_type>
inline
void
subview_elem1<eT,T1>::inplace_op(const eT val)
  {
  arma_debug_sigprint();
  
  Mat<eT>& m_local = const_cast< Mat<eT>& >(m);
  
        eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  if(strip_op_find_default<T1>::do_op_find_default)
    {
    const strip_op_find_default<T1> strip(a.get_ref());
    
    constexpr bool has_sv = strip_op_find_default<T1>::stored_type::has_subview;
    
    if( (has_sv == false) || ((has_sv == true) && (strip.M.is_alias(m_local) == false)) )
      {
      if(has_sv == false)  { arma_debug_print("op_find_default optimisation; has_sv = false"); }
      if(has_sv == true )  { arma_debug_print("op_find_default optimisation; has_sv = true" ); }
      
      bool mi_bad = false;
      
      if(is_same_type<op_type, op_internal_equ  >::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { m_mem[mi] =  val; } else { mi_bad = true; } }; op_find_aux::apply(modifier, strip.M); }
      if(is_same_type<op_type, op_internal_plus >::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { m_mem[mi] += val; } else { mi_bad = true; } }; op_find_aux::apply(modifier, strip.M); }
      if(is_same_type<op_type, op_internal_minus>::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { m_mem[mi] -= val; } else { mi_bad = true; } }; op_find_aux::apply(modifier, strip.M); }
      if(is_same_type<op_type, op_internal_schur>::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { m_mem[mi] *= val; } else { mi_bad = true; } }; op_find_aux::apply(modifier, strip.M); }
      if(is_same_type<op_type, op_internal_div  >::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { m_mem[mi] /= val; } else { mi_bad = true; } }; op_find_aux::apply(modifier, strip.M); }
      
      arma_conform_check_bounds( mi_bad, "Mat::elem(): index out of bounds" );
      
      return;
      }
    }
  
  const unwrap_check_mixed<T1> U(a.get_ref(), m_local);
  const umat& aa = U.M;
  
  if(resolves_to_vector<T1>::no)
    {
    arma_conform_check( ( (aa.is_vec() == false) && (aa.is_empty() == false) ), "Mat::elem(): given object must be a vector" );
    }
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  bool ii_jj_bad = false;
  
  uword iq,jq;
  for(iq=0, jq=1; jq < aa_n_elem; iq+=2, jq+=2)
    {
    const uword ii = aa_mem[iq];
    const uword jj = aa_mem[jq];
    
    if( (ii < m_n_elem) && (jj < m_n_elem) )
      {
      if(is_same_type<op_type, op_internal_equ  >::yes) { m_mem[ii] =  val; m_mem[jj] =  val; }
      if(is_same_type<op_type, op_internal_plus >::yes) { m_mem[ii] += val; m_mem[jj] += val; }
      if(is_same_type<op_type, op_internal_minus>::yes) { m_mem[ii] -= val; m_mem[jj] -= val; }
      if(is_same_type<op_type, op_internal_schur>::yes) { m_mem[ii] *= val; m_mem[jj] *= val; }
      if(is_same_type<op_type, op_internal_div  >::yes) { m_mem[ii] /= val; m_mem[jj] /= val; }
      }
    else
      {
      ii_jj_bad = true;
      }
    }
  
  if(iq < aa_n_elem)
    {
    const uword ii = aa_mem[iq];
    
    if(ii < m_n_elem)
      {
      if(is_same_type<op_type, op_internal_equ  >::yes) { m_mem[ii] =  val; }
      if(is_same_type<op_type, op_internal_plus >::yes) { m_mem[ii] += val; }
      if(is_same_type<op_type, op_internal_minus>::yes) { m_mem[ii] -= val; }
      if(is_same_type<op_type, op_internal_schur>::yes) { m_mem[ii] *= val; }
      if(is_same_type<op_type, op_internal_div  >::yes) { m_mem[ii] /= val; }
      }
    else
      {
      ii_jj_bad = true;
      }
    }
  
  arma_conform_check_bounds( ii_jj_bad, "Mat::elem(): index out of bounds" );
  }



template<typename eT, typename T1>
template<typename op_type, typename T2>
inline
void
subview_elem1<eT,T1>::inplace_op(const Base<eT,T2>& expr)
  {
  arma_debug_sigprint();
  
  Mat<eT>& m_local = const_cast< Mat<eT>& >(m);
  
        eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  typedef typename quasi_unwrap<T2>::stored_type U1_M_type;
  
  const quasi_unwrap<T2       > U1(expr.get_ref());
  const unwrap_check<U1_M_type> U2(U1.M, U1.is_alias(m_local));
  
  const eT*   X_mem    = U2.M.memptr();
  const uword X_n_elem = U2.M.n_elem;
  
  if(strip_op_find_default<T1>::do_op_find_default)
    {
    const strip_op_find_default<T1> strip(a.get_ref());
    
    constexpr bool has_sv = strip_op_find_default<T1>::stored_type::has_subview;
    
    if( (has_sv == false) || ((has_sv == true) && (strip.M.is_alias(m_local) == false)) )
      {
      if(has_sv == false)  { arma_debug_print("op_find_default optimisation; has_sv = false"); }
      if(has_sv == true )  { arma_debug_print("op_find_default optimisation; has_sv = true" ); }
      
      bool mi_bad = false;
      
      uword Xi = 0;
      
      if(is_same_type<op_type, op_internal_equ  >::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { if(Xi < X_n_elem) { m_mem[mi] =  X_mem[Xi]; } } else { mi_bad = true; } ++Xi; }; op_find_aux::apply(modifier, strip.M); }
      if(is_same_type<op_type, op_internal_plus >::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { if(Xi < X_n_elem) { m_mem[mi] += X_mem[Xi]; } } else { mi_bad = true; } ++Xi; }; op_find_aux::apply(modifier, strip.M); }
      if(is_same_type<op_type, op_internal_minus>::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { if(Xi < X_n_elem) { m_mem[mi] -= X_mem[Xi]; } } else { mi_bad = true; } ++Xi; }; op_find_aux::apply(modifier, strip.M); }
      if(is_same_type<op_type, op_internal_schur>::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { if(Xi < X_n_elem) { m_mem[mi] *= X_mem[Xi]; } } else { mi_bad = true; } ++Xi; }; op_find_aux::apply(modifier, strip.M); }
      if(is_same_type<op_type, op_internal_div  >::yes)  { auto modifier = [&](const uword mi) { if(mi < m_n_elem) { if(Xi < X_n_elem) { m_mem[mi] /= X_mem[Xi]; } } else { mi_bad = true; } ++Xi; }; op_find_aux::apply(modifier, strip.M); }
      
      arma_conform_check_bounds( mi_bad, "Mat::elem(): index out of bounds" );
      
      arma_conform_check( (Xi != X_n_elem), "Mat::elem(): size mismatch" );
      
      return;
      }
    }
  
  const unwrap_check_mixed<T1> U(a.get_ref(), m_local);
  const umat& aa = U.M;
  
  if(resolves_to_vector<T1>::no)
    {
    arma_conform_check( ( (aa.is_vec() == false) && (aa.is_empty() == false) ), "Mat::elem(): given object must be a vector" );
    }
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  arma_conform_check( (aa_n_elem != X_n_elem), "Mat::elem(): size mismatch" );
  
  bool ii_jj_bad = false;
  
  uword iq,jq;
  for(iq=0, jq=1; jq < aa_n_elem; iq+=2, jq+=2)
    {
    const uword ii = aa_mem[iq];
    const uword jj = aa_mem[jq];
    
    if( (ii < m_n_elem) && (jj < m_n_elem) )
      {
      if(is_same_type<op_type, op_internal_equ  >::yes) { m_mem[ii] =  X_mem[iq]; m_mem[jj]  = X_mem[jq]; }
      if(is_same_type<op_type, op_internal_plus >::yes) { m_mem[ii] += X_mem[iq]; m_mem[jj] += X_mem[jq]; }
      if(is_same_type<op_type, op_internal_minus>::yes) { m_mem[ii] -= X_mem[iq]; m_mem[jj] -= X_mem[jq]; }
      if(is_same_type<op_type, op_internal_schur>::yes) { m_mem[ii] *= X_mem[iq]; m_mem[jj] *= X_mem[jq]; }
      if(is_same_type<op_type, op_internal_div  >::yes) { m_mem[ii] /= X_mem[iq]; m_mem[jj] /= X_mem[jq]; }
      }
    else
      {
      ii_jj_bad = true;
      }
    }
  
  if(iq < aa_n_elem)
    {
    const uword ii = aa_mem[iq];
    
    if(ii < m_n_elem)
      {
      if(is_same_type<op_type, op_internal_equ  >::yes) { m_mem[ii] =  X_mem[iq]; }
      if(is_same_type<op_type, op_internal_plus >::yes) { m_mem[ii] += X_mem[iq]; }
      if(is_same_type<op_type, op_internal_minus>::yes) { m_mem[ii] -= X_mem[iq]; }
      if(is_same_type<op_type, op_internal_schur>::yes) { m_mem[ii] *= X_mem[iq]; }
      if(is_same_type<op_type, op_internal_div  >::yes) { m_mem[ii] /= X_mem[iq]; }
      }
    else
      {
      ii_jj_bad = true;
      }
    }
  
  arma_conform_check_bounds( ii_jj_bad, "Mat::elem(): index out of bounds" );
  }



//
//



template<typename eT, typename T1>
arma_inline
const Op<subview_elem1<eT,T1>,op_htrans>
subview_elem1<eT,T1>::t() const
  {
  return Op<subview_elem1<eT,T1>,op_htrans>(*this);
  }



template<typename eT, typename T1>
arma_inline
const Op<subview_elem1<eT,T1>,op_htrans>
subview_elem1<eT,T1>::ht() const
  {
  return Op<subview_elem1<eT,T1>,op_htrans>(*this);
  }



template<typename eT, typename T1>
arma_inline
const Op<subview_elem1<eT,T1>,op_strans>
subview_elem1<eT,T1>::st() const
  {
  return Op<subview_elem1<eT,T1>,op_strans>(*this);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::replace(const eT old_val, const eT new_val)
  {
  arma_debug_sigprint();
  
  Mat<eT>& m_local = const_cast< Mat<eT>& >(m);
  
        eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  const unwrap_check_mixed<T1> U(a.get_ref(), m_local);
  const umat& aa = U.M;
  
  if(resolves_to_vector<T1>::no)
    {
    arma_conform_check( ( (aa.is_vec() == false) && (aa.is_empty() == false) ), "Mat::elem(): given object must be a vector" );
    }
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  if(arma_isnan(old_val))
    {
    for(uword iq=0; iq < aa_n_elem; ++iq)
      {
      const uword ii = aa_mem[iq];
      
      arma_conform_check_bounds( (ii >= m_n_elem), "Mat::elem(): index out of bounds" );
      
      eT& val = m_mem[ii];
      
      val = (arma_isnan(val)) ? new_val : val;
      }
    }
  else
    {
    for(uword iq=0; iq < aa_n_elem; ++iq)
      {
      const uword ii = aa_mem[iq];
      
      arma_conform_check_bounds( (ii >= m_n_elem), "Mat::elem(): index out of bounds" );
      
      eT& val = m_mem[ii];
      
      val = (val == old_val) ? new_val : val;
      }
    }
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::clean(const pod_type threshold)
  {
  arma_debug_sigprint();
  
  Mat<eT> tmp(*this);
  
  tmp.clean(threshold);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::clamp(const eT min_val, const eT max_val)
  {
  arma_debug_sigprint();
  
  Mat<eT> tmp(*this);
  
  tmp.clamp(min_val, max_val);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::fill(const eT val)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_equ>(val);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::zeros()
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_equ>(eT(0));
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::ones()
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_equ>(eT(1));
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::randu()
  {
  arma_debug_sigprint();
  
  Mat<eT>& m_local = const_cast< Mat<eT>& >(m);
  
        eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  const unwrap_check_mixed<T1> U(a.get_ref(), m_local);
  const umat& aa = U.M;
  
  if(resolves_to_vector<T1>::no)
    {
    arma_conform_check( ( (aa.is_vec() == false) && (aa.is_empty() == false) ), "Mat::elem(): given object must be a vector" );
    }
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  podarray<eT> tmp(aa_n_elem);
  
  eT* tmp_mem = tmp.memptr();
  
  arma_rng::randu<eT>::fill(tmp_mem, aa_n_elem);
  
  for(uword iq=0; iq < aa_n_elem; ++iq)
    {
    const uword ii = aa_mem[iq];
    
    arma_conform_check_bounds( (ii >= m_n_elem), "Mat::elem(): index out of bounds" );
    
    m_mem[ii] = tmp_mem[iq];
    }
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::randn()
  {
  arma_debug_sigprint();
  
  Mat<eT>& m_local = const_cast< Mat<eT>& >(m);
  
        eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  const unwrap_check_mixed<T1> U(a.get_ref(), m_local);
  const umat& aa = U.M;
  
  if(resolves_to_vector<T1>::no)
    {
    arma_conform_check( ( (aa.is_vec() == false) && (aa.is_empty() == false) ), "Mat::elem(): given object must be a vector" );
    }
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  podarray<eT> tmp(aa_n_elem);
  
  eT* tmp_mem = tmp.memptr();
  
  arma_rng::randn<eT>::fill(tmp_mem, aa_n_elem);
  
  for(uword iq=0; iq < aa_n_elem; ++iq)
    {
    const uword ii = aa_mem[iq];
    
    arma_conform_check_bounds( (ii >= m_n_elem), "Mat::elem(): index out of bounds" );
    
    m_mem[ii] = tmp_mem[iq];
    }
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator+= (const eT val)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_plus>(val);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator-= (const eT val)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_minus>(val);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator*= (const eT val)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_schur>(val);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator/= (const eT val)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_div>(val);
  }



//
//



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator= (const subview_elem1<eT,T2>& x)
  {
  arma_debug_sigprint();
  
  const Mat<eT> tmp(x);
  
  inplace_op<op_internal_equ>(tmp);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator= (const subview_elem1<eT,T1>& x)
  {
  arma_debug_sigprint();
  
  const Mat<eT> tmp(x);
  
  inplace_op<op_internal_equ>(tmp);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator= (const Base<eT,T2>& x)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_equ>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator+= (const Base<eT,T2>& x)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_plus>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator-= (const Base<eT,T2>& x)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_minus>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator%= (const Base<eT,T2>& x)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_schur>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator/= (const Base<eT,T2>& x)
  {
  arma_debug_sigprint();
  
  inplace_op<op_internal_div>(x);
  }



//
//



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::extract_noalias(Mat<eT>& out, const subview_elem1<eT,T1>& in)
  {
  arma_debug_sigprint();
  
  const eT*   m_mem    = in.m.memptr();
  const uword m_n_elem = in.m.n_elem;
  
  if(strip_op_find_default<T1>::do_op_find_default)
    {
    arma_debug_print("op_find_default optimisation");
    
    const strip_op_find_default<T1> strip(in.a.get_ref());
    
    Mat<eT> tmp(m_n_elem, 1, arma_nozeros_indicator());  // worst-case scenario
    
    eT* tmp_mem = tmp.memptr();
    
    uword count = 0;
    
    bool mi_bad = false;
    
    auto element_extractor = [&](const uword mi) { if((mi < m_n_elem) && (count < m_n_elem)) { tmp_mem[count] = m_mem[mi]; ++count; } else { mi_bad = true; } };
    
    op_find_aux::apply(element_extractor, strip.M);
    
    arma_conform_check_bounds( mi_bad, "Mat::elem(): index out of bounds" );
    
    out.steal_mem_col(tmp, count);
    
    return;
    }
  
  const quasi_unwrap<T1> U(in.a.get_ref());
  const umat& aa = U.M;
  
  if(resolves_to_vector<T1>::no)
    {
    arma_conform_check( ( (aa.is_vec() == false) && (aa.is_empty() == false) ), "Mat::elem(): given object must be a vector" );
    }
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  out.set_size(aa_n_elem, 1);
  
  eT* out_mem = out.memptr();
  
  bool ii_jj_bad = false;
  
  uword i,j;
  for(i=0, j=1; j<aa_n_elem; i+=2, j+=2)
    {
    const uword ii = aa_mem[i];
    const uword jj = aa_mem[j];
    
    eT m_ii_val;
    eT m_jj_val;
    
    if( (ii < m_n_elem) && (jj < m_n_elem) )
      {
      m_ii_val = m_mem[ii];
      m_jj_val = m_mem[jj];
      }
    else
      {
      ii_jj_bad = true;
      
      m_ii_val = eT(0);
      m_jj_val = eT(0);
      }
    
    out_mem[i] = m_ii_val;
    out_mem[j] = m_jj_val;
    }
  
  if(i < aa_n_elem)
    {
    const uword ii = aa_mem[i];
    
    eT m_ii_val;
    
    if(ii < m_n_elem)
      {
      m_ii_val = m_mem[ii];
      }
    else
      {
      ii_jj_bad = true;
      
      m_ii_val = eT(0);
      }
    
    out_mem[i] = m_ii_val;
    }
  
  arma_conform_check_bounds( ii_jj_bad, "Mat::elem(): index out of bounds" );
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::extract(Mat<eT>& out, const subview_elem1<eT,T1>& in)
  {
  arma_debug_sigprint();
  
  if(in.is_alias(out))
    {
    Mat<eT> tmp;
    
    subview_elem1<eT,T1>::extract_noalias(tmp, in);
    
    out.steal_mem(tmp);
    }
  else
    {
    subview_elem1<eT,T1>::extract_noalias(out, in);
    }
  }



template<typename eT, typename T1>
template<typename eT2>
inline
bool
subview_elem1<eT,T1>::is_alias(const Mat<eT2>& X) const
  {
  arma_debug_sigprint();
  
  return (m.is_alias(X) || a.get_ref().is_alias(X));
  }



//! @}
