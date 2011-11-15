// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup subview_elem1
//! @{


template<typename eT, typename T1>
inline
subview_elem1<eT,T1>::~subview_elem1()
  {
  arma_extra_debug_sigprint();
  }


template<typename eT, typename T1>
arma_inline
subview_elem1<eT,T1>::subview_elem1(const Mat<eT>& in_m, const Base<uword,T1>& in_a)
  : m(in_m)
  , m_ptr(0)
  , a(in_a)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT, typename T1>
arma_inline
subview_elem1<eT,T1>::subview_elem1(Mat<eT>& in_m, const Base<uword,T1>& in_a)
  : m(in_m)
  , m_ptr(&in_m)
  , a(in_a)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT, typename T1>
template<typename op_type>
inline
void
subview_elem1<eT,T1>::inplace_op(const eT val)
  {
  Mat<eT>& m_local = *m_ptr;
  
        eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  const unwrap_check_mixed<T1> tmp(a.get_ref(), m_local);
  const umat& aa = tmp.M;
  
  arma_debug_check
    (
    ( aa.is_vec() == false ),
    "Mat::elem(): given object is not a vector"
    );
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  uword i,j;
  for(i=0, j=1; j<aa_n_elem; i+=2, j+=2)
    {
    const uword ii = aa_mem[i];
    const uword jj = aa_mem[j];
    
    arma_debug_check( ( (ii >= m_n_elem) || (jj >= m_n_elem) ), "Mat::elem(): index out of bounds" );
    
         if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { m_mem[ii] =  val; m_mem[jj] =  val; }
    else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { m_mem[ii] += val; m_mem[jj] += val; }
    else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { m_mem[ii] -= val; m_mem[jj] -= val; }
    else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { m_mem[ii] *= val; m_mem[jj] *= val; }
    else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { m_mem[ii] /= val; m_mem[jj] /= val; }
    }
  
  if(i < aa_n_elem)
    {
    const uword ii = aa_mem[i];
    
    arma_debug_check( (ii >= m_n_elem) , "Mat::elem(): index out of bounds" ); 
    
         if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { m_mem[ii] =  val; }
    else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { m_mem[ii] += val; }
    else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { m_mem[ii] -= val; }
    else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { m_mem[ii] *= val; }
    else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { m_mem[ii] /= val; }
    }
  }



template<typename eT, typename T1>
template<typename op_type, typename T2>
inline
void
subview_elem1<eT,T1>::inplace_op(const subview_elem1<eT,T2>& x)
  {
  subview_elem1<eT,T1>& t = *this;
  
  if(&(t.m) == &(x.m))
    {
    arma_extra_debug_print("subview_elem1::inplace_op(): aliasing detected");
    
    const Mat<eT> tmp(x);
    
         if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { t.operator= (tmp); }
    else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { t.operator+=(tmp); }
    else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { t.operator-=(tmp); }
    else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { t.operator%=(tmp); }
    else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { t.operator/=(tmp); }
    }
  else
    {
          Mat<eT>& t_m_local = *(t.m_ptr);
    const Mat<eT>& x_m_local = x.m;
    
    const unwrap_check_mixed<T1> t_tmp(t.a.get_ref(), t_m_local);
    const unwrap_check_mixed<T2> x_tmp(x.a.get_ref(), t_m_local);
    
    const umat& t_aa = t_tmp.M;
    const umat& x_aa = x_tmp.M;
    
    arma_debug_check
      (
      ( (t_aa.is_vec() == false) || (x_aa.is_vec() == false) ),
      "Mat::elem(): given object is not a vector"
      );
    
    const uword* t_aa_mem = t_aa.memptr();
    const uword* x_aa_mem = x_aa.memptr();
    
    const uword t_aa_n_elem = t_aa.n_elem;
    
    arma_debug_check( (t_aa_n_elem != x_aa.n_elem), "Mat::elem(): size mismatch" );
    
    
          eT*   t_m_mem    = t_m_local.memptr();
    const uword t_m_n_elem = t_m_local.n_elem;
    
    const eT*   x_m_mem    = x_m_local.memptr();
    const uword x_m_n_elem = x_m_local.n_elem;
    
    uword i,j;
    for(i=0, j=1; j<t_aa_n_elem; i+=2, j+=2)
      {
      const uword t_ii = t_aa_mem[i];
      const uword t_jj = t_aa_mem[j];
      
      const uword x_ii = x_aa_mem[i];
      const uword x_jj = x_aa_mem[j];
      
      arma_debug_check
        (
        (t_ii >= t_m_n_elem) || (t_jj >= t_m_n_elem) || (x_ii >= x_m_n_elem) || (x_jj >= x_m_n_elem),
        "Mat::elem(): index out of bounds"
        );
      
           if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { t_m_mem[t_ii]  = x_m_mem[x_ii]; t_m_mem[t_jj]  = x_m_mem[x_jj]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { t_m_mem[t_ii] += x_m_mem[x_ii]; t_m_mem[t_jj] += x_m_mem[x_jj]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { t_m_mem[t_ii] -= x_m_mem[x_ii]; t_m_mem[t_jj] -= x_m_mem[x_jj]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { t_m_mem[t_ii] *= x_m_mem[x_ii]; t_m_mem[t_jj] *= x_m_mem[x_jj]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { t_m_mem[t_ii] /= x_m_mem[x_ii]; t_m_mem[t_jj] /= x_m_mem[x_jj]; }
      }
    
    if(i < t_aa_n_elem)
      {
      const uword t_ii = t_aa_mem[i];
      const uword x_ii = x_aa_mem[i];
      
      arma_debug_check
        (
        ( (t_ii >= t_m_n_elem) || (x_ii >= x_m_n_elem) ),
        "Mat::elem(): index out of bounds"
        );
      
           if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { t_m_mem[t_ii]  = x_m_mem[x_ii]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { t_m_mem[t_ii] += x_m_mem[x_ii]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { t_m_mem[t_ii] -= x_m_mem[x_ii]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { t_m_mem[t_ii] *= x_m_mem[x_ii]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { t_m_mem[t_ii] /= x_m_mem[x_ii]; }
      }
    }
  }



template<typename eT, typename T1>
template<typename op_type, typename T2>
inline
void
subview_elem1<eT,T1>::inplace_op(const Base<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& m_local = *m_ptr;
  
        eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  const unwrap_check_mixed<T1> tmp(a.get_ref(), m_local);
  const umat& aa = tmp.M;
  
  arma_debug_check
    (
    ( aa.is_vec() == false ),
    "Mat::elem(): given object is not a vector"
    );
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  const Proxy<T2> P(x.get_ref());
  
  arma_debug_check( (aa_n_elem != P.get_n_elem()), "Mat::elem(): size mismatch" );
  
  if(P.is_alias(m) == false)
    {
    typename Proxy<T2>::ea_type X = P.get_ea();
    
    uword i,j;
    for(i=0, j=1; j<aa_n_elem; i+=2, j+=2)
      {
      const uword ii = aa_mem[i];
      const uword jj = aa_mem[j];
      
      arma_debug_check( ( (ii >= m_n_elem) || (jj >= m_n_elem) ), "Mat::elem(): index out of bounds" );
      
           if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { m_mem[ii] =  X[i]; m_mem[jj]  = X[j]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { m_mem[ii] += X[i]; m_mem[jj] += X[j]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { m_mem[ii] -= X[i]; m_mem[jj] -= X[j]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { m_mem[ii] *= X[i]; m_mem[jj] *= X[j]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { m_mem[ii] /= X[i]; m_mem[jj] /= X[j]; }
      }
    
    if(i < aa_n_elem)
      {
      const uword ii = aa_mem[i];
      
      arma_debug_check( (ii >= m_n_elem) , "Mat::elem(): index out of bounds" );
      
           if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { m_mem[ii] =  X[i]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { m_mem[ii] += X[i]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { m_mem[ii] -= X[i]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { m_mem[ii] *= X[i]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { m_mem[ii] /= X[i]; }
      }
    }
  else
    {
    arma_extra_debug_print("subview_elem1::inplace_op(): aliasing detected");
    
    const unwrap_check<typename Proxy<T2>::stored_type> tmp(P.Q, m_local);
    const Mat<eT>& M = tmp.M;
    
    const eT* X = M.memptr();
    
    uword i,j;
    for(i=0, j=1; j<aa_n_elem; i+=2, j+=2)
      {
      const uword ii = aa_mem[i];
      const uword jj = aa_mem[j];
      
      arma_debug_check( ( (ii >= m_n_elem) || (jj >= m_n_elem) ), "Mat::elem(): index out of bounds" );
      
           if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { m_mem[ii] =  X[i]; m_mem[jj]  = X[j]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { m_mem[ii] += X[i]; m_mem[jj] += X[j]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { m_mem[ii] -= X[i]; m_mem[jj] -= X[j]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { m_mem[ii] *= X[i]; m_mem[jj] *= X[j]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { m_mem[ii] /= X[i]; m_mem[jj] /= X[j]; }
      }
    
    if(i < aa_n_elem)
      {
      const uword ii = aa_mem[i];
      
      arma_debug_check( (ii >= m_n_elem) , "Mat::elem(): index out of bounds" );
      
           if(is_same_type<op_type, op_subview_elem_equ          >::value == true) { m_mem[ii] =  X[i]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { m_mem[ii] += X[i]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { m_mem[ii] -= X[i]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { m_mem[ii] *= X[i]; }
      else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { m_mem[ii] /= X[i]; }
      }
    }
  }



//
//



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_equ>(val);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::zeros()
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_equ>(eT(0));
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::ones()
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_equ>(eT(1));
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator+= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_plus>(val);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator-= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_minus>(val);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator*= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_schur>(val);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator/= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_div>(val);
  }



//
//



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator_equ(const subview_elem1<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_equ>(x);
  }




template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator= (const subview_elem1<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  (*this).operator_equ(x);
  }



//! work around compiler bugs
template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::operator= (const subview_elem1<eT,T1>& x)
  {
  arma_extra_debug_sigprint();
  
  (*this).operator_equ(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator+= (const subview_elem1<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_plus>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator-= (const subview_elem1<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_minus>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator%= (const subview_elem1<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_schur>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator/= (const subview_elem1<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_div>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator= (const Base<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_equ>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator+= (const Base<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_plus>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator-= (const Base<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_minus>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator%= (const Base<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_schur>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_elem1<eT,T1>::operator/= (const Base<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_subview_elem_inplace_div>(x);
  }



//
//



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::extract(Mat<eT>& actual_out, const subview_elem1<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_check_mixed<T1> tmp1(in.a.get_ref(), actual_out);
  const umat& aa = tmp1.M;
  
  arma_debug_check
    (
    ( aa.is_vec() == false ),
    "Mat::elem(): given object is not a vector"
    );
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  const Mat<eT>& m_local = in.m;
  
  const eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  const bool alias = (&actual_out == &m_local);
  
  arma_extra_debug_warn(alias, "subview_elem1::extract(): aliasing detected");
  
  Mat<eT>* tmp_out = alias ? new Mat<eT>() : 0;
  Mat<eT>& out     = alias ? *tmp_out      : actual_out;
  
  out.set_size(aa_n_elem, 1);
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j<aa_n_elem; i+=2, j+=2)
    {
    const uword ii = aa_mem[i];
    const uword jj = aa_mem[j];
    
    arma_debug_check( ( (ii >= m_n_elem) || (jj >= m_n_elem) ), "Mat::elem(): index out of bounds" );
    
    out_mem[i] = m_mem[ii];
    out_mem[j] = m_mem[jj];
    }
  
  if(i < aa_n_elem)
    {
    const uword ii = aa_mem[i];
    
    arma_debug_check( (ii >= m_n_elem) , "Mat::elem(): index out of bounds" );
    
    out_mem[i] = m_mem[ii];
    }
  
  if(alias == true)
    {
    actual_out = out;
    delete tmp_out;
    }
  }



template<typename eT, typename T1>
template<typename op_type>
inline
void
subview_elem1<eT,T1>::mat_inplace_op(Mat<eT>& out, const subview_elem1& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(in.a.get_ref());
  const umat& aa = tmp1.M;
  
  arma_debug_check
    (
    ( aa.is_vec() == false ),
    "Mat::elem(): given object is not a vector"
    );
  
  const uword* aa_mem    = aa.memptr();
  const uword  aa_n_elem = aa.n_elem;
  
  const unwrap_check< Mat<eT> > tmp2(in.m, out);
  const Mat<eT>& m_local      = tmp2.M;
  
  const eT*   m_mem    = m_local.memptr();
  const uword m_n_elem = m_local.n_elem;
  
  arma_debug_check( (out.n_elem != aa_n_elem), "Mat::elem(): size mismatch" );
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j<aa_n_elem; i+=2, j+=2)
    {
    const uword ii = aa_mem[i];
    const uword jj = aa_mem[j];
    
    arma_debug_check( ( (ii >= m_n_elem) || (jj >= m_n_elem) ), "Mat::elem(): index out of bounds" );
    
         if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { out_mem[i] += m_mem[ii]; out_mem[j] += m_mem[jj]; }
    else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { out_mem[i] -= m_mem[ii]; out_mem[j] -= m_mem[jj]; }
    else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { out_mem[i] *= m_mem[ii]; out_mem[j] *= m_mem[jj]; }
    else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { out_mem[i] /= m_mem[ii]; out_mem[j] /= m_mem[jj]; }
    }
  
  if(i < aa_n_elem)
    {
    const uword ii = aa_mem[i];
    
    arma_debug_check( (ii >= m_n_elem) , "Mat::elem(): index out of bounds" );
    
         if(is_same_type<op_type, op_subview_elem_inplace_plus >::value == true) { out_mem[i] += m_mem[ii]; }
    else if(is_same_type<op_type, op_subview_elem_inplace_minus>::value == true) { out_mem[i] -= m_mem[ii]; }
    else if(is_same_type<op_type, op_subview_elem_inplace_schur>::value == true) { out_mem[i] *= m_mem[ii]; }
    else if(is_same_type<op_type, op_subview_elem_inplace_div  >::value == true) { out_mem[i] /= m_mem[ii]; }
    }
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::plus_inplace(Mat<eT>& out, const subview_elem1& in)
  {
  arma_extra_debug_sigprint();
  
  mat_inplace_op<op_subview_elem_inplace_plus>(out, in);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::minus_inplace(Mat<eT>& out, const subview_elem1& in)
  {
  arma_extra_debug_sigprint();
  
  mat_inplace_op<op_subview_elem_inplace_minus>(out, in);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::schur_inplace(Mat<eT>& out, const subview_elem1& in)
  {
  arma_extra_debug_sigprint();
  
  mat_inplace_op<op_subview_elem_inplace_schur>(out, in);
  }



template<typename eT, typename T1>
inline
void
subview_elem1<eT,T1>::div_inplace(Mat<eT>& out, const subview_elem1& in)
  {
  arma_extra_debug_sigprint();
  
  mat_inplace_op<op_subview_elem_inplace_div>(out, in);
  }



//! @}
