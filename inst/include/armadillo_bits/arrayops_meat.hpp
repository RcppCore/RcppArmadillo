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


//! \addtogroup arrayops
//! @{



template<typename eT>
arma_hot
arma_inline
void
arrayops::copy(eT* dest, const eT* src, const uword n_elem)
  {
  if(is_cx<eT>::no)
    {
    if(n_elem <= 9)
      {
      arrayops::copy_small(dest, src, n_elem);
      }
    else
      {
      std::memcpy(dest, src, n_elem*sizeof(eT));
      }
    }
  else
    {
    if(n_elem > 0)  { std::memcpy(dest, src, n_elem*sizeof(eT)); }
    }
  }



template<typename eT>
arma_cold
inline
void
arrayops::copy_small(eT* dest, const eT* src, const uword n_elem)
  {
  switch(n_elem)
    {
    case  9:  dest[ 8] = src[ 8];
    // fallthrough
    case  8:  dest[ 7] = src[ 7];
    // fallthrough
    case  7:  dest[ 6] = src[ 6];
    // fallthrough
    case  6:  dest[ 5] = src[ 5];
    // fallthrough
    case  5:  dest[ 4] = src[ 4];
    // fallthrough
    case  4:  dest[ 3] = src[ 3];
    // fallthrough
    case  3:  dest[ 2] = src[ 2];
    // fallthrough
    case  2:  dest[ 1] = src[ 1];
    // fallthrough
    case  1:  dest[ 0] = src[ 0];
    // fallthrough
    default:  ;
    }
  }



template<typename eT>
inline
void
arrayops::fill_zeros(eT* dest, const uword n_elem)
  {
  typedef typename get_pod_type<eT>::result pod_type;
  
  if(std::numeric_limits<eT>::is_integer || std::numeric_limits<pod_type>::is_iec559)
    {
    if(n_elem > 0)  { std::memset((void*)dest, 0, sizeof(eT)*n_elem); }
    }
  else
    {
    arrayops::inplace_set_simple(dest, eT(0), n_elem);
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::replace(eT* mem, const uword n_elem, const eT old_val, const eT new_val)
  {
  if(arma_isnan(old_val))
    {
    for(uword i=0; i<n_elem; ++i)
      {
      eT& val = mem[i];
      
      val = (arma_isnan(val)) ? new_val : val;
      }
    }
  else
    {
    for(uword i=0; i<n_elem; ++i)
      {
      eT& val = mem[i];
      
      val = (val == old_val) ? new_val : val;
      }
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::clean(eT* mem, const uword n_elem, const eT abs_limit, const typename arma_not_cx<eT>::result* junk)
  {
  arma_ignore(junk);
  
  for(uword i=0; i<n_elem; ++i)
    {
    eT& val = mem[i];
    
    val = (std::abs(val) <= abs_limit) ? eT(0) : val;
    }
  }



template<typename T>
arma_hot
inline
void
arrayops::clean(std::complex<T>* mem, const uword n_elem, const T abs_limit)
  {
  typedef typename std::complex<T> eT;
  
  for(uword i=0; i<n_elem; ++i)
    {
    eT& val = mem[i];
    
    T val_real = std::real(val);
    T val_imag = std::imag(val);
    
    if(std::abs(val_real) <= abs_limit)
      {
      val_imag = (std::abs(val_imag) <= abs_limit) ? T(0) : val_imag;
      
      val = std::complex<T>(T(0), val_imag);
      }
    else
    if(std::abs(val_imag) <= abs_limit)
      {
      val = std::complex<T>(val_real, T(0));
      }
    }
  }



template<typename out_eT, typename in_eT>
arma_hot
arma_inline
void
arrayops::convert_cx_scalar
  (
        out_eT& out,
  const in_eT&  in,
  const typename arma_not_cx<out_eT>::result* junk1,
  const typename arma_not_cx< in_eT>::result* junk2
  )
  {
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  out = out_eT(in);
  }



template<typename out_eT, typename in_T>
arma_hot
arma_inline
void
arrayops::convert_cx_scalar
  (
        out_eT&             out,
  const std::complex<in_T>& in,
  const typename arma_not_cx<out_eT>::result* junk
  )
  {
  arma_ignore(junk);
  
  out = out_eT( in.real() );
  }



template<typename out_T, typename in_T>
arma_hot
arma_inline
void
arrayops::convert_cx_scalar
  (
        std::complex<out_T>& out,
  const std::complex< in_T>& in
  )
  {
  typedef std::complex<out_T> out_eT;
  
  out = out_eT(in);
  }



template<typename out_eT, typename in_eT>
arma_hot
inline
void
arrayops::convert(out_eT* dest, const in_eT* src, const uword n_elem)
  {
  if(is_same_type<out_eT,in_eT>::value)
    {
    const out_eT* src2 = (const out_eT*)src;
    
    if(dest != src2)  { arrayops::copy(dest, src2, n_elem); }
    
    return;
    }
  
  
  uword j;
  
  for(j=1; j<n_elem; j+=2)
    {
    const in_eT tmp_i = (*src);  src++;
    const in_eT tmp_j = (*src);  src++;
    
    // dest[i] = out_eT( tmp_i );
    // dest[j] = out_eT( tmp_j );
    
    (*dest) = (is_signed<out_eT>::value)
              ? out_eT( tmp_i )
              : ( cond_rel< is_signed<in_eT>::value >::lt(tmp_i, in_eT(0)) ? out_eT(0) : out_eT(tmp_i) );
    
    dest++;
    
    (*dest) = (is_signed<out_eT>::value)
              ? out_eT( tmp_j )
              : ( cond_rel< is_signed<in_eT>::value >::lt(tmp_j, in_eT(0)) ? out_eT(0) : out_eT(tmp_j) );
    dest++;
    }
  
  if((j-1) < n_elem)
    {
    const in_eT tmp_i = (*src);
    
    // dest[i] = out_eT( tmp_i );
    
    (*dest) = (is_signed<out_eT>::value)
              ? out_eT( tmp_i )
              : ( cond_rel< is_signed<in_eT>::value >::lt(tmp_i, in_eT(0)) ? out_eT(0) : out_eT(tmp_i) );
    }
  }



template<typename out_eT, typename in_eT>
arma_hot
inline
void
arrayops::convert_cx(out_eT* dest, const in_eT* src, const uword n_elem)
  {
  uword j;
  
  for(j=1; j<n_elem; j+=2)
    {
    arrayops::convert_cx_scalar( (*dest), (*src) );  dest++; src++;
    arrayops::convert_cx_scalar( (*dest), (*src) );  dest++; src++;
    }
  
  if((j-1) < n_elem)
    {
    arrayops::convert_cx_scalar( (*dest), (*src) );
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_plus(eT* dest, const eT* src, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      
      arrayops::inplace_plus_base(dest, src, n_elem);
      }
    else
      {
      arrayops::inplace_plus_base(dest, src, n_elem);
      }
    }
  else
    {
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      
      arrayops::inplace_plus_base(dest, src, n_elem);
      }
    else
      {
      arrayops::inplace_plus_base(dest, src, n_elem);
      }
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_minus(eT* dest, const eT* src, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      
      arrayops::inplace_minus_base(dest, src, n_elem);
      }
    else
      {
      arrayops::inplace_minus_base(dest, src, n_elem);
      }
    }
  else
    {
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      
      arrayops::inplace_minus_base(dest, src, n_elem);
      }
    else
      {
      arrayops::inplace_minus_base(dest, src, n_elem);
      }
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_mul(eT* dest, const eT* src, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      
      arrayops::inplace_mul_base(dest, src, n_elem);
      }
    else
      {
      arrayops::inplace_mul_base(dest, src, n_elem);
      }
    }
  else
    {
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      
      arrayops::inplace_mul_base(dest, src, n_elem);
      }
    else
      {
      arrayops::inplace_mul_base(dest, src, n_elem);
      }
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_div(eT* dest, const eT* src, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      
      arrayops::inplace_div_base(dest, src, n_elem);
      }
    else
      {
      arrayops::inplace_div_base(dest, src, n_elem);
      }
    }
  else
    {
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      
      arrayops::inplace_div_base(dest, src, n_elem);
      }
    else
      {
      arrayops::inplace_div_base(dest, src, n_elem);
      }
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_plus_base(eT* dest, const eT* src, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] += src[i];
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = src[i];
      const eT tmp_j = src[j];
      
      dest[i] += tmp_i;
      dest[j] += tmp_j;
      }
    
    if(i < n_elem)
      {
      dest[i] += src[i];
      }
    }
  #endif
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_minus_base(eT* dest, const eT* src, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] -= src[i];
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = src[i];
      const eT tmp_j = src[j];
      
      dest[i] -= tmp_i;
      dest[j] -= tmp_j;
      }
    
    if(i < n_elem)
      {
      dest[i] -= src[i];
      }
    }
  #endif
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_mul_base(eT* dest, const eT* src, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] *= src[i];
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = src[i];
      const eT tmp_j = src[j];
      
      dest[i] *= tmp_i;
      dest[j] *= tmp_j;
      }
    
    if(i < n_elem)
      {
      dest[i] *= src[i];
      }
    }
  #endif
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_div_base(eT* dest, const eT* src, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] /= src[i];
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = src[i];
      const eT tmp_j = src[j];
      
      dest[i] /= tmp_i;
      dest[j] /= tmp_j;
      }
    
    if(i < n_elem)
      {
      dest[i] /= src[i];
      }
    }
  #endif
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_set(eT* dest, const eT val, const uword n_elem)
  {
  if(val == eT(0))
    {
    arrayops::fill_zeros(dest, n_elem);
    }
  else
    {
    if( (n_elem <= 9) && (is_cx<eT>::no) )
      {
      arrayops::inplace_set_small(dest, val, n_elem);
      }
    else
      {
      arrayops::inplace_set_simple(dest, val, n_elem);
      }
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_set_simple(eT* dest, const eT val, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    arrayops::inplace_set_base(dest, val, n_elem);
    }
  else
    {
    arrayops::inplace_set_base(dest, val, n_elem);
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_set_base(eT* dest, const eT val, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] = val;
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      dest[i] = val;
      dest[j] = val;
      }
    
    if(i < n_elem)
      {
      dest[i] = val;
      }
    }
  #endif
  }



template<typename eT>
arma_cold
inline
void
arrayops::inplace_set_small(eT* dest, const eT val, const uword n_elem)
  {
  switch(n_elem)
    {
    case  9: dest[ 8] = val;
    // fallthrough
    case  8: dest[ 7] = val;
    // fallthrough
    case  7: dest[ 6] = val;
    // fallthrough
    case  6: dest[ 5] = val;
    // fallthrough
    case  5: dest[ 4] = val;
    // fallthrough
    case  4: dest[ 3] = val;
    // fallthrough
    case  3: dest[ 2] = val;
    // fallthrough
    case  2: dest[ 1] = val;
    // fallthrough
    case  1: dest[ 0] = val;
    // fallthrough
    default:;
    }
  }



template<typename eT, const uword n_elem>
arma_hot
inline
void
arrayops::inplace_set_fixed(eT* dest, const eT val)
  {
  for(uword i=0; i<n_elem; ++i)
    {
    dest[i] = val;
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_plus(eT* dest, const eT val, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    arrayops::inplace_plus_base(dest, val, n_elem);
    }
  else
    {
    arrayops::inplace_plus_base(dest, val, n_elem);
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_minus(eT* dest, const eT val, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    arrayops::inplace_minus_base(dest, val, n_elem);
    }
  else
    {
    arrayops::inplace_minus_base(dest, val, n_elem);
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_mul(eT* dest, const eT val, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    arrayops::inplace_mul_base(dest, val, n_elem);
    }
  else
    {
    arrayops::inplace_mul_base(dest, val, n_elem);
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_div(eT* dest, const eT val, const uword n_elem)
  {
  if(memory::is_aligned(dest))
    {
    memory::mark_as_aligned(dest);
    
    arrayops::inplace_div_base(dest, val, n_elem);
    }
  else
    {
    arrayops::inplace_div_base(dest, val, n_elem);
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_plus_base(eT* dest, const eT val, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] += val;
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      dest[i] += val;
      dest[j] += val;
      }
    
    if(i < n_elem)
      {
      dest[i] += val;
      }
    }
  #endif
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_minus_base(eT* dest, const eT val, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] -= val;
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      dest[i] -= val;
      dest[j] -= val;
      }
    
    if(i < n_elem)
      {
      dest[i] -= val;
      }
    }
  #endif
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_mul_base(eT* dest, const eT val, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] *= val;
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      dest[i] *= val;
      dest[j] *= val;
      }
    
    if(i < n_elem)
      {
      dest[i] *= val;
      }
    }
  #endif
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_div_base(eT* dest, const eT val, const uword n_elem)
  {
  #if defined(ARMA_SIMPLE_LOOPS)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      dest[i] /= val;
      }
    }
  #else
    {
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      dest[i] /= val;
      dest[j] /= val;
      }
    
    if(i < n_elem)
      {
      dest[i] /= val;
      }
    }
  #endif
  }



template<typename eT>
arma_hot
inline
eT
arrayops::accumulate(const eT* src, const uword n_elem)
  {
  #if defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)
    {
    eT acc = eT(0);
    
    if(memory::is_aligned(src))
      {
      memory::mark_as_aligned(src);
      for(uword i=0; i<n_elem; ++i)  { acc += src[i]; }
      }
    else
      {
      for(uword i=0; i<n_elem; ++i)  { acc += src[i]; }
      }
    
    return acc;
    }
  #else
    {
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    
    uword j;
    
    for(j=1; j<n_elem; j+=2)
      {
      acc1 += (*src);  src++;
      acc2 += (*src);  src++;
      }
    
    if((j-1) < n_elem)
      {
      acc1 += (*src);
      }
    
    return acc1 + acc2;
    }
  #endif
  }



template<typename eT>
arma_hot
inline
eT
arrayops::product(const eT* src, const uword n_elem)
  {
  eT val1 = eT(1);
  eT val2 = eT(1);
  
  uword i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    val1 *= src[i];
    val2 *= src[j];
    }
  
  if(i < n_elem)
    {
    val1 *= src[i];
    }
  
  return val1 * val2;
  }



template<typename eT>
arma_hot
inline
bool
arrayops::is_zero(const eT* mem, const uword n_elem, const eT abs_limit, const typename arma_not_cx<eT>::result* junk)
  {
  arma_ignore(junk);
  
  if(n_elem == 0)  { return false; }
  
  if(abs_limit == eT(0))
    {
    for(uword i=0; i<n_elem; ++i)
      {
      if(mem[i] != eT(0))  { return false; }
      }
    }
  else
    {
    for(uword i=0; i<n_elem; ++i)
      {
      if(std::abs(mem[i]) > abs_limit)  { return false; }
      }
    }
  
  return true;
  }



template<typename T>
arma_hot
inline
bool
arrayops::is_zero(const std::complex<T>* mem, const uword n_elem, const T abs_limit)
  {
  typedef typename std::complex<T> eT;
  
  if(n_elem == 0)  { return false; }
  
  if(abs_limit == T(0))
    {
    for(uword i=0; i<n_elem; ++i)
      {
      const eT& val = mem[i];
      
      if(std::real(val) != T(0))  { return false; }
      if(std::imag(val) != T(0))  { return false; }
      }
    }
  else
    {
    for(uword i=0; i<n_elem; ++i)
      {
      const eT& val = mem[i];
      
      if(std::abs(std::real(val)) > abs_limit)  { return false; }
      if(std::abs(std::imag(val)) > abs_limit)  { return false; }
      }
    }
  
  return true;
  }



template<typename eT>
arma_hot
inline
bool
arrayops::is_finite(const eT* src, const uword n_elem)
  {
  uword j;
  
  for(j=1; j<n_elem; j+=2)
    {
    const eT val_i = (*src);  src++;
    const eT val_j = (*src);  src++;
    
    if( (arma_isfinite(val_i) == false) || (arma_isfinite(val_j) == false) )
      {
      return false;
      }
    }
  
  if((j-1) < n_elem)
    {
    if(arma_isfinite(*src) == false)
      {
      return false;
      }
    }
  
  return true;
  }



template<typename eT>
arma_hot
inline
bool
arrayops::has_inf(const eT* src, const uword n_elem)
  {
  uword j;
  
  for(j=1; j<n_elem; j+=2)
    {
    const eT val_i = (*src);  src++;
    const eT val_j = (*src);  src++;
    
    if( arma_isinf(val_i) || arma_isinf(val_j) )  { return true; }
    }
  
  if((j-1) < n_elem)
    {
    if(arma_isinf(*src))  { return true; }
    }
  
  return false;
  }



template<typename eT>
arma_hot
inline
bool
arrayops::has_nan(const eT* src, const uword n_elem)
  {
  uword j;
  
  for(j=1; j<n_elem; j+=2)
    {
    const eT val_i = (*src);  src++;
    const eT val_j = (*src);  src++;
    
    if( arma_isnan(val_i) || arma_isnan(val_j) )  { return true; }
    }
  
  if((j-1) < n_elem)
    {
    if(arma_isnan(*src))  { return true; }
    }
  
  return false;
  }



//! @}
