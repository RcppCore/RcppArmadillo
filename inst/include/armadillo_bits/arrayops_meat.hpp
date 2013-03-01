// Copyright (C) 2011-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2011-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup arrayops
//! @{



template<typename eT>
arma_hot
arma_inline
void
arrayops::copy(eT* dest, const eT* src, const uword n_elem)
  {
  switch(n_elem)
    {
    default:
      arrayops::copy_big(dest, src, n_elem);
      break;
    case 8:
      dest[7] = src[7];
    case 7:
      dest[6] = src[6];
    case 6:
      dest[5] = src[5];
    case 5:
      dest[4] = src[4];
    case 4:
      dest[3] = src[3];
    case 3:
      dest[2] = src[2];
    case 2:
      dest[1] = src[1];
    case 1:
      dest[0] = src[0];
    }
  }



template<typename eT>
inline
void
arrayops::copy_big(eT* dest, const eT* src, const uword n_elem)
  {
  switch(n_elem)
    {
    default:
      std::memcpy(dest, src, n_elem*sizeof(eT));
      break;
    case 32:
      dest[31] = src[31];
    case 31:
      dest[30] = src[30];
    case 30:
      dest[29] = src[29];
    case 29:
      dest[28] = src[28];
    case 28:
      dest[27] = src[27];
    case 27:
      dest[26] = src[26];
    case 26:
      dest[25] = src[25];
    case 25:
      dest[24] = src[24];
    case 24:
      dest[23] = src[23];
    case 23:
      dest[22] = src[22];
    case 22:
      dest[21] = src[21];
    case 21:
      dest[20] = src[20];
    case 20:
      dest[19] = src[19];
    case 19:
      dest[18] = src[18];
    case 18:
      dest[17] = src[17];
    case 17:
      dest[16] = src[16];
    case 16:
      dest[15] = src[15];
    case 15:
      dest[14] = src[14];
    case 14:
      dest[13] = src[13];
    case 13:
      dest[12] = src[12];
    case 12:
      dest[11] = src[11];
    case 11:
      dest[10] = src[10];
    case 10:
      dest[9] = src[9];
    case 9:
      dest[8] = src[8];
    case 8:
      dest[7] = src[7];
    case 7:
      dest[6] = src[6];
    case 6:
      dest[5] = src[5];
    case 5:
      dest[4] = src[4];
    case 4:
      dest[3] = src[3];
    case 3:
      dest[2] = src[2];
    case 2:
      dest[1] = src[1];
    case 1:
      dest[0] = src[0];
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::copy_forwards(eT* dest, const eT* src, const uword n_elem)
  {
  // can't use std::memcpy(), as we don't know how it copies data
  uword i,j;
  
  for(i=0, j=1; j < n_elem; i+=2, j+=2)
    {
    dest[i] = src[i];
    dest[j] = src[j];
    }
  
  if(i < n_elem)
    {
    dest[i] = src[i];
    }
  }



// template<typename eT>
// arma_hot
// inline
// void
// arrayops::copy_backwards(eT* dest, const eT* src, const uword n_elem)
//   {
//   // can't use std::memcpy(), as we don't know how it copies data
//   
//   switch(n_elem)
//     {
//     default:
//       for(uword i = (n_elem-1); i >= 1; --i) { dest[i] = src[i]; }
//       // NOTE: the 'break' statement has been deliberately omitted
//     case 1:
//       dest[0] = src[0];
//     case 0:
//       ;
//     }
//   }



template<typename eT>
arma_hot
inline
void
arrayops::copy_backwards(eT* dest, const eT* src, const uword n_elem)
  {
  // can't use std::memcpy(), as we don't know how it copies data
  
  switch(n_elem)
    {
    default:
      {
      uword i, j;
      for(i = (n_elem-1), j = (n_elem-2); j >= 2; i-=2, j-=2)
        {
        const eT tmp_i = src[i];
        const eT tmp_j = src[j];
        
        dest[i] = tmp_i;
        dest[j] = tmp_j;
        }
      
      // j is less than 2:  it can be 1 or 0
      // i is j+1, ie. less than 3: it can be 2 or 1
      
      if(i == 2)
        {
        dest[2] = src[2];
        }
      }
      // NOTE: the 'break' statement has been deliberately omitted
    case 2:
      dest[1] = src[1];
    case 1:
      dest[0] = src[0];
    case 0:
      ;
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
  uword i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    const in_eT tmp_i = src[i];
    const in_eT tmp_j = src[j];
    
    // dest[i] = out_eT( tmp_i );
    // dest[j] = out_eT( tmp_j );
    
    dest[i] = (is_signed<out_eT>::value)
              ? out_eT( tmp_i )
              : ( cond_rel< is_signed<in_eT>::value >::lt(tmp_i, in_eT(0)) ? out_eT(0) : out_eT(tmp_i) );
              
    dest[j] = (is_signed<out_eT>::value)
              ? out_eT( tmp_j )
              : ( cond_rel< is_signed<in_eT>::value >::lt(tmp_j, in_eT(0)) ? out_eT(0) : out_eT(tmp_j) );
    }
  
  if(i < n_elem)
    {
    const in_eT tmp_i = src[i];
    
    // dest[i] = out_eT( tmp_i );
    
    dest[i] = (is_signed<out_eT>::value)
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
  uword i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    arrayops::convert_cx_scalar( dest[i], src[i] );
    arrayops::convert_cx_scalar( dest[j], src[j] );
    }
  
  if(i < n_elem)
    {
    arrayops::convert_cx_scalar( dest[i], src[i] );
    }
  }



template<typename eT>
arma_hot
inline
void
arrayops::inplace_plus(eT* dest, const eT* src, const uword n_elem)
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



template<typename eT>
arma_hot
inline
void
arrayops::inplace_minus(eT* dest, const eT* src, const uword n_elem)
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



template<typename eT>
arma_hot
inline
void
arrayops::inplace_mul(eT* dest, const eT* src, const uword n_elem)
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



template<typename eT>
arma_hot
inline
void
arrayops::inplace_div(eT* dest, const eT* src, const uword n_elem)
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



template<typename eT>
arma_hot
inline
void
arrayops::inplace_set(eT* dest, const eT val, const uword n_elem)
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



template<typename eT>
arma_hot
inline
void
arrayops::inplace_minus(eT* dest, const eT val, const uword n_elem)
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



template<typename eT>
arma_hot
inline
void
arrayops::inplace_mul(eT* dest, const eT val, const uword n_elem)
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



template<typename eT>
arma_hot
inline
void
arrayops::inplace_div(eT* dest, const eT val, const uword n_elem)
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



template<typename eT>
arma_hot
arma_pure
inline
eT
arrayops::accumulate(const eT* src, const uword n_elem)
  {
  uword i,j;
  
  eT acc1 = eT(0);
  eT acc2 = eT(0);
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    acc1 += src[i];
    acc2 += src[j];
    }
  
  if(i < n_elem)
    {
    acc1 += src[i];
    }
  
  return acc1 + acc2;
  }



template<typename eT>
arma_hot
arma_pure
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
arma_pure
inline
bool
arrayops::is_finite(const eT* src, const uword n_elem)
  {
  uword i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    const eT val_i = src[i];
    const eT val_j = src[j];
    
    if( (arma_isfinite(val_i) == false) || (arma_isfinite(val_j) == false) )
      {
      return false;
      }
    }
  
  if(i < n_elem)
    {
    if(arma_isfinite(src[i]) == false)
      {
      return false;
      }
    }
  
  return true;
  }



// TODO: this function is currently not used
template<typename eT>
arma_hot
arma_pure
inline
typename get_pod_type<eT>::result
arrayops::norm_1(const eT* src, const uword n_elem)
  {
  typedef typename get_pod_type<eT>::result T;
  
  T acc = T(0);
  
  uword i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    acc += std::abs(src[i]);
    acc += std::abs(src[j]);
    }
  
  if(i < n_elem)
    {
    acc += std::abs(src[i]);
    }
  
  return acc;
  }



// TODO: this function is currently not used
template<typename eT>
arma_hot
arma_pure
inline
eT
arrayops::norm_2(const eT* src, const uword n_elem, const typename arma_not_cx<eT>::result* junk)
  {
  arma_ignore(junk);
  
  eT acc = eT(0);
  
  uword i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    const eT tmp_i = src[i];
    const eT tmp_j = src[j];
    
    acc += tmp_i * tmp_i;
    acc += tmp_j * tmp_j;
    }
  
  if(i < n_elem)
    {
    const eT tmp_i = src[i];
    
    acc += tmp_i * tmp_i;
    }
  
  return std::sqrt(acc);
  }



// TODO: this function is currently not used
template<typename T>
arma_hot
arma_pure
inline
T
arrayops::norm_2(const std::complex<T>* src, const uword n_elem)
  {
  T acc = T(0);
  
  uword i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    const T tmp_i = std::abs(src[i]);
    const T tmp_j = std::abs(src[j]);
    
    acc += tmp_i * tmp_i;
    acc += tmp_j * tmp_j;
    }
  
  if(i < n_elem)
    {
    const T tmp_i = std::abs(src[i]);
    
    acc += tmp_i * tmp_i;
    }
  
  return std::sqrt(acc);
  }



// TODO: this function is currently not used
template<typename eT>
arma_hot
arma_pure
inline
typename get_pod_type<eT>::result
arrayops::norm_k(const eT* src, const uword n_elem, const int k)
  {
  typedef typename get_pod_type<eT>::result T;
  
  T acc = T(0);
  
  uword i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    acc += std::pow(std::abs(src[i]), k);
    acc += std::pow(std::abs(src[j]), k);
    }
  
  if(i < n_elem)
    {
    acc += std::pow(std::abs(src[i]), k);
    }
  
  return std::pow(acc, T(1)/T(k));
  }



// TODO: this function is currently not used
template<typename eT>
arma_hot
arma_pure
inline
typename get_pod_type<eT>::result
arrayops::norm_max(const eT* src, const uword n_elem)
  {
  typedef typename get_pod_type<eT>::result T;
  
  T max_val = std::abs(src[0]);
  
  uword i,j;
  
  for(i=1, j=2; j<n_elem; i+=2, j+=2)
    {
    const T tmp_i = std::abs(src[i]);
    const T tmp_j = std::abs(src[j]);
    
    if(max_val < tmp_i) { max_val = tmp_i; }
    if(max_val < tmp_j) { max_val = tmp_j; }
    }
  
  if(i < n_elem)
    {
    const T tmp_i = std::abs(src[i]);
    
    if(max_val < tmp_i) { max_val = tmp_i; }
    }
  
  return max_val;
  }



// TODO: this function is currently not used
template<typename eT>
arma_hot
arma_pure
inline
typename get_pod_type<eT>::result
arrayops::norm_min(const eT* src, const uword n_elem)
  {
  typedef typename get_pod_type<eT>::result T;
  
  T min_val = std::abs(src[0]);
  
  uword i,j;
  
  for(i=1, j=2; j<n_elem; i+=2, j+=2)
    {
    const T tmp_i = std::abs(src[i]);
    const T tmp_j = std::abs(src[j]);
    
    if(min_val > tmp_i) { min_val = tmp_i; }
    if(min_val > tmp_j) { min_val = tmp_j; }
    }
  
  if(i < n_elem)
    {
    const T tmp_i = std::abs(src[i]);
    
    if(min_val > tmp_i) { min_val = tmp_i; }
    }
  
  return min_val;
  }



//! @}
