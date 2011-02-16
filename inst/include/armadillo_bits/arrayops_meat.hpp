// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup arrayops
//! @{



template<typename eT>
arma_hot
inline
void
arrayops::copy(eT* dest, const eT* src, const u32 n_elem)
  {
  #if !defined(__OPTIMIZE__)
    {
    std::memcpy(dest, src, n_elem*sizeof(eT));
    }
  #else
    {
    switch(n_elem)
      {
      case 0:
        break;
      
      case 1:
        *dest = *src;
        break;
        
      case 2:
        dest[0] = src[0];
        dest[1] = src[1];
        break;
      
      case 3:
        dest[0] = src[0];
        dest[1] = src[1];
        dest[2] = src[2];
        break;
      
      case 4:
        dest[0] = src[0];
        dest[1] = src[1];
        dest[2] = src[2];
        dest[3] = src[3];
        break;
      
      default:
        if( n_elem <= (128/sizeof(eT)) )
          {
          u32 i,j;
          
          for(i=0, j=1; j<n_elem; i+=2, j+=2)
            {
            dest[i] = src[i];
            dest[j] = src[j];
            }
          
          if(i < n_elem)
            {
            dest[i] = src[i];
            }
          }
        else
          {
          std::memcpy(dest, src, n_elem*sizeof(eT));
          }
      }
    }
  #endif
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
arrayops::convert(out_eT* dest, const in_eT* src, const u32 n_elem)
  {
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    dest[i] = out_eT( src[i] );
    dest[j] = out_eT( src[j] );
    }
  
  if(i < n_elem)
    {
    dest[i] = out_eT( src[i] );
    }
  }



template<typename out_eT, typename in_eT>
arma_hot
inline
void
arrayops::convert_cx(out_eT* dest, const in_eT* src, const u32 n_elem)
  {
  u32 i,j;
  
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
arrayops::inplace_plus(eT* dest, const eT* src, const u32 n_elem)
  {
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    dest[i] += src[i];
    dest[j] += src[j];
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
arrayops::inplace_minus(eT* dest, const eT* src, const u32 n_elem)
  {
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    dest[i] -= src[i];
    dest[j] -= src[j];
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
arrayops::inplace_mul(eT* dest, const eT* src, const u32 n_elem)
  {
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    dest[i] *= src[i];
    dest[j] *= src[j];
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
arrayops::inplace_div(eT* dest, const eT* src, const u32 n_elem)
  {
  u32 i,j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    dest[i] /= src[i];
    dest[j] /= src[j];
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
arrayops::inplace_set(eT* dest, const eT val, const u32 n_elem)
  {
  u32 i,j;
  
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



template<typename eT>
arma_hot
inline
void
arrayops::inplace_plus(eT* dest, const eT val, const u32 n_elem)
  {
  u32 i,j;
  
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
arrayops::inplace_minus(eT* dest, const eT val, const u32 n_elem)
  {
  u32 i,j;
  
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
arrayops::inplace_mul(eT* dest, const eT val, const u32 n_elem)
  {
  u32 i,j;
  
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
arrayops::inplace_div(eT* dest, const eT val, const u32 n_elem)
  {
  u32 i,j;
  
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
arrayops::accumulate(const eT* src, const u32 n_elem)
  {
  u32 i,j;
  
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
arrayops::product(const eT* src, const u32 n_elem)
  {
  eT val1 = eT(1);
  eT val2 = eT(1);
  
  u32 i,j;
  
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
typename get_pod_type<eT>::result
arrayops::norm_1(const eT* src, const u32 n_elem)
  {
  typedef typename get_pod_type<eT>::result T;
  
  T acc = T(0);
  
  u32 i,j;
  
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



template<typename eT>
arma_hot
arma_pure
inline
eT
arrayops::norm_2(const eT* src, const u32 n_elem, const typename arma_not_cx<eT>::result* junk)
  {
  arma_ignore(junk);
  
  eT acc = eT(0);
  
  u32 i,j;
  
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



template<typename T>
arma_hot
arma_pure
inline
T
arrayops::norm_2(const std::complex<T>* src, const u32 n_elem)
  {
  T acc = T(0);
  
  u32 i,j;
  
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



template<typename eT>
arma_hot
arma_pure
inline
typename get_pod_type<eT>::result
arrayops::norm_k(const eT* src, const u32 n_elem, const int k)
  {
  typedef typename get_pod_type<eT>::result T;
  
  T acc = T(0);
  
  u32 i,j;
  
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



template<typename eT>
arma_hot
arma_pure
inline
typename get_pod_type<eT>::result
arrayops::norm_max(const eT* src, const u32 n_elem)
  {
  typedef typename get_pod_type<eT>::result T;
  
  T max_val = std::abs(src[0]);
  
  u32 i,j;
  
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



template<typename eT>
arma_hot
arma_pure
inline
typename get_pod_type<eT>::result
arrayops::norm_min(const eT* src, const u32 n_elem)
  {
  typedef typename get_pod_type<eT>::result T;
  
  T min_val = std::abs(src[0]);
  
  u32 i,j;
  
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
