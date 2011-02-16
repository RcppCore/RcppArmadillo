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


class arrayops
  {
  public:
  
  template<typename eT>
  arma_hot inline static void
  copy(eT* dest, const eT* src, const u32 n_elem);
  
  
  // 
  // array = convert(array)
  
  template<typename out_eT, typename in_eT>
  arma_hot arma_inline static void
  convert_cx_scalar(out_eT& out, const in_eT&  in, const typename arma_not_cx<out_eT>::result* junk1 = 0, const typename arma_not_cx< in_eT>::result* junk2 = 0);
  
  template<typename out_eT, typename in_T>
  arma_hot arma_inline static void
  convert_cx_scalar(out_eT& out, const std::complex<in_T>& in, const typename arma_not_cx<out_eT>::result* junk = 0);
  
  template<typename out_T, typename in_T>
  arma_hot arma_inline static void
  convert_cx_scalar(std::complex<out_T>& out, const std::complex< in_T>& in);
  
  template<typename out_eT, typename in_eT>
  arma_hot inline static void
  convert(out_eT* dest, const in_eT* src, const u32 n_elem);
  
  template<typename out_eT, typename in_eT>
  arma_hot inline static void
  convert_cx(out_eT* dest, const in_eT* src, const u32 n_elem);
  
  
  // 
  // array op= array
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_plus(eT* dest, const eT* src, const u32 n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_minus(eT* dest, const eT* src, const u32 n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_mul(eT* dest, const eT* src, const u32 n_elem);
   
  template<typename eT>
  arma_hot inline static
  void
  inplace_div(eT* dest, const eT* src, const u32 n_elem);
  
  
  // 
  // array op= scalar
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_set(eT* dest, const eT val, const u32 n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_plus(eT* dest, const eT val, const u32 n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_minus(eT* dest, const eT val, const u32 n_elem);
  
  template<typename eT>
  arma_hot inline static void
  inplace_mul(eT* dest, const eT val, const u32 n_elem);
  
  template<typename eT>
  arma_hot inline static
  void
  inplace_div(eT* dest, const eT val, const u32 n_elem);
  
  
  // 
  // scalar = op(array)
  
  template<typename eT>
  arma_hot arma_pure inline static
  eT
  accumulate(const eT* src, const u32 n_elem);
  
  template<typename eT>
  arma_hot arma_pure inline static
  eT
  product(const eT* src, const u32 n_elem);
  
  template<typename eT>
  arma_hot arma_pure inline static
  typename get_pod_type<eT>::result
  norm_1(const eT* src, const u32 n_elem);
  
  template<typename eT>
  arma_hot arma_pure inline static
  eT
  norm_2(const eT* src, const u32 n_elem, const typename arma_not_cx<eT>::result* junk = 0);
  
  template<typename T>
  arma_hot arma_pure inline static
  T
  norm_2(const std::complex<T>* src, const u32 n_elem);
  
  template<typename eT>
  arma_hot arma_pure inline static
  typename get_pod_type<eT>::result
  norm_k(const eT* src, const u32 n_elem, const int k);
  
  template<typename eT>
  arma_hot arma_pure inline static
  typename get_pod_type<eT>::result
  norm_max(const eT* src, const u32 n_elem);
  
  template<typename eT>
  arma_hot arma_pure inline static
  typename get_pod_type<eT>::result
  norm_min(const eT* src, const u32 n_elem);
  
  };



//! @}
