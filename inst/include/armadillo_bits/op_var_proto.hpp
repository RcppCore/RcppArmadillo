// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_var
//! @{



//! Class for finding variance values of a matrix
class op_var
  {
  public:
  
  template<typename eT>
  inline static eT direct_var(const eT* const X, const u32 N, const u32 norm_type = 0);
  
  template<typename T>
  inline static T direct_var(const std::complex<T>* const X, const u32 N, const u32 norm_type = 0);
  
  
  template<typename eT>
  inline static typename get_pod_type<eT>::result direct_var(const subview_row<eT>& X, const u32 norm_type = 0);
  
  template<typename eT>
  inline static typename get_pod_type<eT>::result direct_var(const subview_col<eT>& X, const u32 norm_type = 0);
  
  template<typename eT>
  inline static typename get_pod_type<eT>::result direct_var(const diagview<eT>& X, const u32 norm_type = 0);
  
  
  template<typename eT>
  inline static void apply(Mat< typename get_pod_type<eT>::result >& out, const Mat<eT>& X, const u32 norm_type, const u32 dim);
  
  };



//! @}
