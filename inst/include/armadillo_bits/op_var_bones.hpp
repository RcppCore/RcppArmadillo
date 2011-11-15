// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
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
  inline static eT direct_var(const eT* const X, const uword N, const uword norm_type = 0);
  
  template<typename T>
  inline static  T direct_var(const std::complex<T>* const X, const uword N, const uword norm_type = 0);
  
  
  template<typename eT>
  inline static typename get_pod_type<eT>::result direct_var(const subview_row<eT>& X, const uword norm_type = 0);
  
  template<typename eT>
  inline static typename get_pod_type<eT>::result direct_var(const subview_col<eT>& X, const uword norm_type = 0);
  
  template<typename eT>
  inline static typename get_pod_type<eT>::result direct_var(const diagview<eT>& X, const uword norm_type = 0);
  
  
  template<typename T1>
  inline static void apply(Mat<typename T1::pod_type>& out, const mtOp<typename T1::pod_type, T1, op_var>& in);
  
  
  template<typename eT>
  inline static eT direct_var_robust(const eT* const X, const uword N, const uword norm_type = 0);
  
  template<typename T>
  inline static  T direct_var_robust(const std::complex<T>* const X, const uword N, const uword norm_type = 0);
  };



//! @}
