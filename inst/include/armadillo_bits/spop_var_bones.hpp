// Copyright (C) 2012 Ryan Curtin
// Copyright (C) 2012 Conrad Sanderson
//
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup spop_var
//! @{



//! Class for finding variance values of a sparse matrix
class spop_var
  {
  public:

  template<typename T1>
  inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_var>& in);

  template<typename T1>
  inline static void apply_noalias(SpMat<typename T1::pod_type>& out, const SpProxy<T1>& p, const uword norm_type, const uword dim);
  
  // Calculate variance of a sparse vector, where we can directly use the memory.
  template<typename T1>
  inline static typename T1::pod_type var_vec(const T1& X, const uword norm_type = 0);
  
  // Calculate the variance directly.  Because this is for sparse matrices, we
  // specify both the number of elements in the array (the length of the array)
  // as well as the actual number of elements when zeros are included.
  template<typename eT>
  inline static eT direct_var(const eT* const X, const uword length, const uword N, const uword norm_type = 0);

  // For complex numbers.

  template<typename T>
  inline static T direct_var(const std::complex<T>* const X, const uword length, const uword N, const uword norm_type = 0);

  // Calculate the variance using iterators, for non-complex numbers.
  template<typename T1, typename eT>
  inline static eT iterator_var(T1& it, const T1& end, const uword n_zero, const uword norm_type, const eT junk1, const typename arma_not_cx<eT>::result* junk2 = 0);

  // Calculate the variance using iterators, for complex numbers.
  template<typename T1, typename eT>
  inline static typename get_pod_type<eT>::result iterator_var(T1& it, const T1& end, const uword n_zero, const uword norm_type, const eT junk1, const typename arma_cx_only<eT>::result* junk2 = 0);

  };



//! @}
  
