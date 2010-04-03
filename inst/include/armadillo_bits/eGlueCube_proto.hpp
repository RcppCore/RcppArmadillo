// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup eGlueCube
//! @{


template<typename T1, typename T2, typename eglue_type>
class eGlueCube : public BaseCube<typename T1::elem_type, eGlueCube<T1, T2, eglue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  const ProxyCube<T1> P1;
  const ProxyCube<T2> P2;
  
  arma_inline ~eGlueCube();
  arma_inline  eGlueCube(const T1& in_A, const T2& in_B);
  
  };



//! @}
