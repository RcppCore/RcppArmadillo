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


//! \addtogroup spglue_minus
//! @{



class spglue_minus
  {
  public:
  
  template<typename T1, typename T2>
  arma_hot inline static void apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_minus>& X);
  
  template<typename eT, typename T1, typename T2>
  arma_hot inline static void apply_noalias(SpMat<eT>& result, const SpProxy<T1>& pa, const SpProxy<T2>& pb);
  };



class spglue_minus2
  {
  public:
  
  template<typename T1, typename T2>
  arma_hot inline static void apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_minus2>& X);
  };



//! @}

