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


//! \addtogroup spop_sum
//! @{


class spop_sum
  {
  public:
  
  template<typename T1>
  arma_hot inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_sum>& in);
  
  template<typename T1>
  arma_hot inline static void apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& p, const uword dim);
  };


//! @}
