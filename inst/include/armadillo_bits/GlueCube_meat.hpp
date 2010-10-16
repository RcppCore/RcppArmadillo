// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup GlueCube
//! @{



template<typename T1, typename T2, typename glue_type>
inline
GlueCube<T1,T2,glue_type>::GlueCube(const BaseCube<typename T1::elem_type, T1>& in_A, const BaseCube<typename T1::elem_type, T2>& in_B)
  : A(in_A.get_ref())
  , B(in_B.get_ref())
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename glue_type>
inline
GlueCube<T1,T2,glue_type>::~GlueCube()
  {
  arma_extra_debug_sigprint();
  }



//! @}
