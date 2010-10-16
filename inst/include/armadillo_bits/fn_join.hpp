// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_join
//! @{



template<typename T1, typename T2>
inline
const Glue<T1, T2, glue_join>
join_cols(const Base<typename T1::elem_type,T1>& A, const Base<typename T1::elem_type,T2>& B)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_join>(A.get_ref(), B.get_ref(), 0);
  }



template<typename T1, typename T2>
inline
const Glue<T1, T2, glue_join>
join_rows(const Base<typename T1::elem_type,T1>& A, const Base<typename T1::elem_type,T2>& B)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_join>(A.get_ref(), B.get_ref(), 1);
  }



template<typename T1, typename T2>
inline
const GlueCube<T1, T2, glue_join>
join_slices(const BaseCube<typename T1::elem_type,T1>& A, const BaseCube<typename T1::elem_type,T2>& B)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_join>(A.get_ref(), B.get_ref());
  }



//! @}
