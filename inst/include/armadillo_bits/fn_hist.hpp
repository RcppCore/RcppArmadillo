// Copyright (C) 2012 NICTA (www.nicta.com.au)
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



template<typename T1>
inline
const mtOp<uword,T1,op_hist>
hist
  (
  const Base<typename T1::elem_type,T1>& A,
  const uword n_bins = 10,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename arma_not_cx<typename T1::elem_type>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return mtOp<uword,T1,op_hist>( A.get_ref(), n_bins, 0 );
  }



template<typename T1, typename T2>
inline
const mtGlue<uword,T1,T2,glue_hist>
hist
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const uword dim = 0,
  const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return mtGlue<uword,T1,T2,glue_hist>( A.get_ref(), B.get_ref(), dim );
  }
