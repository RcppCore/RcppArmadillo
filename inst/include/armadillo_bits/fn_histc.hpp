// Copyright (C) 2015 Conrad Sanderson
// Copyright (C) 2015 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  is_not_complex<typename T1::elem_type>::value,
  const mtGlue<uword,T1,T2,glue_histc>
  >::result
histc
  (
  const Base<typename T1::elem_type,T1>& X,
  const Base<typename T1::elem_type,T2>& Y,
  const uword dim = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<uword,T1,T2,glue_histc>( X.get_ref(), Y.get_ref(), dim );
  }
