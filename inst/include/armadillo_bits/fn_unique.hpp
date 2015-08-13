// Copyright (C) 2012-2015 Conrad Sanderson
// Copyright (C) 2012-2015 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Arnold Wiliem
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



template<typename T1>
inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1,op_unique>
  >::result
unique(const T1& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1,op_unique>(A);
  }
