// Copyright (C) 2012 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Conrad Sanderson
// Copyright (C) 2012 Arnold Wiliem
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
const Op<T1,op_unique>
unique
  (
  const Base<typename T1::elem_type,T1>& A,
  const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1,op_unique>( A.get_ref() );
  }
