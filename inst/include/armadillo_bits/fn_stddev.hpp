// Copyright (C) 2009-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_stddev
//! @{



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_stddev>
stddev
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword norm_type = 0,
  const uword dim = 0,
  const typename enable_if<resolves_to_vector<T1>::value == false>::result* junk1 = 0,
  const typename enable_if<is_basevec<T1>::value == false>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return mtOp<typename T1::pod_type, T1, op_stddev>(X.get_ref(), norm_type, dim);
  }



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_stddev>
stddev
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword norm_type,
  const uword dim,
  const typename enable_if<resolves_to_vector<T1>::value == true>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return mtOp<typename T1::pod_type, T1, op_stddev>(X.get_ref(), norm_type, dim);
  }



template<typename T1>
inline
arma_warn_unused
typename T1::pod_type
stddev
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword norm_type = 0,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename enable_if<resolves_to_vector<T1>::value == true>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return std::sqrt( op_var::var_vec(X.get_ref(), norm_type) );
  }



template<typename T1>
inline
arma_warn_unused
typename T1::pod_type
stddev
  (
  const T1& X,
  const uword norm_type = 0,
  const arma_empty_class junk1 = arma_empty_class(),
  const typename enable_if<is_basevec<T1>::value == true>::result* junk2 = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  return std::sqrt( op_var::var_vec(X, norm_type) );
  }



//! @}
