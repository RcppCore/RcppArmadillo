// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------



//! \addtogroup fn_shift
//! @{


template<typename T1>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value && resolves_to_vector<T1>::yes,
  const Op<T1, op_shift_vec>
  >::result
shift
  (
  const T1&   X,
  const sword N
  )
  {
  arma_extra_debug_sigprint();
  
  const uword len = (N < 0) ? uword(-N) : uword(N);
  const uword neg = (N < 0) ? uword( 1) : uword(0);
  
  return Op<T1, op_shift_vec>(X, len, neg);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_type<T1>::value && resolves_to_vector<T1>::no,
  Mat<typename T1::elem_type>
  >::result
shift
  (
  const T1&   X,
  const sword N
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword len = (N < 0) ? uword(-N) : uword(N);
  const uword neg = (N < 0) ? uword( 1) : uword(0);
  
  quasi_unwrap<T1> U(X);
  
  Mat<eT> out;
  
  op_shift::apply_noalias(out, U.M, len, neg, 0);
  
  return out;
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value),
  Mat<typename T1::elem_type>
  >::result
shift
  (
  const T1&   X,
  const sword N,
  const uword dim
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (dim > 1), "shift(): parameter 'dim' must be 0 or 1" );
  
  const uword len = (N < 0) ? uword(-N) : uword(N);
  const uword neg = (N < 0) ? uword( 1) : uword(0);
  
  quasi_unwrap<T1> U(X);
  
  Mat<eT> out;
  
  op_shift::apply_noalias(out, U.M, len, neg, dim);
  
  return out;
  }



//



template<typename T1>
arma_warn_unused
inline
SpMat<typename T1::elem_type>
shift
  (
  const SpBase<typename T1::elem_type,T1>& expr,
  const sword N,
  const uword dim = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (dim > 1), "shift(): parameter 'dim' must be 0 or 1" );
  
  const uword len = (N < 0) ? uword(-N) : uword(N);
  const uword neg = (N < 0) ? uword( 1) : uword(0);
  
  unwrap_spmat<T1> U(expr.get_ref());
  
  SpMat<eT> out;
  
  spop_shift::apply_noalias(out, U.M, len, neg, dim);
  
  return out;
  }



//! @}
