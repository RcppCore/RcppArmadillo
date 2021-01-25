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


//! \addtogroup fn_sort_index
//! @{



template<typename T1>
arma_warn_unused
arma_inline
const mtOp<uword,T1,op_sort_index>
sort_index
  (
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword,T1,op_sort_index>(X.get_ref(), uword(0), uword(0));
  }



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  ( (is_arma_type<T1>::value) && (is_same_type<T2, char>::value) ),
  const mtOp<uword,T1,op_sort_index>
  >::result
sort_index
  (
  const T1& X,
  const T2* sort_direction
  )
  {
  arma_extra_debug_sigprint();
  
  const char sig = (sort_direction != nullptr) ? sort_direction[0] : char(0);
  
  arma_debug_check( ((sig != 'a') && (sig != 'd')), "sort_index(): unknown sort direction" );
  
  return mtOp<uword,T1,op_sort_index>(X, ((sig == 'a') ? uword(0) : uword(1)), uword(0));
  }



//



template<typename T1>
arma_warn_unused
arma_inline
const mtOp<uword,T1,op_stable_sort_index>
stable_sort_index
  (
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<uword,T1,op_stable_sort_index>(X.get_ref(), uword(0), uword(0));
  }



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  ( (is_arma_type<T1>::value) && (is_same_type<T2, char>::value) ),
  const mtOp<uword,T1,op_stable_sort_index>
  >::result
stable_sort_index
  (
  const T1& X,
  const T2* sort_direction
  )
  {
  arma_extra_debug_sigprint();
  
  const char sig = (sort_direction != nullptr) ? sort_direction[0] : char(0);
  
  arma_debug_check( ((sig != 'a') && (sig != 'd')), "stable_sort_index(): unknown sort direction" );
  
  return mtOp<uword,T1,op_stable_sort_index>(X, ((sig == 'a') ? uword(0) : uword(1)), uword(0));
  }



//! @}
