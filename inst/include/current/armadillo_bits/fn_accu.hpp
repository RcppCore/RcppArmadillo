// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup fn_accu
//! @{



template<typename T1>
arma_warn_unused
inline
auto
accu(const T1& expr, const typename enable_if< is_arma_type<T1>::value >::result* junk = nullptr)
  {
  arma_debug_sigprint();
  
  arma_ignore(junk);
  
  return op_accu_mat::apply(expr);
  }



template<typename T1>
arma_warn_unused
inline
auto
accu(const T1& expr, const typename enable_if< is_arma_cube_type<T1>::value >::result* junk = nullptr)
  {
  arma_debug_sigprint();
  
  arma_ignore(junk);
  
  return op_accu_cube::apply(expr);
  }



template<typename T1>
arma_warn_unused
inline
auto
accu(const T1& expr, const typename enable_if< is_arma_sparse_type<T1>::value >::result* junk = nullptr)
  {
  arma_debug_sigprint();
  
  arma_ignore(junk);
  
  return spop_accu::apply(expr);
  }



template<typename T>
arma_warn_unused
inline
typename arma_scalar_only<T>::result
accu(const T& x)
  {
  return x;
  }



//! @}
