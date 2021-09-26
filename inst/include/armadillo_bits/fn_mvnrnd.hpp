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


//! \addtogroup fn_mvnrnd
//! @{



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  const Glue<T1, T2, glue_mvnrnd_vec>
  >::result
mvnrnd(const Base<typename T1::elem_type, T1>& M, const Base<typename T1::elem_type, T2>& C)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_mvnrnd_vec>(M.get_ref(), C.get_ref());
  }



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  const Glue<T1, T2, glue_mvnrnd>
  >::result
mvnrnd(const Base<typename T1::elem_type, T1>& M, const Base<typename T1::elem_type, T2>& C, const uword N)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_mvnrnd>(M.get_ref(), C.get_ref(), N);
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  bool
  >::result
mvnrnd(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& M, const Base<typename T1::elem_type, T2>& C)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_mvnrnd::apply_direct(out, M.get_ref(), C.get_ref(), uword(1));
  
  if(status == false)
    {
    out.soft_reset();
    arma_debug_warn_level(3, "mvnrnd(): given covariance matrix is not symmetric positive semi-definite");
    }
  
  return status;
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  bool
  >::result
mvnrnd(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& M, const Base<typename T1::elem_type, T2>& C, const uword N)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_mvnrnd::apply_direct(out, M.get_ref(), C.get_ref(), N);
  
  if(status == false)
    {
    out.soft_reset();
    arma_debug_warn_level(3, "mvnrnd(): given covariance matrix is not symmetric positive semi-definite");
    }
  
  return status;
  }



//! @}
