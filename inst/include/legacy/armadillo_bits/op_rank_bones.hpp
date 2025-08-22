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



//! \addtogroup op_rank
//! @{



class op_rank
  : public traits_op_default
  {
  public:
  
  template<typename T1> inline static bool apply(uword& out, const Base<typename T1::elem_type,T1>& expr, const typename T1::pod_type tol);
  
  template<typename eT> inline static bool  apply_gen(uword& out, Mat<eT>& A, typename get_pod_type<eT>::result tol);
  
  template<typename eT> inline static bool  apply_sym(uword& out, Mat<eT>& A, typename get_pod_type<eT>::result tol);
  
  template<typename eT> inline static bool apply_diag(uword& out, Mat<eT>& A, typename get_pod_type<eT>::result tol);
  };



//! @}
