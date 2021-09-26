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


//! \addtogroup op_dotext
//! @{



class op_dotext
  : public traits_op_default
  {
  public:
  
  
  template<typename eT>
  inline static eT direct_rowvec_mat_colvec       (const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_transmat_colvec  (const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_diagmat_colvec   (const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_invdiagmat_colvec(const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_invdiagvec_colvec(const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  };



//! @}

