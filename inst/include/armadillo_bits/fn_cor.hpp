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


//! \addtogroup fn_cor
//! @{



template<typename T1>
arma_warn_unused
inline
const Op<T1, op_cor>
cor(const Base<typename T1::elem_type,T1>& X, const uword norm_type = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (norm_type > 1), "cor(): parameter 'norm_type' must be 0 or 1" );
  
  return Op<T1, op_cor>(X.get_ref(), norm_type, 0);
  }



template<typename T1, typename T2>
arma_warn_unused
inline
const Glue<T1,T2,glue_cor>
cor(const Base<typename T1::elem_type,T1>& A, const Base<typename T1::elem_type,T2>& B, const uword norm_type = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (norm_type > 1), "cor(): parameter 'norm_type' must be 0 or 1" );
  
  return Glue<T1, T2, glue_cor>(A.get_ref(), B.get_ref(), norm_type);
  }



//! @}
