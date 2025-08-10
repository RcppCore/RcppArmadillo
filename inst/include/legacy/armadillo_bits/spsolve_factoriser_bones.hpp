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


//! \addtogroup spsolve_factoriser
//! @{



class spsolve_factoriser
  {
  private:
  
  void_ptr worker_ptr          = nullptr;
  uword    elem_type_indicator = 0;
  uword    n_rows              = 0;
  double   rcond_value         = double(0);
  
  template<typename worker_type> inline void delete_worker();
  
  inline void cleanup();
  
  
  public:
  
  inline ~spsolve_factoriser();
  inline  spsolve_factoriser();
  
  inline void reset();
  
  inline double rcond() const;
  
  template<typename T1> inline bool factorise(const SpBase<typename T1::elem_type,T1>& A_expr, const spsolve_opts_base& settings = spsolve_opts_none(), const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr);
  
  template<typename T1> inline bool solve(Mat<typename T1::elem_type>& X, const Base<typename T1::elem_type,T1>& B_expr, const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr);
  
  inline      spsolve_factoriser(const spsolve_factoriser&) = delete;
  inline void operator=         (const spsolve_factoriser&) = delete;
  };



//! @}
