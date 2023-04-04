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


#if defined(ARMA_USE_FFTW3)


namespace fftw3
  {
  template<typename eT>
  arma_inline
  void_ptr
  plan_dft_1d(int N, eT* input, eT* output, int fftw3_sign, unsigned int fftw3_flags)
    {
    arma_type_check((is_cx<eT>::value == false));
    
    if(is_cx_float<eT>::value)
      {
      return fftwf_plan_dft_1d(N, (cx_float*)input, (cx_float*)output, fftw3_sign, fftw3_flags);
      }
    else
    if(is_cx_double<eT>::value)
      {
      return fftw_plan_dft_1d(N, (cx_double*)input, (cx_double*)output, fftw3_sign, fftw3_flags);
      }
    
    return nullptr;
    }
  
  
  
  template<typename eT>
  arma_inline
  void
  execute(void_ptr plan)
    {
    arma_type_check((is_cx<eT>::value == false));
    
    if(is_cx_float<eT>::value)
      {
      fftwf_execute(plan);
      }
    else
    if(is_cx_double<eT>::value)
      {
      fftw_execute(plan);
      }
    }
  
  
  
  template<typename eT>
  arma_inline
  void
  destroy_plan(void_ptr plan)
    {
    arma_type_check((is_cx<eT>::value == false));
    
    if(is_cx_float<eT>::value)
      {
      fftwf_destroy_plan(plan);
      }
    else
    if(is_cx_double<eT>::value)
      {
      fftw_destroy_plan(plan);
      }
    }
  
  
  
  template<typename eT>
  arma_inline
  void
  cleanup()
    {
    arma_type_check((is_cx<eT>::value == false));
    
    if(is_cx_float<eT>::value)
      {
      fftwf_cleanup();
      }
    else
    if(is_cx_double<eT>::value)
      {
      fftw_cleanup();
      }
    }
  }


#endif
