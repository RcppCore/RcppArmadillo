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


extern "C"
  {
  // function prefix for single precision: fftwf_
  // function prefix for double precision: fftw_
  
  
  // single precision (float)
  
  void_ptr fftwf_plan_dft_1d(int N, void* input, void* output, int fftw3_sign, unsigned int fftw3_flags);
  
  void      fftwf_execute(void_ptr plan);
  void fftwf_destroy_plan(void_ptr plan);
  
  void fftwf_cleanup();
  
  
  // double precision (double)
  
  void_ptr fftw_plan_dft_1d(int N, void* input, void* output, int fftw3_sign, unsigned int fftw3_flags);
  
  void      fftw_execute(void_ptr plan);
  void fftw_destroy_plan(void_ptr plan);
  
  void fftw_cleanup();
  }


#endif
