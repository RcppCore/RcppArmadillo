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


#if defined(ARMA_USE_FFTW3) && !defined(FFTW3_H)


// prefix for single precision: fftwf_
// prefix for double precision: fftw_


typedef void fftwf_complex;
typedef void  fftw_complex;

typedef void_ptr fftwf_plan;
typedef void_ptr  fftw_plan;


extern "C"
  {
  // single precision (float)
  
  fftwf_plan fftwf_plan_dft_1d(int N, fftwf_complex* input, fftwf_complex* output, int fftw3_sign, unsigned int fftw3_flags);
  
  void      fftwf_execute(fftwf_plan plan);
  void fftwf_destroy_plan(fftwf_plan plan);
  
  void fftwf_cleanup();
  
  
  // double precision (double)
  
  fftw_plan fftw_plan_dft_1d(int N, fftw_complex* input, fftw_complex* output, int fftw3_sign, unsigned int fftw3_flags);
  
  void      fftw_execute(fftw_plan plan);
  void fftw_destroy_plan(fftw_plan plan);
  
  void fftw_cleanup();
  }


#endif
