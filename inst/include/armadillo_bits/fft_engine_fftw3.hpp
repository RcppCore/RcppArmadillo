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
// 
// ------------------------------------------------------------------------


//! \addtogroup fft_engine_fftw3
//! @{


#if defined(ARMA_USE_FFTW3)

struct fft_engine_fftw3_aux
  {
  #if (!defined(ARMA_DONT_USE_STD_MUTEX))
  static inline std::mutex& get_plan_mutex() { static std::mutex plan_mutex; return plan_mutex; }
  #endif
  };

template<typename cx_type, bool inverse>
class fft_engine_fftw3
  {
  public:
  
  constexpr static int fftw3_sign_forward  = -1;
  constexpr static int fftw3_sign_backward = +1;
  
  constexpr static unsigned int fftw3_flag_destroy  = (1u << 0);
  constexpr static unsigned int fftw3_flag_preserve = (1u << 4);
  constexpr static unsigned int fftw3_flag_estimate = (1u << 6);
  
  const uword N;
  
  void_ptr fftw3_plan;
  
  podarray<cx_type> X_work;  // for storing copy of input (can be overwritten by FFTW3)
  podarray<cx_type> Y_work;  // for storing output
  
  inline
  ~fft_engine_fftw3()
    {
    arma_extra_debug_sigprint();
    
    if(fftw3_plan != nullptr)  { fftw3::destroy_plan<cx_type>(fftw3_plan); }
    
    // fftw3::cleanup<cx_type>();  // NOTE: this also removes any wisdom acquired by FFTW3 
    }
  
  inline
  fft_engine_fftw3(const uword in_N)
    : N         (in_N   )
    , fftw3_plan(nullptr)
    {
    arma_extra_debug_sigprint();
    
    if(N == 0)  { return; }
    
    if(N > uword(std::numeric_limits<int>::max()))
      {
      arma_stop_runtime_error("integer overflow: FFT size too large for integer type used by FFTW3");
      }
    
    arma_extra_debug_print("fft_engine_fftw3::constructor: allocating work arrays");
    X_work.set_size(N);
    Y_work.set_size(N);
    
    const int fftw3_sign  = (inverse) ? fftw3_sign_backward : fftw3_sign_forward;
    const int fftw3_flags = fftw3_flag_destroy | fftw3_flag_estimate;
    
    arma_extra_debug_print("fft_engine_fftw3::constructor: generating 1D plan");
    
    // only fftw3::execute() is thread safe, as per FFTW docs:
    // https://www.fftw.org/fftw3_doc/Thread-safety.html
    
    #if defined(ARMA_USE_OPENMP)
      {
      #pragma omp critical (arma_fft_engine_fftw3)
        {
        fftw3_plan = fftw3::plan_dft_1d<cx_type>(N, X_work.memptr(), Y_work.memptr(), fftw3_sign, fftw3_flags);
        }
      }
    #elif (!defined(ARMA_DONT_USE_STD_MUTEX))
      {
      std::mutex& plan_mutex = fft_engine_fftw3_aux::get_plan_mutex();
      
      const std::lock_guard<std::mutex> lock(plan_mutex);
      
      fftw3_plan = fftw3::plan_dft_1d<cx_type>(N, X_work.memptr(), Y_work.memptr(), fftw3_sign, fftw3_flags);
      }
    #else
      {
      fftw3_plan = fftw3::plan_dft_1d<cx_type>(N, X_work.memptr(), Y_work.memptr(), fftw3_sign, fftw3_flags);
      }
    #endif
    
    if(fftw3_plan == nullptr)  { arma_stop_runtime_error("fft_engine_fftw3::constructor: failed to create plan"); }
    }
  
  inline
  void
  run(cx_type* Y, const cx_type* X)
    {
    arma_extra_debug_sigprint();
    
    if(fftw3_plan == nullptr)  { return; }
    
    arma_extra_debug_print("fft_engine_fftw3::run(): copying input array");
    arrayops::copy(X_work.memptr(), X, N);
    
    arma_extra_debug_print("fft_engine_fftw3::run(): executing plan");
    fftw3::execute<cx_type>(fftw3_plan);
    
    arma_extra_debug_print("fft_engine_fftw3::run(): copying output array");
    arrayops::copy(Y, Y_work.memptr(), N);
    }
  };

#endif


//! @}
