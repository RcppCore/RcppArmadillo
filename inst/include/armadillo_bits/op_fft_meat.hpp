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



//! \addtogroup op_fft
//! @{


#if defined(ARMA_USE_FFTW3)

template<typename cx_type, bool inverse>
class fft_engine_wrapper
  {
  public:
  
  static constexpr uword threshold = 512;
  
  fft_engine_kissfft<cx_type,inverse>* worker_kissfft = nullptr;
  fft_engine_fftw3  <cx_type,inverse>* worker_fftw3   = nullptr;
  
  inline
  ~fft_engine_wrapper()
    {
    arma_extra_debug_sigprint();
    
    if(worker_kissfft != nullptr)  { delete worker_kissfft; }
    if(worker_fftw3   != nullptr)  { delete worker_fftw3;   }
    }
  
  inline
  fft_engine_wrapper(const uword N_samples, const uword N_exec)
    {
    arma_extra_debug_sigprint();
    
    const bool use_fftw3 = N_samples >= (threshold / N_exec);
    
    worker_kissfft = (use_fftw3 == false) ? new fft_engine_kissfft<cx_type,inverse>(N_samples) : nullptr;
    worker_fftw3   = (use_fftw3 == true ) ? new fft_engine_fftw3  <cx_type,inverse>(N_samples) : nullptr;
    }
  
  inline
  void
  run(cx_type* Y, const cx_type* X)
    {
    arma_extra_debug_sigprint();
    
         if(worker_kissfft != nullptr)  { (*worker_kissfft).run(Y,X); }
    else if(worker_fftw3   != nullptr)  {   (*worker_fftw3).run(Y,X); }
    }
  };

#endif


//
// op_fft_real


template<typename T1>
inline
void
op_fft_real::apply( Mat< std::complex<typename T1::pod_type> >& out, const mtOp<std::complex<typename T1::pod_type>,T1,op_fft_real>& in )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type         in_eT;
  typedef typename std::complex<in_eT> out_eT;
  
  // no need to worry about aliasing, as we're going from a real object to complex complex, which by definition cannot alias
  
  const quasi_unwrap<T1> U(in.m);
  const Mat<in_eT>& X  = U.M;
  
  const uword n_rows = X.n_rows;
  const uword n_cols = X.n_cols;
  const uword n_elem = X.n_elem;
  
  const bool is_vec = ( (n_rows == 1) || (n_cols == 1) );
  
  const uword N_orig = (is_vec)              ? n_elem         : n_rows;
  const uword N_user = (in.aux_uword_b == 0) ? in.aux_uword_a : N_orig;
  
  #if defined(ARMA_USE_FFTW3)
    const uword N_exec = (is_vec) ? uword(1) : n_cols;
    fft_engine_wrapper<out_eT,false> worker(N_user, N_exec);
  #else
    fft_engine_kissfft<out_eT,false> worker(N_user);
  #endif
  
  if(is_vec)
    {
    (n_cols == 1) ? out.set_size(N_user, 1) : out.set_size(1, N_user);
    
    if( (out.n_elem == 0) || (N_orig == 0) )  { out.zeros(); return; }
    
    if( (N_user == 1) && (N_orig >= 1) )  { out[0] = out_eT( X[0] ); return; }
    
    podarray<out_eT> data(N_user, arma_zeros_indicator());
    
          out_eT* data_mem = data.memptr();
    const  in_eT*    X_mem =    X.memptr();
    
    const uword N = (std::min)(N_user, N_orig);
    
    for(uword i=0; i < N; ++i)  { data_mem[i].real(X_mem[i]); }
    
    worker.run( out.memptr(), data_mem );
    }
  else
    {
    // process each column seperately
    
    out.set_size(N_user, n_cols);
    
    if( (out.n_elem == 0) || (N_orig == 0) )  { out.zeros(); return; }
    
    if( (N_user == 1) && (N_orig >= 1) )
      {
      for(uword col=0; col < n_cols; ++col)  { out.at(0,col).real( X.at(0,col) ); }
      
      return;
      }
    
    podarray<out_eT> data(N_user, arma_zeros_indicator());
    
    out_eT* data_mem = data.memptr();
    
    const uword N = (std::min)(N_user, N_orig);
    
    for(uword col=0; col < n_cols; ++col)
      {
      for(uword i=0; i < N; ++i)  { data_mem[i].real( X.at(i, col) ); }
      
      worker.run( out.colptr(col), data_mem );
      }
    }
  }



//
// op_fft_cx


template<typename T1>
inline
void
op_fft_cx::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_fft_cx>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_fft_cx::apply_noalias<eT,false>(tmp, U.M, in.aux_uword_a, in.aux_uword_b);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_fft_cx::apply_noalias<eT,false>(out, U.M, in.aux_uword_a, in.aux_uword_b);
    }
  }
  


template<typename eT, bool inverse>
inline
void
op_fft_cx::apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword a, const uword b)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = X.n_rows;
  const uword n_cols = X.n_cols;
  const uword n_elem = X.n_elem;
  
  const bool is_vec = ( (n_rows == 1) || (n_cols == 1) );
  
  const uword N_orig = (is_vec) ? n_elem : n_rows;
  const uword N_user = (b == 0) ? a      : N_orig;
  
  #if defined(ARMA_USE_FFTW3)
    const uword N_exec = (is_vec) ? uword(1) : n_cols;
    fft_engine_wrapper<eT,inverse> worker(N_user, N_exec);
  #else
    fft_engine_kissfft<eT,inverse> worker(N_user);
  #endif
  
  if(is_vec)
    {
    (n_cols == 1) ? out.set_size(N_user, 1) : out.set_size(1, N_user);
    
    if( (out.n_elem == 0) || (N_orig == 0) )  { out.zeros(); return; }
    
    if( (N_user == 1) && (N_orig >= 1) )  { out[0] = X[0]; return; }
    
    if(N_user > N_orig)
      {
      podarray<eT> data(N_user);
      
      eT* data_mem = data.memptr();
      
      arrayops::fill_zeros( &data_mem[N_orig], (N_user - N_orig) );
      
      arrayops::copy(data_mem, X.memptr(), (std::min)(N_user, N_orig));
      
      worker.run( out.memptr(), data_mem );
      }
    else
      {
      worker.run( out.memptr(), X.memptr() );
      }
    }
  else
    {
    // process each column seperately
    
    out.set_size(N_user, n_cols);
    
    if( (out.n_elem == 0) || (N_orig == 0) )  { out.zeros(); return; }
    
    if( (N_user == 1) && (N_orig >= 1) )
      {
      for(uword col=0; col < n_cols; ++col)  { out.at(0,col) = X.at(0,col); }
      
      return;
      }
    
    if(N_user > N_orig)
      {
      podarray<eT> data(N_user);
      
      eT* data_mem = data.memptr();
      
      arrayops::fill_zeros( &data_mem[N_orig], (N_user - N_orig) );
      
      const uword N = (std::min)(N_user, N_orig);
      
      for(uword col=0; col < n_cols; ++col)
        {
        arrayops::copy(data_mem, X.colptr(col), N);
        
        worker.run( out.colptr(col), data_mem );
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
        {
        worker.run( out.colptr(col), X.colptr(col) );
        }
      }
    }
    
  
  // correct the scaling for the inverse transform
  if(inverse)
    {
    typedef typename get_pod_type<eT>::result T;
    
    const T k = T(1) / T(N_user);
    
    eT* out_mem = out.memptr();
    
    const uword out_n_elem = out.n_elem;
    
    for(uword i=0; i < out_n_elem; ++i)  { out_mem[i] *= k; }
    }
  }



//
// op_ifft_cx


template<typename T1>
inline
void
op_ifft_cx::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_ifft_cx>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_fft_cx::apply_noalias<eT,true>(tmp, U.M, in.aux_uword_a, in.aux_uword_b);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_fft_cx::apply_noalias<eT,true>(out, U.M, in.aux_uword_a, in.aux_uword_b);
    }
  }
  


//! @}
