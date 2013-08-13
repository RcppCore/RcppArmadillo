// Copyright (C) 2009-2011 Conrad Sanderson
// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup running_stat
//! @{



template<typename eT>
inline
arma_counter<eT>::~arma_counter()
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
inline
arma_counter<eT>::arma_counter()
  : d_count(   eT(0))
  , i_count(uword(0))
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
inline
const arma_counter<eT>& 
arma_counter<eT>::operator++()
  {
  if(i_count < ARMA_MAX_UWORD)
    {
    i_count++;
    }
  else
    {
    d_count += eT(ARMA_MAX_UWORD);
    i_count  = 1;
    }
  
  return *this;
  }



template<typename eT>
inline
void
arma_counter<eT>::operator++(int)
  {
  operator++();
  }



template<typename eT>
inline
void
arma_counter<eT>::reset()
  {
  d_count =    eT(0);
  i_count = uword(0);
  }



template<typename eT>
inline
eT
arma_counter<eT>::value() const
  {
  return d_count + eT(i_count);
  }



template<typename eT>
inline
eT
arma_counter<eT>::value_plus_1() const
  {
  if(i_count < ARMA_MAX_UWORD)
    {
    return d_count + eT(i_count + 1);
    }
  else
    {
    return d_count + eT(ARMA_MAX_UWORD) + eT(1);
    }
  }



template<typename eT>
inline
eT
arma_counter<eT>::value_minus_1() const
  {
  if(i_count > 0)
    {
    return d_count + eT(i_count - 1);
    }
  else
    {
    return d_count - eT(1);
    }
  }



//



template<typename eT>
running_stat<eT>::~running_stat()
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
running_stat<eT>::running_stat()
  : r_mean      (                          eT(0))
  , r_var       (typename running_stat<eT>::T(0))
  , min_val     (                          eT(0))
  , max_val     (                          eT(0))
  , min_val_norm(typename running_stat<eT>::T(0))
  , max_val_norm(typename running_stat<eT>::T(0))
  {
  arma_extra_debug_sigprint_this(this);
  }



//! update statistics to reflect new sample
template<typename eT>
inline
void
running_stat<eT>::operator() (const typename running_stat<eT>::T sample)
  {
  arma_extra_debug_sigprint();
  
  if( arma_isfinite(sample) == false )
    {
    arma_warn(true, "running_stat: sample ignored as it is non-finite" );
    return;
    }
  
  running_stat_aux::update_stats(*this, sample);
  }



//! update statistics to reflect new sample (version for complex numbers)
template<typename eT>
inline
void
running_stat<eT>::operator() (const std::complex< typename running_stat<eT>::T >& sample)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type<eT, std::complex< typename running_stat<eT>::T > >::no ));
  
  if( arma_isfinite(sample) == false )
    {
    arma_warn(true, "running_stat: sample ignored as it is non-finite" );
    return;
    }
  
  running_stat_aux::update_stats(*this, sample);
  }



//! set all statistics to zero
template<typename eT>
inline
void
running_stat<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  // typedef typename running_stat<eT>::T T;
  
  counter.reset();
  
  r_mean       = eT(0);
  r_var        =  T(0);
  
  min_val      = eT(0);
  max_val      = eT(0);
  
  min_val_norm =  T(0);
  max_val_norm =  T(0);
  }



//! mean or average value
template<typename eT>
inline
eT
running_stat<eT>::mean() const
  {
  arma_extra_debug_sigprint();
  
  return r_mean;
  }



//! variance
template<typename eT>
inline
typename running_stat<eT>::T
running_stat<eT>::var(const uword norm_type) const
  {
  arma_extra_debug_sigprint();
  
  const T N = counter.value();
  
  if(N > T(1))
    {
    if(norm_type == 0)
      {
      return r_var;
      }
    else
      {
      const T N_minus_1 = counter.value_minus_1();
      return (N_minus_1/N) * r_var;
      }
    }
  else
    {
    return T(0);
    }
  }



//! standard deviation
template<typename eT>
inline
typename running_stat<eT>::T
running_stat<eT>::stddev(const uword norm_type) const
  {
  arma_extra_debug_sigprint();
  
  return std::sqrt( (*this).var(norm_type) );
  }



//! minimum value
template<typename eT>
inline
eT
running_stat<eT>::min() const
  {
  arma_extra_debug_sigprint();
  
  return min_val;
  }



//! maximum value
template<typename eT>
inline
eT
running_stat<eT>::max() const
  {
  arma_extra_debug_sigprint();
  
  return max_val;
  }



//! number of samples so far
template<typename eT>
inline
typename get_pod_type<eT>::result
running_stat<eT>::count() const
  {
  arma_extra_debug_sigprint();
  
  return counter.value();
  }



//! update statistics to reflect new sample
template<typename eT>
inline
void
running_stat_aux::update_stats(running_stat<eT>& x, const eT sample)
  {
  arma_extra_debug_sigprint();
  
  typedef typename running_stat<eT>::T T;
  
  const T N = x.counter.value();
  
  if(N > T(0))
    {
    if(sample < x.min_val)
      {
      x.min_val = sample;
      }
    
    if(sample > x.max_val)
      {
      x.max_val = sample;
      }
    
    const T  N_plus_1   = x.counter.value_plus_1();
    const T  N_minus_1  = x.counter.value_minus_1();
    
    // note: variance has to be updated before the mean
    
    const eT tmp = sample - x.r_mean;
    
    x.r_var  = N_minus_1/N * x.r_var + (tmp*tmp)/N_plus_1;
    
    x.r_mean = x.r_mean + (sample - x.r_mean)/N_plus_1;
    //x.r_mean = (N/N_plus_1)*x.r_mean + sample/N_plus_1;
    //x.r_mean = (x.r_mean + sample/N) * N/N_plus_1;
    }
  else
    {
    x.r_mean  = sample;
    x.min_val = sample;
    x.max_val = sample;
    
    // r_var is initialised to zero
    // in the constructor and reset()
    }
  
  x.counter++;
  }



//! update statistics to reflect new sample (version for complex numbers)
template<typename T>
inline
void
running_stat_aux::update_stats(running_stat< std::complex<T> >& x, const T sample)
  {
  arma_extra_debug_sigprint();

  running_stat_aux::update_stats(x, std::complex<T>(sample));
  }



//! alter statistics to reflect new sample (version for complex numbers)
template<typename T>
inline
void
running_stat_aux::update_stats(running_stat< std::complex<T> >& x, const std::complex<T>& sample)
  {
  arma_extra_debug_sigprint();
  
  const T sample_norm = std::norm(sample);
  const T N           = x.counter.value();
  
  if(N > T(0))
    {
    if(sample_norm < x.min_val_norm)
      {
      x.min_val_norm = sample_norm;
      x.min_val      = sample;
      }
    
    if(sample_norm > x.max_val_norm)
      {
      x.max_val_norm = sample_norm;
      x.max_val      = sample;
      }
    
    const T  N_plus_1   = x.counter.value_plus_1();
    const T  N_minus_1  = x.counter.value_minus_1();
    
    x.r_var = N_minus_1/N * x.r_var + std::norm(sample - x.r_mean)/N_plus_1;
    
    x.r_mean = x.r_mean + (sample - x.r_mean)/N_plus_1;
    //x.r_mean = (N/N_plus_1)*x.r_mean + sample/N_plus_1;
    //x.r_mean = (x.r_mean + sample/N) * N/N_plus_1;
    }
  else
    {
    x.r_mean       = sample;
    x.min_val      = sample;
    x.max_val      = sample;
    x.min_val_norm = sample_norm;
    x.max_val_norm = sample_norm;
    
    // r_var is initialised to zero
    // in the constructor and reset()
    }
  
  x.counter++;
  }



//! @}
