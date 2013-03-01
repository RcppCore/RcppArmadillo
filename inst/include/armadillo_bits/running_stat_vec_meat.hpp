// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup running_stat_vec
//! @{



template<typename eT>
running_stat_vec<eT>::~running_stat_vec()
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
running_stat_vec<eT>::running_stat_vec(const bool in_calc_cov)
  : calc_cov(in_calc_cov)
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
running_stat_vec<eT>::running_stat_vec(const running_stat_vec<eT>& in_rsv)
  : calc_cov    (in_rsv.calc_cov)
  , counter     (in_rsv.counter)
  , r_mean      (in_rsv.r_mean)
  , r_var       (in_rsv.r_var)
  , r_cov       (in_rsv.r_cov)
  , min_val     (in_rsv.min_val)
  , max_val     (in_rsv.max_val)
  , min_val_norm(in_rsv.min_val_norm)
  , max_val_norm(in_rsv.max_val_norm)
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
const running_stat_vec<eT>&
running_stat_vec<eT>::operator=(const running_stat_vec<eT>& in_rsv)
  {
  arma_extra_debug_sigprint();
  
  access::rw(calc_cov) = in_rsv.calc_cov;
  
  counter      = in_rsv.counter;
  r_mean       = in_rsv.r_mean;
  r_var        = in_rsv.r_var;
  r_cov        = in_rsv.r_cov;
  min_val      = in_rsv.min_val;
  max_val      = in_rsv.max_val;
  min_val_norm = in_rsv.min_val_norm;
  max_val_norm = in_rsv.max_val_norm;
  
  return *this;
  }



//! update statistics to reflect new sample
template<typename eT>
template<typename T1>
arma_hot
inline
void
running_stat_vec<eT>::operator() (const Base<typename get_pod_type<eT>::result, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  //typedef typename get_pod_type<eT>::result T;
  
  const unwrap<T1>        tmp(X.get_ref());
  const Mat<eT>& sample = tmp.M;
  
  if( sample.is_empty() )
    {
    return;
    }
  
  if( sample.is_finite() == false )
    {
    arma_warn(true, "running_stat_vec: sample ignored as it has non-finite elements");
    return;
    }
  
  running_stat_vec_aux::update_stats(*this, sample);
  }



//! update statistics to reflect new sample (version for complex numbers)
template<typename eT>
template<typename T1>
arma_hot
inline
void
running_stat_vec<eT>::operator() (const Base<std::complex<typename get_pod_type<eT>::result>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  //typedef typename std::complex<typename get_pod_type<eT>::result> eT;
  
  const unwrap<T1>        tmp(X.get_ref());
  const Mat<eT>& sample = tmp.M;
  
  if( sample.is_empty() )
    {
    return;
    }
  
  if( sample.is_finite() == false )
    {
    arma_warn(true, "running_stat_vec: sample ignored as it has non-finite elements");
    return;
    }
  
  running_stat_vec_aux::update_stats(*this, sample);
  }



//! set all statistics to zero
template<typename eT>
inline
void
running_stat_vec<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  counter.reset();
  
  r_mean.reset();
  r_var.reset();
  r_cov.reset();
  
  min_val.reset();
  max_val.reset();
  
  min_val_norm.reset();
  max_val_norm.reset();
  
  r_var_dummy.reset();
  r_cov_dummy.reset();
  
  tmp1.reset();
  tmp2.reset();
  }



//! mean or average value
template<typename eT>
inline
const Mat<eT>&
running_stat_vec<eT>::mean() const
  {
  arma_extra_debug_sigprint();
  
  return r_mean;
  }



//! variance
template<typename eT>
inline
const Mat<typename get_pod_type<eT>::result>&
running_stat_vec<eT>::var(const uword norm_type)
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
      
      r_var_dummy = (N_minus_1/N) * r_var;
      
      return r_var_dummy;
      }
    }
  else
    {
    r_var_dummy.zeros(r_mean.n_rows, r_mean.n_cols);
    
    return r_var_dummy;
    }
  
  }



//! standard deviation
template<typename eT>
inline
Mat<typename get_pod_type<eT>::result>
running_stat_vec<eT>::stddev(const uword norm_type) const
  {
  arma_extra_debug_sigprint();
  
  const T N = counter.value();
  
  if(N > T(1))
    {
    if(norm_type == 0)
      {
      return sqrt(r_var);
      }
    else
      {
      const T N_minus_1 = counter.value_minus_1();
      
      return sqrt( (N_minus_1/N) * r_var );
      }
    }
  else
    {
    return Mat<T>();
    }
  }



//! covariance
template<typename eT>
inline
const Mat<eT>&
running_stat_vec<eT>::cov(const uword norm_type)
  {
  arma_extra_debug_sigprint();
  
  if(calc_cov == true)
    {
    const T N = counter.value();
    
    if(N > T(1))
      {
      if(norm_type == 0)
        {
        return r_cov;
        }
      else
        {
        const T N_minus_1 = counter.value_minus_1();
        
        r_cov_dummy = (N_minus_1/N) * r_cov;
        
        return r_cov_dummy;
        }
      }
    else
      {
      r_cov_dummy.zeros(r_mean.n_rows, r_mean.n_cols);
      
      return r_cov_dummy;
      }
    }
  else
    {
    r_cov_dummy.reset();
    
    return r_cov_dummy;
    }
  
  }



//! vector with minimum values
template<typename eT>
inline
const Mat<eT>&
running_stat_vec<eT>::min() const
  {
  arma_extra_debug_sigprint();
  
  return min_val;
  }



//! vector with maximum values
template<typename eT>
inline
const Mat<eT>&
running_stat_vec<eT>::max() const
  {
  arma_extra_debug_sigprint();
  
  return max_val;
  }



//! number of samples so far
template<typename eT>
inline
typename get_pod_type<eT>::result
running_stat_vec<eT>::count() const
  {
  arma_extra_debug_sigprint();
  
  return counter.value();
  }



//



//! update statistics to reflect new sample
template<typename eT>
inline
void
running_stat_vec_aux::update_stats(running_stat_vec<eT>& x, const Mat<eT>& sample)
  {
  arma_extra_debug_sigprint();
  
  typedef typename running_stat_vec<eT>::T T;
  
  const T N = x.counter.value();
  
  if(N > T(0))
    {
    arma_debug_assert_same_size(x.r_mean, sample, "running_stat_vec(): dimensionality mismatch");
    
    const uword n_elem      = sample.n_elem;
    const eT*   sample_mem  = sample.memptr();
          eT*   r_mean_mem  = x.r_mean.memptr();
           T*   r_var_mem   = x.r_var.memptr();
          eT*   min_val_mem = x.min_val.memptr();
          eT*   max_val_mem = x.max_val.memptr();
    
    const T  N_plus_1   = x.counter.value_plus_1();
    const T  N_minus_1  = x.counter.value_minus_1();
    
    if(x.calc_cov == true)
      {
      Mat<eT>& tmp1 = x.tmp1;
      Mat<eT>& tmp2 = x.tmp2;
      
      tmp1 = sample - x.r_mean;
      
      if(sample.n_cols == 1)
        {
        tmp2 = tmp1*trans(tmp1);
        }
      else
        {
        tmp2 = trans(tmp1)*tmp1;
        }
      
      x.r_cov *= (N_minus_1/N);
      x.r_cov += tmp2 / N_plus_1;
      }
    
    
    for(uword i=0; i<n_elem; ++i)
      {
      const eT val = sample_mem[i];
      
      if(val < min_val_mem[i])
        {
        min_val_mem[i] = val;
        }
      
      if(val > max_val_mem[i])
        {
        max_val_mem[i] = val;
        }
        
      const eT r_mean_val = r_mean_mem[i];
      const eT tmp        = val - r_mean_val;
    
      r_var_mem[i] = N_minus_1/N * r_var_mem[i] + (tmp*tmp)/N_plus_1;
      
      r_mean_mem[i] = r_mean_val + (val - r_mean_val)/N_plus_1;
      }
    }
  else
    {
    arma_debug_check( (sample.is_vec() == false), "running_stat_vec(): given sample is not a vector");
    
    x.r_mean.set_size(sample.n_rows, sample.n_cols);
    
    x.r_var.zeros(sample.n_rows, sample.n_cols);
    
    if(x.calc_cov == true)
      {
      x.r_cov.zeros(sample.n_elem, sample.n_elem);
      }
    
    x.min_val.set_size(sample.n_rows, sample.n_cols);
    x.max_val.set_size(sample.n_rows, sample.n_cols);
    
    
    const uword n_elem      = sample.n_elem;
    const eT*   sample_mem  = sample.memptr();
          eT*   r_mean_mem  = x.r_mean.memptr();
          eT*   min_val_mem = x.min_val.memptr();
          eT*   max_val_mem = x.max_val.memptr();
          
    
    for(uword i=0; i<n_elem; ++i)
      {
      const eT val = sample_mem[i];
      
      r_mean_mem[i]  = val;
      min_val_mem[i] = val;
      max_val_mem[i] = val;
      }
    }
  
  x.counter++;
  }



//! update statistics to reflect new sample (version for complex numbers)
template<typename T>
inline
void
running_stat_vec_aux::update_stats(running_stat_vec< std::complex<T> >& x, const Mat<T>& sample)
  {
  arma_extra_debug_sigprint();
  
  const Mat< std::complex<T> > tmp = conv_to< Mat< std::complex<T> > >::from(sample);
  
  running_stat_vec_aux::update_stats(x, tmp);
  }



//! alter statistics to reflect new sample (version for complex numbers)
template<typename T>
inline
void
running_stat_vec_aux::update_stats(running_stat_vec< std::complex<T> >& x, const Mat< std::complex<T> >& sample)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const T N = x.counter.value();
  
  if(N > T(0))
    {
    arma_debug_assert_same_size(x.r_mean, sample, "running_stat_vec(): dimensionality mismatch");
    
    const uword n_elem           = sample.n_elem;
    const eT*   sample_mem       = sample.memptr();
          eT*   r_mean_mem       = x.r_mean.memptr();
           T*   r_var_mem        = x.r_var.memptr();
          eT*   min_val_mem      = x.min_val.memptr();
          eT*   max_val_mem      = x.max_val.memptr();
           T*   min_val_norm_mem = x.min_val_norm.memptr();
           T*   max_val_norm_mem = x.max_val_norm.memptr();
    
    const T  N_plus_1   = x.counter.value_plus_1();
    const T  N_minus_1  = x.counter.value_minus_1();
    
    if(x.calc_cov == true)
      {
      Mat<eT>& tmp1 = x.tmp1;
      Mat<eT>& tmp2 = x.tmp2;
      
      tmp1 = sample - x.r_mean;
      
      if(sample.n_cols == 1)
        {
        tmp2 = arma::conj(tmp1)*strans(tmp1);
        }
      else
        {
        tmp2 = trans(tmp1)*tmp1;  //tmp2 = strans(conj(tmp1))*tmp1;
        }
      
      x.r_cov *= (N_minus_1/N);
      x.r_cov += tmp2 / N_plus_1;
      }
    
    
    for(uword i=0; i<n_elem; ++i)
      {
      const eT& val      = sample_mem[i];
      const  T  val_norm = std::norm(val);
      
      if(val_norm < min_val_norm_mem[i])
        {
        min_val_norm_mem[i] = val_norm;
        min_val_mem[i]      = val;
        }
      
      if(val_norm > max_val_norm_mem[i])
        {
        max_val_norm_mem[i] = val_norm;
        max_val_mem[i]      = val;
        }
      
      const eT& r_mean_val = r_mean_mem[i];
      
      r_var_mem[i] = N_minus_1/N * r_var_mem[i] + std::norm(val - r_mean_val)/N_plus_1;
      
      r_mean_mem[i] = r_mean_val + (val - r_mean_val)/N_plus_1;
      }
    
    }
  else
    {
    arma_debug_check( (sample.is_vec() == false), "running_stat_vec(): given sample is not a vector");
    
    x.r_mean.set_size(sample.n_rows, sample.n_cols);
    
    x.r_var.zeros(sample.n_rows, sample.n_cols);
    
    if(x.calc_cov == true)
      {
      x.r_cov.zeros(sample.n_elem, sample.n_elem);
      }
    
    x.min_val.set_size(sample.n_rows, sample.n_cols);
    x.max_val.set_size(sample.n_rows, sample.n_cols);
    
    x.min_val_norm.set_size(sample.n_rows, sample.n_cols);
    x.max_val_norm.set_size(sample.n_rows, sample.n_cols);
    
    
    const uword n_elem           = sample.n_elem;
    const eT*   sample_mem       = sample.memptr();
          eT*   r_mean_mem       = x.r_mean.memptr();
          eT*   min_val_mem      = x.min_val.memptr();
          eT*   max_val_mem      = x.max_val.memptr();
           T*   min_val_norm_mem = x.min_val_norm.memptr();
           T*   max_val_norm_mem = x.max_val_norm.memptr();
    
    for(uword i=0; i<n_elem; ++i)
      {
      const eT& val      = sample_mem[i];
      const  T  val_norm = std::norm(val);
      
      r_mean_mem[i]  = val;
      min_val_mem[i] = val;
      max_val_mem[i] = val;
      
      min_val_norm_mem[i] = val_norm;
      max_val_norm_mem[i] = val_norm;
      }
    }
  
  x.counter++;
  }



//! @}
