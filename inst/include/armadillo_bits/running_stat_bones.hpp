// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup running_stat
//! @{



template<typename eT>
class arma_counter
  {
  public:
  
  inline ~arma_counter();
  inline  arma_counter();
  
  inline const arma_counter& operator++();
  inline void                operator++(int);
  
  inline void reset();
  inline eT   value()         const;
  inline eT   value_plus_1()  const;
  inline eT   value_minus_1() const;
  
  
  private:
  
  arma_aligned eT    d_count;
  arma_aligned uword i_count;
  };



//! Class for keeping statistics of a continuously sampled process / signal.
//! Useful if the storage of individual samples is not necessary or desired.
//! Also useful if the number of samples is not known beforehand or exceeds 
//! available memory.
template<typename eT>
class running_stat
  {
  public:
  
  typedef typename get_pod_type<eT>::result T;
  
  
  inline ~running_stat();
  inline  running_stat();
  
  inline void operator() (const T sample);
  inline void operator() (const std::complex<T>& sample);
  
  inline void reset();
  
  inline eT mean() const;
  
  inline  T var   (const uword norm_type = 0) const;
  inline  T stddev(const uword norm_type = 0) const;
  
  inline eT min()  const;
  inline eT max()  const;
  
  inline T count() const;
  
  //
  //
  
  private:
  
  arma_aligned arma_counter<T> counter;
  
  arma_aligned eT r_mean;
  arma_aligned  T r_var;
  
  arma_aligned eT min_val;
  arma_aligned eT max_val;
  
  arma_aligned  T min_val_norm;
  arma_aligned  T max_val_norm;
  
  
  friend class running_stat_aux;
  };



class running_stat_aux
  {
  public:
  
  template<typename eT>
  inline static void update_stats(running_stat<eT>&               x,  const eT               sample);
  
  template<typename T>
  inline static void update_stats(running_stat< std::complex<T> >& x, const T                sample);
  
  template<typename T>
  inline static void update_stats(running_stat< std::complex<T> >& x, const std::complex<T>& sample);
  
  };



//! @}
