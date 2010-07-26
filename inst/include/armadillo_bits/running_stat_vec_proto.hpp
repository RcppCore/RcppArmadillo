// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup running_stat_vec
//! @{



//! Class for keeping statistics of a continuously sampled process / signal.
//! Useful if the storage of individual samples is not necessary or desired.
//! Also useful if the number of samples is not known beforehand or exceeds 
//! available memory.
template<typename eT>
class running_stat_vec
  {
  public:
  
  typedef typename get_pod_type<eT>::result T;
  
  inline ~running_stat_vec();
  inline  running_stat_vec(const bool in_calc_cov = false);
  
  inline running_stat_vec(const running_stat_vec& in_rsv);
  
  inline const running_stat_vec& operator=(const running_stat_vec& in_rsv);
  
  template<typename T1> arma_hot inline void operator() (const Base<               T,  T1>& X);
  template<typename T1> arma_hot inline void operator() (const Base< std::complex<T>,  T1>& X);
  
  inline void reset();
  
  inline const Mat<eT>&  mean() const;
  
  inline const Mat< T>&  var   (const u32 norm_type = 0);
  inline       Mat< T>   stddev(const u32 norm_type = 0) const;
  inline const Mat<eT>&  cov   (const u32 norm_type = 0);
  
  inline const Mat<eT>& min() const;
  inline const Mat<eT>& max() const;
  
  //
  //
  
  private:
  
  const bool calc_cov;
  
  arma_aligned arma_counter<T> counter;
  
  arma_aligned Mat<eT> r_mean;
  arma_aligned Mat< T> r_var;
  arma_aligned Mat<eT> r_cov;
  
  arma_aligned Mat<eT> min_val;
  arma_aligned Mat<eT> max_val;
  
  arma_aligned Mat< T> min_val_norm;
  arma_aligned Mat< T> max_val_norm;
  
  arma_aligned Mat< T> r_var_dummy;
  arma_aligned Mat<eT> r_cov_dummy;
  
  arma_aligned Mat<eT> tmp1;
  arma_aligned Mat<eT> tmp2;
  
  friend class running_stat_vec_aux;
  };



class running_stat_vec_aux
  {
  public:
  
  template<typename eT>
  inline static void update_stats(running_stat_vec< eT >&              x, const Mat<eT>&                sample);
  
  template<typename T>
  inline static void update_stats(running_stat_vec< std::complex<T> >& x, const Mat< T>&                sample);
  
  template<typename T>
  inline static void update_stats(running_stat_vec< std::complex<T> >& x, const Mat< std::complex<T> >& sample);
  
  //
  
  template<typename eT>
  inline static Mat<eT> var(const running_stat_vec< eT >&              x, const u32 norm_type = 0);
  
  template<typename T>
  inline static Mat< T> var(const running_stat_vec< std::complex<T> >& x, const u32 norm_type = 0);
  
  //
  
  template<typename eT>
  inline static Mat<              eT > cov(const running_stat_vec< eT >&              x, const u32 norm_type = 0);
  
  template<typename T>
  inline static Mat< std::complex<T> > cov(const running_stat_vec< std::complex<T> >& x, const u32 norm_type = 0);
  };



//! @}
