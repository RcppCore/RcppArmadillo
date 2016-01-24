// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_kmeans
//! @{


#if defined(ARMA_BAD_COMPILER)


template<typename T1>
inline
bool
kmeans(Mat<typename T1::elem_type>& means, const Base<typename T1::elem_type,T1>&, const uword, const gmm_seed_mode&, const uword, const bool)
  {
  arma_extra_debug_sigprint();
  
  arma_stop("kmeans(): unsupported/inadequate compiler");
  
  means.reset();
  
  return false;
  }


#else


template<typename T1>
inline
typename enable_if2<is_real<typename T1::elem_type>::value, bool>::result
kmeans
  (
         Mat<typename T1::elem_type>&    means,
  const Base<typename T1::elem_type,T1>& data,
  const uword                            k,
  const gmm_seed_mode&                   seed_mode,
  const uword                            n_iter,
  const bool                             print_mode
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  gmm_priv::gmm_diag<eT> model;
  
  const bool status = model.kmeans_wrapper(means, data.get_ref(), k, seed_mode, n_iter, print_mode);
  
  if(status == true)
    {
    means = model.means;
    }
  else
    {
    means.reset();
    }
  
  return status;
  }


#endif




//! @}
