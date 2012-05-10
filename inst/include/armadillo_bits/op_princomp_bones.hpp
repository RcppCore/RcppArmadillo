// Copyright (C) 2010-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2012 Conrad Sanderson
// Copyright (C) 2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_princomp
//! @{



class op_princomp
  {
  public:
  
  //
  // real element versions
  
  template<typename T1>
  inline static bool
  direct_princomp
    (
           Mat<typename T1::elem_type>&     coeff_out,
    const Base<typename T1::elem_type, T1>& X,
    const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
    );
  
  template<typename T1>
  inline static bool
  direct_princomp
    (
           Mat<typename T1::elem_type>&     coeff_out,
           Mat<typename T1::elem_type>&     score_out,
    const Base<typename T1::elem_type, T1>& X,
    const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
    );
  
  template<typename T1>
  inline static bool
  direct_princomp
    (
           Mat<typename T1::elem_type>&     coeff_out,
           Mat<typename T1::elem_type>&     score_out,
           Col<typename T1::elem_type>&     latent_out,
    const Base<typename T1::elem_type, T1>& X,
    const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
    );
  
  template<typename T1>
  inline static bool
  direct_princomp
    (
           Mat<typename T1::elem_type>&     coeff_out,
           Mat<typename T1::elem_type>&     score_out,
           Col<typename T1::elem_type>&     latent_out,
           Col<typename T1::elem_type>&     tsquared_out,
    const Base<typename T1::elem_type, T1>& X,
    const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
    );
  
  
  //
  // complex element versions
  
  template<typename T1>
  inline static bool
  direct_princomp
    (
           Mat< std::complex<typename T1::pod_type> >&     coeff_out,
    const Base< std::complex<typename T1::pod_type>, T1 >& X,
    const typename arma_cx_only<typename T1::elem_type>::result* junk = 0
    );
  
  template<typename T1>
  inline static bool
  direct_princomp
    (
           Mat< std::complex<typename T1::pod_type> >&     coeff_out,
           Mat< std::complex<typename T1::pod_type> >&     score_out,
    const Base< std::complex<typename T1::pod_type>, T1 >& X,
    const typename arma_cx_only<typename T1::elem_type>::result* junk = 0
    );
  
  template<typename T1>
  inline static bool
  direct_princomp
    (
           Mat< std::complex<typename T1::pod_type> >&     coeff_out,
           Mat< std::complex<typename T1::pod_type> >&     score_out,
           Col<              typename T1::pod_type  >&     latent_out,
    const Base< std::complex<typename T1::pod_type>, T1 >& X,
    const typename arma_cx_only<typename T1::elem_type>::result* junk = 0
    );
  
  template<typename T1>
  inline static bool
  direct_princomp
    (
           Mat< std::complex<typename T1::pod_type> >&     coeff_out,
           Mat< std::complex<typename T1::pod_type> >&     score_out,
           Col<              typename T1::pod_type  >&     latent_out,
           Col< std::complex<typename T1::pod_type> >&     tsquared_out,
    const Base< std::complex<typename T1::pod_type>, T1 >& X,
    const typename arma_cx_only<typename T1::elem_type>::result* junk = 0
    );
  
  
  template<typename T1>
  inline static void
  apply(Mat<typename T1::elem_type>& out, const Op<T1,op_princomp>& in);
  
  };



//! @}
