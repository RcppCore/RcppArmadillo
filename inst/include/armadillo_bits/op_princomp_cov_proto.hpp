// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
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


//! \addtogroup op_princomp_cov
//! @{



class op_princomp_cov
  {
  public:
  
  // real element versions
  
  template<typename eT>
  inline static void
  direct_princomp_cov
    (      
          Mat<eT>& coeff_out,                                                     
    const Mat<eT>& in
    );
                                                           
  template<typename eT>
  inline static void
  direct_princomp_cov
    (      
          Mat<eT>& coeff_out, 
          Col<eT>& latent_out,
    const Mat<eT>& in
    );
                                                           
  template<typename eT>
  inline static void
  direct_princomp_cov
    (
          Mat<eT>& coeff_out, 
          Col<eT>& latent_out,
          Col<eT>& explained_out,                    
    const Mat<eT>& in
    );
  
  // complex element versions
  
  template<typename T>
  inline static void
  direct_princomp_cov
    (
          Mat< std::complex<T> >& coeff_out,
    const Mat< std::complex<T> >& in
    );
    
  template<typename T>
  inline static void
  direct_princomp_cov
    (
          Mat< std::complex<T> >& coeff_out, 
                          Col<T>& latent_out,                                                
    const Mat< std::complex<T> >& in
    );
    
  template<typename T>
  inline static void
  direct_princomp_cov
    (
          Mat< std::complex<T> >& coeff_out, 
                          Col<T>& latent_out,
                          Col<T>& explained_out,                                   
    const Mat< std::complex<T> >& in
    );
  
  template<typename T1>
  inline static void
  apply(Mat<typename T1::elem_type>& out, const Op<T1,op_princomp_cov>& in);
  
  };



//! @}
