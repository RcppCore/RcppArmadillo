// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
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
  
  // real element versions
  
  template<typename eT>
  inline static bool
  direct_princomp
    (      
          Mat<eT>& coeff_out,                                                     
    const Mat<eT>& in
    );
                                                           
  template<typename eT>
  inline static bool
  direct_princomp
    (      
          Mat<eT>& coeff_out, 
          Mat<eT>& score_out,
    const Mat<eT>& in
    );
                                                           
  template<typename eT>
  inline static bool
  direct_princomp
    (
          Mat<eT>& coeff_out, 
          Mat<eT>& score_out, 
          Col<eT>& latent_out,                    
    const Mat<eT>& in
    );
    
  template<typename eT>
  inline static bool
  direct_princomp
    (
          Mat<eT>& coeff_out, 
          Mat<eT>& score_out, 
          Col<eT>& latent_out, 
          Col<eT>& tsquared_out, 
    const Mat<eT>& in
    );  
  
  
  // complex element versions
  
  template<typename T>
  inline static bool
  direct_princomp
    (
          Mat< std::complex<T> >& coeff_out,
    const Mat< std::complex<T> >& in
    );
    
  template<typename T>
  inline static bool
  direct_princomp
    (
          Mat< std::complex<T> >& coeff_out, 
          Mat< std::complex<T> >& score_out,                                                   
    const Mat< std::complex<T> >& in
    );
    
  template<typename T>
  inline static bool
  direct_princomp
    (
          Mat< std::complex<T> >& coeff_out, 
          Mat< std::complex<T> >& score_out, 
                          Col<T>& latent_out,                                   
    const Mat< std::complex<T> >& in
    );
    
  template<typename T>
  inline static bool
  direct_princomp
    (
          Mat< std::complex<T> >& coeff_out, 
          Mat< std::complex<T> >& score_out, 
                          Col<T>& latent_out, 
          Col< std::complex<T> >& tsquared_out, 
    const Mat< std::complex<T> >& in    
    );
  
  
  template<typename T1>
  inline static void
  apply(Mat<typename T1::elem_type>& out, const Op<T1,op_princomp>& in);
  
  };



//! @}
