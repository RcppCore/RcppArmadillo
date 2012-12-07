// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup SpGlue
//! @{



template<typename T1, typename T2, typename spglue_type>
class SpGlue : public SpBase<typename T1::elem_type, SpGlue<T1, T2, spglue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool is_row = ( (T1::is_row || T2::is_row) && is_spglue_elem<spglue_type>::value ) || ( (is_spglue_times<spglue_type>::value || is_spglue_times2<spglue_type>::value) ? T1::is_row : false );
  static const bool is_col = ( (T1::is_col || T2::is_col) && is_spglue_elem<spglue_type>::value ) || ( (is_spglue_times<spglue_type>::value || is_spglue_times2<spglue_type>::value) ? T2::is_col : false );
  
  arma_inline  SpGlue(const T1& in_A, const T2& in_B);
  arma_inline  SpGlue(const T1& in_A, const T2& in_B, const elem_type in_aux);
  arma_inline ~SpGlue();
  
  const T1&       A;    //!< first operand
  const T2&       B;    //!< second operand
        elem_type aux;
  };



//! @}
