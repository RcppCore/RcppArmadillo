// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup traits
//! @{


template<typename T1>
struct get_pod_type
  { typedef T1 result; };

template<typename T2>
struct get_pod_type< std::complex<T2> >
  { typedef T2 result; };



template<typename T>
struct is_Mat_only
  { static const bool value = false; };

template<typename eT>
struct is_Mat_only< Mat<eT> >
  { static const bool value = true; };



template<typename T>
struct is_Mat
  { static const bool value = false; };

template<typename eT>
struct is_Mat< Mat<eT> >
  { static const bool value = true; };

template<typename eT>
struct is_Mat< Row<eT> >
  { static const bool value = true; };

template<typename eT>
struct is_Mat< Col<eT> >
  { static const bool value = true; };



template<typename T>
struct is_Row
  { static const bool value = false; };

template<typename eT>
struct is_Row< Row<eT> >
  { static const bool value = true; };



template<typename T>
struct is_Col
  { static const bool value = false; };

template<typename eT>
struct is_Col< Col<eT> >
  { static const bool value = true; };






template<typename T>
struct is_subview
  { static const bool value = false; };

template<typename eT>
struct is_subview< subview<eT> >
  { static const bool value = true; };


template<typename T>
struct is_diagview
  { static const bool value = false; };

template<typename eT>
struct is_diagview< diagview<eT> >
  { static const bool value = true; };


//
//
//



template<typename T>
struct is_Cube
  { static const bool value = false; };

template<typename eT>
struct is_Cube< Cube<eT> >
  { static const bool value = true; };

template<typename T>
struct is_subview_cube
  { static const bool value = false; };

template<typename eT>
struct is_subview_cube< subview_cube<eT> >
  { static const bool value = true; };



//
//
//


template<typename T>
struct is_Gen
  { static const bool value = false; };
 
template<typename eT, typename gen_type>
struct is_Gen< Gen<eT,gen_type> >
  { static const bool value = true; };
 

template<typename T>
struct is_Op
  { static const bool value = false; };
 
template<typename T1, typename op_type>
struct is_Op< Op<T1,op_type> >
  { static const bool value = true; };
 

template<typename T>
struct is_eOp
  { static const bool value = false; };
 
template<typename T1, typename eop_type>
struct is_eOp< eOp<T1,eop_type> >
  { static const bool value = true; };
 

template<typename T>
struct is_mtOp
  { static const bool value = false; };
 
template<typename eT, typename T1, typename op_type>
struct is_mtOp< mtOp<eT, T1, op_type> >
  { static const bool value = true; };
 

template<typename T>
struct is_Glue
  { static const bool value = false; };
 
template<typename T1, typename T2, typename glue_type>
struct is_Glue< Glue<T1,T2,glue_type> >
  { static const bool value = true; };


template<typename T>
struct is_eGlue
  { static const bool value = false; };
 
template<typename T1, typename T2, typename eglue_type>
struct is_eGlue< eGlue<T1,T2,eglue_type> >
  { static const bool value = true; };


template<typename T>
struct is_mtGlue
  { static const bool value = false; };
 
template<typename eT, typename T1, typename T2, typename glue_type>
struct is_mtGlue< mtGlue<eT, T1, T2, glue_type> >
  { static const bool value = true; };


//
//


template<typename T>
struct is_glue_times
  { static const bool value = false; };

template<typename T1, typename T2>
struct is_glue_times< Glue<T1,T2,glue_times> >
  { static const bool value = true; };


template<typename T>
struct is_glue_times_diag
  { static const bool value = false; };

template<typename T1, typename T2>
struct is_glue_times_diag< Glue<T1,T2,glue_times_diag> >
  { static const bool value = true; };


template<typename T>
struct is_op_diagmat
  { static const bool value = false; };
 
template<typename T1>
struct is_op_diagmat< Op<T1,op_diagmat> >
  { static const bool value = true; };


//
//


template<typename T>
struct is_GenCube
  { static const bool value = false; };
 
template<typename eT, typename gen_type>
struct is_GenCube< GenCube<eT,gen_type> >
  { static const bool value = true; };
 

template<typename T>
struct is_OpCube
  { static const bool value = false; };
 
template<typename T1, typename op_type>
struct is_OpCube< OpCube<T1,op_type> >
  { static const bool value = true; };


template<typename T>
struct is_eOpCube
  { static const bool value = false; };
 
template<typename T1, typename eop_type>
struct is_eOpCube< eOpCube<T1,eop_type> >
  { static const bool value = true; };
 

template<typename T>
struct is_mtOpCube
  { static const bool value = false; };
 
template<typename eT, typename T1, typename op_type>
struct is_mtOpCube< mtOpCube<eT, T1, op_type> >
  { static const bool value = true; };
 

template<typename T>
struct is_GlueCube
  { static const bool value = false; };
 
template<typename T1, typename T2, typename glue_type>
struct is_GlueCube< GlueCube<T1,T2,glue_type> >
  { static const bool value = true; };


template<typename T>
struct is_eGlueCube
  { static const bool value = false; };
 
template<typename T1, typename T2, typename eglue_type>
struct is_eGlueCube< eGlueCube<T1,T2,eglue_type> >
  { static const bool value = true; };


template<typename T>
struct is_mtGlueCube
  { static const bool value = false; };
 
template<typename eT, typename T1, typename T2, typename glue_type>
struct is_mtGlueCube< mtGlueCube<eT, T1, T2, glue_type> >
  { static const bool value = true; };


//
//
//


template<typename T>
struct is_op_rel
  { static const bool value = false; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_lt_pre> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_lt_post> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_gt_pre> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_gt_post> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_lteq_pre> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_lteq_post> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_gteq_pre> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_gteq_post> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_eq> >
  { static const bool value = true; };

template<typename out_eT, typename T1>
struct is_op_rel< mtOp<out_eT, T1, op_rel_noteq> >
  { static const bool value = true; };



//
//
//



template<typename T1>
struct is_arma_type
  {
  static const bool value
  =  is_Mat<T1>::value
  || is_Gen<T1>::value
  || is_Op<T1>::value
  || is_eOp<T1>::value
  || is_mtOp<T1>::value
  || is_Glue<T1>::value
  || is_eGlue<T1>::value
  || is_mtGlue<T1>::value
  || is_subview<T1>::value
  || is_diagview<T1>::value
  ;
  };



template<typename T1>
struct is_arma_cube_type
  {
  static const bool value
  =  is_Cube<T1>::value
  || is_GenCube<T1>::value
  || is_OpCube<T1>::value
  || is_eOpCube<T1>::value
  || is_mtOpCube<T1>::value
  || is_GlueCube<T1>::value
  || is_eGlueCube<T1>::value
  || is_mtGlueCube<T1>::value
  || is_subview_cube<T1>::value
  ;
  };



//
//
//


template<typename T1, typename T2>
struct is_same_type
  { static const bool value = false; };


template<typename T1>
struct is_same_type<T1,T1>
  { static const bool value = true; };



//
//
//


template<typename T1>
struct is_u8
  { static const bool value = false; };

template<>
struct is_u8<u8>
  { static const bool value = true; };



template<typename T1>
struct is_s8
  { static const bool value = false; };

template<>
struct is_s8<s8>
  { static const bool value = true; };



template<typename T1>
struct is_u16
  { static const bool value = false; };

template<>
struct is_u16<u16>
  { static const bool value = true; };



template<typename T1>
struct is_s16
  { static const bool value = false; };

template<>
struct is_s16<s16>
  { static const bool value = true; };



template<typename T1>
struct is_u32
  { static const bool value = false; };

template<>
struct is_u32<u32>
  { static const bool value = true; };



template<typename T1>
struct is_s32
  { static const bool value = false; };

template<>
struct is_s32<s32>
  { static const bool value = true; };



#if defined(ARMA_64BIT_WORD)
  template<typename T1>
  struct is_u64
    { static const bool value = false; };

  template<>
  struct is_u64<u64>
    { static const bool value = true; };
  
  
  template<typename T1>
  struct is_s64
    { static const bool value = false; };

  template<>
  struct is_s64<s64>
    { static const bool value = true; };
#endif



template<typename T1>
struct is_uword
  { static const bool value = false; };

template<>
struct is_uword<uword>
  { static const bool value = true; };



template<typename T1>
struct is_sword
  { static const bool value = false; };

template<>
struct is_sword<sword>
  { static const bool value = true; };



template<typename T1>
struct is_float
  { static const bool value = false; };

template<>
struct is_float<float>
  { static const bool value = true; };



template<typename T1>
struct is_double
  { static const bool value = false; };

template<>
struct is_double<double>
  { static const bool value = true; };



template<typename T1>
struct is_complex
  { static const bool value = false; };

// template<>
template<typename eT>
struct is_complex< std::complex<eT> >
  { static const bool value = true; };



template<typename T1>
struct is_complex_float
  { static const bool value = false; };

template<>
struct is_complex_float< std::complex<float> >
  { static const bool value = true; };



template<typename T1>
struct is_complex_double
  { static const bool value = false; };

template<>
struct is_complex_double< std::complex<double> >
  { static const bool value = true; };




//! check for a weird implementation of the std::complex class
template<typename T1>
struct is_supported_complex
  { static const bool value = false; };

//template<>
template<typename eT>
struct is_supported_complex< std::complex<eT> >
  { static const bool value = ( sizeof(std::complex<eT>) == 2*sizeof(eT) ); };



template<typename T1>
struct is_supported_complex_float
  { static const bool value = false; };

template<>
struct is_supported_complex_float< std::complex<float> >
  { static const bool value = ( sizeof(std::complex<float>) == 2*sizeof(float) ); };



template<typename T1>
struct is_supported_complex_double
  { static const bool value = false; };

template<>
struct is_supported_complex_double< std::complex<double> >
  { static const bool value = ( sizeof(std::complex<double>) == 2*sizeof(double) ); };



template<typename T1>
struct is_supported_elem_type
  {
  static const bool value = \
    is_u8<T1>::value ||
    is_s8<T1>::value ||
    is_u16<T1>::value ||
    is_s16<T1>::value ||
    is_u32<T1>::value ||
    is_s32<T1>::value ||
#if defined(ARMA_64BIT_WORD)
    is_u64<T1>::value ||
    is_s64<T1>::value ||
#endif
    is_float<T1>::value ||
    is_double<T1>::value ||
    is_supported_complex_float<T1>::value ||
    is_supported_complex_double<T1>::value;
  };



template<typename T1>
struct is_supported_blas_type
  {
  static const bool value = \
    is_float<T1>::value ||
    is_double<T1>::value ||
    is_supported_complex_float<T1>::value ||
    is_supported_complex_double<T1>::value;
  };



template<typename T>
struct is_signed
  {
  static const bool value = true;
  };


template<> struct is_signed<u8>  { static const bool value = false; };
template<> struct is_signed<u16> { static const bool value = false; };
template<> struct is_signed<u32> { static const bool value = false; };
#if defined(ARMA_64BIT_WORD)
template<> struct is_signed<u64> { static const bool value = false; };
#endif


template<typename T>
struct is_non_integral
  {
  static const bool value = false;
  };


template<> struct is_non_integral<              float   > { static const bool value = true; };
template<> struct is_non_integral<              double  > { static const bool value = true; };
template<> struct is_non_integral< std::complex<float>  > { static const bool value = true; };
template<> struct is_non_integral< std::complex<double> > { static const bool value = true; };




//

class arma_junk_class;

template<typename T1, typename T2>
struct force_different_type
  {
  typedef T1 T1_result;
  typedef T2 T2_result;
  };
  

template<typename T1>
struct force_different_type<T1,T1>
  {
  typedef T1              T1_result;
  typedef arma_junk_class T2_result;
  };
  
  

//! @}
