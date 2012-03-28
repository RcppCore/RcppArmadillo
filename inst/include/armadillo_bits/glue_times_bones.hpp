// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup glue_times
//! @{



//! \brief
//! Template metaprogram depth_lhs
//! calculates the number of Glue<Tx,Ty, glue_type> instances on the left hand side argument of Glue<Tx,Ty, glue_type>
//! i.e. it recursively expands each Tx, until the type of Tx is not "Glue<..,.., glue_type>"  (i.e the "glue_type" changes)

template<typename glue_type, typename T1>
struct depth_lhs
  {
  static const uword num = 0;
  };

template<typename glue_type, typename T1, typename T2>
struct depth_lhs< glue_type, Glue<T1,T2,glue_type> >
  {
  static const uword num = 1 + depth_lhs<glue_type, T1>::num;
  };



template<uword N>
struct glue_times_redirect
  {
  template<typename T1, typename T2>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X);
  };


template<>
struct glue_times_redirect<2>
  {
  template<typename T1, typename T2>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X, const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0);
  
  template<typename T1, typename T2>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X, const typename arma_not_blas_type<typename T1::elem_type>::result* junk = 0);
  };


template<>
struct glue_times_redirect<3>
  {
  template<typename T1, typename T2, typename T3>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Glue< Glue<T1,T2,glue_times>,T3,glue_times>& X);
  };


template<>
struct glue_times_redirect<4>
  {
  template<typename T1, typename T2, typename T3, typename T4>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Glue< Glue< Glue<T1,T2,glue_times>, T3, glue_times>, T4, glue_times>& X);
  };



//! Class which implements the immediate multiplication of two or more matrices
class glue_times
  {
  public:
  
  
  template<typename T1, typename T2>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X);
  
  
  template<typename T1>
  arma_hot inline static void apply_inplace(Mat<typename T1::elem_type>& out, const T1& X);
  
  template<typename T1, typename T2>
  arma_hot inline static void apply_inplace_plus(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_times>& X, const sword sign);
  
  template<typename eT1, typename eT2>
  inline static void apply_mixed(Mat<typename promote_type<eT1,eT2>::result>& out, const Mat<eT1>& X, const Mat<eT2>& Y);
  
  
  template<typename eT, const bool do_trans_A, const bool do_trans_B>
  arma_inline static uword  mul_storage_cost(const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename eT, const bool do_trans_A, const bool do_trans_B, const bool do_scalar_times>
  arma_hot inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const eT val);
  
  template<typename eT, const bool do_trans_A, const bool do_trans_B, const bool do_trans_C, const bool do_scalar_times>
  arma_hot inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C, const eT val);
  
  template<typename eT, const bool do_trans_A, const bool do_trans_B, const bool do_trans_C, const bool do_trans_D, const bool do_scalar_times>
  arma_hot inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C, const Mat<eT>& D, const eT val);
  
  };



class glue_times_diag
  {
  public:
  
  template<typename T1, typename T2>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_times_diag>& X);
  
  };



//! @}

