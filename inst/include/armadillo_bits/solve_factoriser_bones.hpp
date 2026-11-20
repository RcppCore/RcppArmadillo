// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup solve_factoriser
//! @{



// internal use only
template<typename eT>
struct solve_factoriser_trimat_helper
  {
  inline static bool train(Mat<eT>& AA, typename get_pod_type<eT>::result& out_rcond, const Mat<eT>& A, const uword layout);
  
  inline static bool solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA, const uword layout);
  };



// internal use only
template<typename eT>
struct solve_factoriser_sympd_helper
  {
  inline static bool train(Mat<eT>& AA, typename get_pod_type<eT>::result& out_rcond, const Mat<              typename get_pod_type<eT>::result  >& A);
  inline static bool train(Mat<eT>& AA, typename get_pod_type<eT>::result& out_rcond, const Mat<std::complex< typename get_pod_type<eT>::result> >& A);
  
  inline static bool solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA);
  };



// internal use only
template<typename eT>
struct solve_factoriser_square_helper
  {
  inline static bool train(Mat<eT>& AA, podarray<blas_int>& ipiv, typename get_pod_type<eT>::result& out_rcond, const Mat<eT>& A);
  
  inline static bool solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA, const podarray<blas_int>& ipiv);
  };



// internal use only
template<typename eT>
struct solve_factoriser_rect_helper
  {
  inline static bool train(Mat<eT>& AA, typename get_pod_type<eT>::result& out_rcond, const Mat<eT>& A);
  
  inline static bool solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA);
  };



// internal use only
template<typename eT>
struct solve_factoriser_approx_helper
  {
  inline static bool train(Mat<eT>& AA, const Mat<eT>& A);
  
  inline static bool solve(Mat<eT>& out, const Mat<eT>& B, const Mat<eT>& AA);
  };



//



class solve_factoriser
  {
  private:
  
  enum class elem_id_enum : uword { elem_none, elem_float, elem_double, elem_cx_float, elem_cx_double };
  enum class fact_id_enum : uword { fact_none, fact_trimat, fact_sympd, fact_square, fact_rect, fact_approx };
  
  template<typename eT> inline void   set_elem_id();
  template<typename eT> inline bool check_elem_id() const;
  
  template<typename eT> inline bool acquire_AA_ptr();
  template<typename eT> inline void release_AA_ptr();
  
  void_ptr AA_ptr_raw       = nullptr;
  double   stored_rcond_val = double(0);
  uword    layout           = 0;  // 0: upper triangular; 1: lower triangular
  uword    n_rows           = 0;  // number of rows in given A matrix
  
  elem_id_enum elem_id = elem_id_enum::elem_none;
  fact_id_enum fact_id = fact_id_enum::fact_none;
  
  podarray<blas_int> ipiv;
  
  
  public:
  
  inline ~solve_factoriser();
  inline  solve_factoriser();
  
  inline void reset();
  
  inline double rcond() const;
  
  template<typename T1> inline bool factorise(const Base<typename T1::elem_type,T1>& A_expr, const solve_opts::opts& opts = solve_opts::none, const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr);
  
  template<typename T1> inline bool solve(Mat<typename T1::elem_type>& X, const Base<typename T1::elem_type,T1>& B_expr, const typename arma_blas_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr) const;
  
  inline      solve_factoriser(const solve_factoriser&) = delete;
  inline void operator=       (const solve_factoriser&) = delete;
  };



//! @}
