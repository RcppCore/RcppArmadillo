// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup gmm_full
//! @{


namespace gmm_priv
{

template<typename eT>
class gmm_full
  {
  public:
  
  arma_aligned const Mat <eT> means;
  arma_aligned const Cube<eT> fcovs;
  arma_aligned const Row <eT> hefts;
  
  //
  //
  
  inline ~gmm_full();
  inline  gmm_full();
  
  inline            gmm_full(const gmm_full& x);
  inline gmm_full& operator=(const gmm_full& x);
  
  inline explicit            gmm_full(const gmm_diag<eT>& x);
  inline          gmm_full& operator=(const gmm_diag<eT>& x);
  
  inline      gmm_full(const uword in_n_dims, const uword in_n_gaus);
  inline void    reset(const uword in_n_dims, const uword in_n_gaus);
  inline void    reset();
  
  template<typename T1, typename T2, typename T3>
  inline void set_params(const Base<eT,T1>& in_means, const BaseCube<eT,T2>& in_fcovs, const Base<eT,T3>& in_hefts);
  
  template<typename T1> inline void set_means(const Base    <eT,T1>& in_means);
  template<typename T1> inline void set_fcovs(const BaseCube<eT,T1>& in_fcovs);
  template<typename T1> inline void set_hefts(const Base    <eT,T1>& in_hefts);
  
  inline uword n_dims() const;
  inline uword n_gaus() const;
  
  inline bool load(const std::string name);
  inline bool save(const std::string name) const;
  
  inline Col<eT> generate()              const;
  inline Mat<eT> generate(const uword N) const;
  
  template<typename T1> inline eT      log_p(const T1& expr, const gmm_empty_arg& junk1 = gmm_empty_arg(), typename enable_if<((is_arma_type<T1>::value) && (resolves_to_colvector<T1>::value == true ))>::result* junk2 = nullptr) const;
  template<typename T1> inline eT      log_p(const T1& expr, const uword gaus_id,                          typename enable_if<((is_arma_type<T1>::value) && (resolves_to_colvector<T1>::value == true ))>::result* junk2 = nullptr) const;
  
  template<typename T1> inline Row<eT> log_p(const T1& expr, const gmm_empty_arg& junk1 = gmm_empty_arg(), typename enable_if<((is_arma_type<T1>::value) && (resolves_to_colvector<T1>::value == false))>::result* junk2 = nullptr) const;
  template<typename T1> inline Row<eT> log_p(const T1& expr, const uword gaus_id,                          typename enable_if<((is_arma_type<T1>::value) && (resolves_to_colvector<T1>::value == false))>::result* junk2 = nullptr) const;
  
  template<typename T1> inline eT  sum_log_p(const Base<eT,T1>& expr)                      const;
  template<typename T1> inline eT  sum_log_p(const Base<eT,T1>& expr, const uword gaus_id) const;
  
  template<typename T1> inline eT  avg_log_p(const Base<eT,T1>& expr)                      const;
  template<typename T1> inline eT  avg_log_p(const Base<eT,T1>& expr, const uword gaus_id) const;
  
  template<typename T1> inline uword   assign(const T1& expr, const gmm_dist_mode& dist, typename enable_if<((is_arma_type<T1>::value) && (resolves_to_colvector<T1>::value == true ))>::result* junk = nullptr) const;
  template<typename T1> inline urowvec assign(const T1& expr, const gmm_dist_mode& dist, typename enable_if<((is_arma_type<T1>::value) && (resolves_to_colvector<T1>::value == false))>::result* junk = nullptr) const;
  
  template<typename T1> inline urowvec  raw_hist(const Base<eT,T1>& expr, const gmm_dist_mode& dist_mode) const;
  template<typename T1> inline Row<eT> norm_hist(const Base<eT,T1>& expr, const gmm_dist_mode& dist_mode) const;
  
  template<typename T1>
  inline
  bool
  learn
    (
    const Base<eT,T1>&    data,
    const uword           n_gaus,
    const gmm_dist_mode&  dist_mode,
    const gmm_seed_mode&  seed_mode,
    const uword           km_iter,
    const uword           em_iter,
    const eT              var_floor,
    const bool            print_mode
    );
  
  
  //
  
  protected:
  
  
  arma_aligned Cube<eT> inv_fcovs;
  arma_aligned Row<eT>  log_det_etc;
  arma_aligned Row<eT>  log_hefts;
  arma_aligned Col<eT>  mah_aux;
  arma_aligned Cube<eT> chol_fcovs;
  
  //
  
  inline void init(const gmm_full&     x);
  inline void init(const gmm_diag<eT>& x);
  
  inline void init(const uword in_n_dim, const uword in_n_gaus);
  
  inline void init_constants(const bool calc_chol = true);
  
  inline umat internal_gen_boundaries(const uword N) const;
  
  inline eT internal_scalar_log_p(const eT* x                     ) const;
  inline eT internal_scalar_log_p(const eT* x, const uword gaus_id) const;
  
  inline Row<eT> internal_vec_log_p(const Mat<eT>& X                     ) const;
  inline Row<eT> internal_vec_log_p(const Mat<eT>& X, const uword gaus_id) const;
  
  inline eT internal_sum_log_p(const Mat<eT>& X                     ) const;
  inline eT internal_sum_log_p(const Mat<eT>& X, const uword gaus_id) const;
  
  inline eT internal_avg_log_p(const Mat<eT>& X                     ) const;
  inline eT internal_avg_log_p(const Mat<eT>& X, const uword gaus_id) const;
  
  inline uword internal_scalar_assign(const Mat<eT>& X, const gmm_dist_mode& dist_mode) const;
  
  inline void internal_vec_assign(urowvec& out, const Mat<eT>& X, const gmm_dist_mode& dist_mode) const;
  
  inline void internal_raw_hist(urowvec& hist, const Mat<eT>& X, const gmm_dist_mode& dist_mode) const;
  
  //
  
  template<uword dist_id> inline void generate_initial_means(const Mat<eT>& X, const gmm_seed_mode& seed);
  
  template<uword dist_id> inline void generate_initial_params(const Mat<eT>& X, const eT var_floor);
  
  template<uword dist_id> inline bool km_iterate(const Mat<eT>& X, const uword max_iter, const bool verbose);
  
  //
  
  inline bool em_iterate(const Mat<eT>& X, const uword max_iter, const eT var_floor, const bool verbose);
  
  inline void em_update_params(const Mat<eT>& X, const umat& boundaries, field< Mat<eT> >& t_acc_means, field< Cube<eT> >& t_acc_fcovs, field< Col<eT> >& t_acc_norm_lhoods, field< Col<eT> >& t_gaus_log_lhoods, Col<eT>& t_progress_log_lhoods, const eT var_floor);
  
  inline void em_generate_acc(const Mat<eT>& X, const uword start_index, const uword end_index, Mat<eT>& acc_means, Cube<eT>& acc_fcovs, Col<eT>& acc_norm_lhoods, Col<eT>& gaus_log_lhoods, eT& progress_log_lhood) const;
  
  inline void em_fix_params(const eT var_floor);
  };

}


typedef gmm_priv::gmm_full<double>  gmm_full;
typedef gmm_priv::gmm_full<float>  fgmm_full;


//! @}
