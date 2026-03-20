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


//! \addtogroup eop_core
//! @{



template<typename eop_type>
struct eop_core
  {
  // matrices
  
  template<typename outT, typename T1> arma_hot inline static void apply(outT& out, const eOp<T1, eop_type>& x);
  
  template<typename T1> arma_hot inline static void apply_inplace_plus (Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_minus(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_schur(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_div  (Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x);
  
  
  // cubes
  
  template<typename T1> arma_hot inline static void apply(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x);
  
  template<typename T1> arma_hot inline static void apply_inplace_plus (Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_minus(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_schur(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_div  (Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_type>& x);
  
  
  // common
  
  template<typename eT> arma_inline static eT process(const eT val, const eT k);
  };


struct eop_use_mp_true  { static constexpr bool use_mp = true;  };
struct eop_use_mp_false { static constexpr bool use_mp = false; };


struct eop_neg               : public eop_core<eop_neg>               , public eop_use_mp_false {};
struct eop_scalar_plus       : public eop_core<eop_scalar_plus>       , public eop_use_mp_false {};
struct eop_scalar_minus_pre  : public eop_core<eop_scalar_minus_pre>  , public eop_use_mp_false {};
struct eop_scalar_minus_post : public eop_core<eop_scalar_minus_post> , public eop_use_mp_false {};
struct eop_scalar_times      : public eop_core<eop_scalar_times>      , public eop_use_mp_false {};
struct eop_scalar_div_pre    : public eop_core<eop_scalar_div_pre>    , public eop_use_mp_false {};
struct eop_scalar_div_post   : public eop_core<eop_scalar_div_post>   , public eop_use_mp_false {};
struct eop_square            : public eop_core<eop_square>            , public eop_use_mp_false {};
struct eop_sqrt              : public eop_core<eop_sqrt>              , public eop_use_mp_true  {};
struct eop_pow               : public eop_core<eop_pow>               , public eop_use_mp_false {};  // for pow(), use_mp is selectively enabled in eop_core_meat.hpp
struct eop_log               : public eop_core<eop_log>               , public eop_use_mp_true  {};
struct eop_log2              : public eop_core<eop_log2>              , public eop_use_mp_true  {};
struct eop_log10             : public eop_core<eop_log10>             , public eop_use_mp_true  {};
struct eop_trunc_log         : public eop_core<eop_trunc_log>         , public eop_use_mp_true  {};
struct eop_log1p             : public eop_core<eop_log1p>             , public eop_use_mp_true  {};
struct eop_exp               : public eop_core<eop_exp>               , public eop_use_mp_true  {};
struct eop_exp2              : public eop_core<eop_exp2>              , public eop_use_mp_true  {};
struct eop_exp10             : public eop_core<eop_exp10>             , public eop_use_mp_true  {};
struct eop_trunc_exp         : public eop_core<eop_trunc_exp>         , public eop_use_mp_true  {};
struct eop_expm1             : public eop_core<eop_expm1>             , public eop_use_mp_true  {};
struct eop_cos               : public eop_core<eop_cos>               , public eop_use_mp_true  {};
struct eop_sin               : public eop_core<eop_sin>               , public eop_use_mp_true  {};
struct eop_tan               : public eop_core<eop_tan>               , public eop_use_mp_true  {};
struct eop_acos              : public eop_core<eop_acos>              , public eop_use_mp_true  {};
struct eop_asin              : public eop_core<eop_asin>              , public eop_use_mp_true  {};
struct eop_atan              : public eop_core<eop_atan>              , public eop_use_mp_true  {};
struct eop_cosh              : public eop_core<eop_cosh>              , public eop_use_mp_true  {};
struct eop_sinh              : public eop_core<eop_sinh>              , public eop_use_mp_true  {};
struct eop_tanh              : public eop_core<eop_tanh>              , public eop_use_mp_true  {};
struct eop_acosh             : public eop_core<eop_acosh>             , public eop_use_mp_true  {};
struct eop_asinh             : public eop_core<eop_asinh>             , public eop_use_mp_true  {};
struct eop_atanh             : public eop_core<eop_atanh>             , public eop_use_mp_true  {};
struct eop_sinc              : public eop_core<eop_sinc>              , public eop_use_mp_true  {};
struct eop_abs               : public eop_core<eop_abs>               , public eop_use_mp_false {};
struct eop_arg               : public eop_core<eop_arg>               , public eop_use_mp_false {};
struct eop_conj              : public eop_core<eop_conj>              , public eop_use_mp_false {};
struct eop_floor             : public eop_core<eop_floor>             , public eop_use_mp_false {};
struct eop_ceil              : public eop_core<eop_ceil>              , public eop_use_mp_false {};
struct eop_round             : public eop_core<eop_round>             , public eop_use_mp_false {};
struct eop_trunc             : public eop_core<eop_trunc>             , public eop_use_mp_false {};
struct eop_sign              : public eop_core<eop_sign>              , public eop_use_mp_false {};
struct eop_cbrt              : public eop_core<eop_cbrt>              , public eop_use_mp_true  {};
struct eop_erf               : public eop_core<eop_erf>               , public eop_use_mp_true  {};
struct eop_erfc              : public eop_core<eop_erfc>              , public eop_use_mp_true  {};
struct eop_lgamma            : public eop_core<eop_lgamma>            , public eop_use_mp_true  {};
struct eop_tgamma            : public eop_core<eop_tgamma>            , public eop_use_mp_true  {};



// the classes below are currently not used; reserved for potential future use
struct eop_log_approx {};
struct eop_exp_approx {};
struct eop_approx_log {};
struct eop_approx_exp {};



//! @}
