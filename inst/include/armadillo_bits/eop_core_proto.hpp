// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
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


//! \addtogroup eop_core
//! @{



template<typename eop_type>
class eop_core
  {
  public:
  
  // matrices
  
  template<typename T1> arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const eOp<T1, eop_type>& x);
  
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
  
  template<typename eT> arma_hot arma_pure arma_inline static eT process(const eT val, const eT k);
  };



class eop_neg               : public eop_core<eop_neg>               {};
class eop_scalar_plus       : public eop_core<eop_scalar_plus>       {};
class eop_scalar_minus_pre  : public eop_core<eop_scalar_minus_pre>  {};
class eop_scalar_minus_post : public eop_core<eop_scalar_minus_post> {};
class eop_scalar_times      : public eop_core<eop_scalar_times>      {};
class eop_scalar_div_pre    : public eop_core<eop_scalar_div_pre>    {};
class eop_scalar_div_post   : public eop_core<eop_scalar_div_post>   {};
class eop_square            : public eop_core<eop_square>            {};
class eop_sqrt              : public eop_core<eop_sqrt>              {};
class eop_log               : public eop_core<eop_log>               {};
class eop_log2              : public eop_core<eop_log2>              {};
class eop_log10             : public eop_core<eop_log10>             {};
class eop_trunc_log         : public eop_core<eop_trunc_log>         {};
class eop_exp               : public eop_core<eop_exp>               {};
class eop_exp2              : public eop_core<eop_exp2>              {};
class eop_exp10             : public eop_core<eop_exp10>             {};
class eop_trunc_exp         : public eop_core<eop_trunc_exp>         {};
class eop_cos               : public eop_core<eop_cos>               {};
class eop_sin               : public eop_core<eop_sin>               {};
class eop_tan               : public eop_core<eop_tan>               {};
class eop_acos              : public eop_core<eop_acos>              {};
class eop_asin              : public eop_core<eop_asin>              {};
class eop_atan              : public eop_core<eop_atan>              {};
class eop_cosh              : public eop_core<eop_cosh>              {};
class eop_sinh              : public eop_core<eop_sinh>              {};
class eop_tanh              : public eop_core<eop_tanh>              {};
class eop_acosh             : public eop_core<eop_acosh>             {};
class eop_asinh             : public eop_core<eop_asinh>             {};
class eop_atanh             : public eop_core<eop_atanh>             {};
class eop_eps               : public eop_core<eop_eps>               {};
class eop_abs               : public eop_core<eop_abs>               {};
class eop_conj              : public eop_core<eop_conj>              {};
class eop_pow               : public eop_core<eop_pow>               {};
class eop_pow_int           : public eop_core<eop_pow_int>           {};
class eop_randu             : public eop_core<eop_randu>             {};
class eop_randn             : public eop_core<eop_randn>             {};
class eop_zeros             : public eop_core<eop_zeros>             {};
class eop_ones_full         : public eop_core<eop_ones_full>         {};
class eop_ones_diag         : public eop_core<eop_ones_diag>         {};

template<typename T1> struct is_generator                { static const bool value = false; };
template<>            struct is_generator<eop_randu>     { static const bool value = true;  };
template<>            struct is_generator<eop_randn>     { static const bool value = true;  };
template<>            struct is_generator<eop_zeros>     { static const bool value = true;  };
template<>            struct is_generator<eop_ones_full> { static const bool value = true;  };
template<>            struct is_generator<eop_ones_diag> { static const bool value = true;  };



//! @}
