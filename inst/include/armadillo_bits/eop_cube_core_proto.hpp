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


//! \addtogroup eop_cube_core
//! @{



template<typename eop_cube_type>
class eop_cube_core
  {
  public:
  
  template<typename T1> arma_hot arma_inline static typename T1::elem_type get_elem(const eOpCube<T1, eop_cube_type>& x, const u32 i);
  template<typename T1> arma_hot arma_inline static typename T1::elem_type get_elem(const eOpCube<T1, eop_cube_type>& x, const u32 row, const u32 col, const u32 slice);
  
  template<typename T1> arma_hot arma_inline static typename T1::elem_type process(const eOpCube<T1, eop_cube_type>& x, const typename T1::elem_type val);
  
  template<typename T1> arma_hot arma_inline static void apply(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x);
  
  template<typename T1> arma_hot inline static void apply_proxy (Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x);
  template<typename T1> arma_hot inline static void apply_unwrap(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x);
  
  template<typename T1> arma_hot inline static void apply_inplace_plus (Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_minus(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_schur(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x);
  template<typename T1> arma_hot inline static void apply_inplace_div  (Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x);

  };



class eop_cube_neg               : public eop_cube_core<eop_cube_neg>               {};
class eop_cube_scalar_plus       : public eop_cube_core<eop_cube_scalar_plus>       {};
class eop_cube_scalar_minus_pre  : public eop_cube_core<eop_cube_scalar_minus_pre>  {};
class eop_cube_scalar_minus_post : public eop_cube_core<eop_cube_scalar_minus_post> {};
class eop_cube_scalar_times      : public eop_cube_core<eop_cube_scalar_times>      {};
class eop_cube_scalar_div_pre    : public eop_cube_core<eop_cube_scalar_div_pre>    {};
class eop_cube_scalar_div_post   : public eop_cube_core<eop_cube_scalar_div_post>   {};
class eop_cube_square            : public eop_cube_core<eop_cube_square>            {};
class eop_cube_sqrt              : public eop_cube_core<eop_cube_sqrt>              {};
class eop_cube_log10             : public eop_cube_core<eop_cube_log10>             {};
class eop_cube_log               : public eop_cube_core<eop_cube_log>               {};
class eop_cube_trunc_log         : public eop_cube_core<eop_cube_trunc_log>         {};
class eop_cube_exp               : public eop_cube_core<eop_cube_exp>               {};
class eop_cube_trunc_exp         : public eop_cube_core<eop_cube_trunc_exp>         {};
class eop_cube_cos               : public eop_cube_core<eop_cube_cos>               {};
class eop_cube_cosh              : public eop_cube_core<eop_cube_cosh>              {};
class eop_cube_acos              : public eop_cube_core<eop_cube_acos>              {};
class eop_cube_acosh             : public eop_cube_core<eop_cube_acosh>             {};
class eop_cube_sin               : public eop_cube_core<eop_cube_sin>               {};
class eop_cube_sinh              : public eop_cube_core<eop_cube_sinh>              {};
class eop_cube_asin              : public eop_cube_core<eop_cube_asin>              {};
class eop_cube_asinh             : public eop_cube_core<eop_cube_asinh>             {};
class eop_cube_tan               : public eop_cube_core<eop_cube_tan>               {};
class eop_cube_tanh              : public eop_cube_core<eop_cube_tanh>              {};
class eop_cube_atan              : public eop_cube_core<eop_cube_atan>              {};
class eop_cube_atanh             : public eop_cube_core<eop_cube_atanh>             {};
class eop_cube_eps               : public eop_cube_core<eop_cube_eps>               {};
class eop_cube_abs               : public eop_cube_core<eop_cube_abs>               {};
class eop_cube_conj              : public eop_cube_core<eop_cube_conj>              {};
class eop_cube_pow               : public eop_cube_core<eop_cube_pow>               {};
class eop_cube_pow_int           : public eop_cube_core<eop_cube_pow_int>           {};



class eop_cube_ones_full : public eop_cube_core<eop_cube_ones_full> {};
class eop_cube_rand      : public eop_cube_core<eop_cube_rand>      {};
class eop_cube_randn     : public eop_cube_core<eop_cube_randn>     {};
class eop_cube_zeros     : public eop_cube_core<eop_cube_zeros>     {};

template<typename T1> struct is_cube_generator                     { static const bool value = false; };
template<>            struct is_cube_generator<eop_cube_ones_full> { static const bool value = true;  };
template<>            struct is_cube_generator<eop_cube_rand>      { static const bool value = true;  };
template<>            struct is_cube_generator<eop_cube_randn>     { static const bool value = true;  };
template<>            struct is_cube_generator<eop_cube_zeros>     { static const bool value = true;  };



//! @}
