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
template<typename T1>
arma_hot
arma_inline
typename T1::elem_type
eop_cube_core<eop_cube_type>::get_elem(const eOpCube<T1, eop_cube_type>& x, const u32 i)
  {
  typedef typename T1::elem_type eT;

  if(is_cube_generator<eop_cube_type>::value == true) { return eop_aux::generate<eT,eop_cube_type>();            }
  else                                                { return eop_cube_core<eop_cube_type>::process(x, x.P[i]); }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
arma_inline
typename T1::elem_type
eop_cube_core<eop_cube_type>::get_elem(const eOpCube<T1, eop_cube_type>& x, const u32 row, const u32 col, const u32 slice)
  {
  typedef typename T1::elem_type eT;

  if(is_cube_generator<eop_cube_type>::value == true) { return eop_aux::generate<eT,eop_cube_type>();                           }
  else                                                { return eop_cube_core<eop_cube_type>::process(x, x.P.at(row,col,slice)); }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
arma_inline
typename T1::elem_type
eop_cube_core<eop_cube_type>::process(const eOpCube<T1, eop_cube_type>& x, const typename T1::elem_type val)
  {
  typedef typename T1::elem_type eT;

  // the optimiser will keep only one return statement

       if(is_same_type<eop_cube_type, eop_cube_neg              >::value == true) { return -val;                    }
  else if(is_same_type<eop_cube_type, eop_cube_scalar_plus      >::value == true) { return val + x.aux;             }
  else if(is_same_type<eop_cube_type, eop_cube_scalar_minus_pre >::value == true) { return x.aux - val;             }
  else if(is_same_type<eop_cube_type, eop_cube_scalar_minus_post>::value == true) { return val - x.aux;             }
  else if(is_same_type<eop_cube_type, eop_cube_scalar_times     >::value == true) { return val * x.aux;             }
  else if(is_same_type<eop_cube_type, eop_cube_scalar_div_pre   >::value == true) { return x.aux / val;             }
  else if(is_same_type<eop_cube_type, eop_cube_scalar_div_post  >::value == true) { return val / x.aux;             }
  else if(is_same_type<eop_cube_type, eop_cube_square           >::value == true) { return val*val;                 }
  else if(is_same_type<eop_cube_type, eop_cube_sqrt             >::value == true) { return eop_aux::sqrt(val);      }
  else if(is_same_type<eop_cube_type, eop_cube_log10            >::value == true) { return eop_aux::log10(val);     }
  else if(is_same_type<eop_cube_type, eop_cube_log              >::value == true) { return eop_aux::log(val);       }
  else if(is_same_type<eop_cube_type, eop_cube_trunc_log        >::value == true) { return    arma::trunc_log(val); }
  else if(is_same_type<eop_cube_type, eop_cube_exp              >::value == true) { return eop_aux::exp(val);       }
  else if(is_same_type<eop_cube_type, eop_cube_trunc_exp        >::value == true) { return    arma::trunc_exp(val); }
  else if(is_same_type<eop_cube_type, eop_cube_cos              >::value == true) { return eop_aux::cos(val);       }
  else if(is_same_type<eop_cube_type, eop_cube_cosh             >::value == true) { return eop_aux::cosh(val);      }
  else if(is_same_type<eop_cube_type, eop_cube_acos             >::value == true) { return eop_aux::acos(val);      }
#if !defined(ARMA_OLD_MINGW)
  else if(is_same_type<eop_cube_type, eop_cube_acosh            >::value == true) { return eop_aux::acosh(val);     }
#endif
  else if(is_same_type<eop_cube_type, eop_cube_sin              >::value == true) { return eop_aux::sin(val);       }
  else if(is_same_type<eop_cube_type, eop_cube_sinh             >::value == true) { return eop_aux::sinh(val);      }
  else if(is_same_type<eop_cube_type, eop_cube_asin             >::value == true) { return eop_aux::asin(val);      }
#if !defined(ARMA_OLD_MINGW)
  else if(is_same_type<eop_cube_type, eop_cube_asinh            >::value == true) { return eop_aux::asinh(val);     }
#endif
  else if(is_same_type<eop_cube_type, eop_cube_tan              >::value == true) { return eop_aux::tan(val);       }
  else if(is_same_type<eop_cube_type, eop_cube_tanh             >::value == true) { return eop_aux::tanh(val);      }
  else if(is_same_type<eop_cube_type, eop_cube_atan             >::value == true) { return eop_aux::atan(val);      }
#if !defined(ARMA_OLD_MINGW)
  else if(is_same_type<eop_cube_type, eop_cube_atanh            >::value == true) { return eop_aux::atanh(val);     }
#endif
  else if(is_same_type<eop_cube_type, eop_cube_eps              >::value == true) { return eop_aux::direct_eps(val);}
  else if(is_same_type<eop_cube_type, eop_cube_abs              >::value == true) { return eop_aux::arma_abs(val);  }
  else if(is_same_type<eop_cube_type, eop_cube_conj             >::value == true) { return eop_aux::conj(val);      }
  else if(is_same_type<eop_cube_type, eop_cube_pow              >::value == true) { return eop_aux::pow(val, x.aux);}
  else if(is_same_type<eop_cube_type, eop_cube_pow_int          >::value == true)
    {
    const int exponent = (x.aux_u32_b == 0) ? int(x.aux_u32_a) : -int(x.aux_u32_a);

    return eop_aux::pow_int(val, exponent);
    }
  else
    {
    arma_stop("eop_cube_core::process(): unhandled eop_cube_type");
    return eT(0);
    }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
inline
void
eop_cube_core<eop_cube_type>::apply(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x)
  {
  arma_extra_debug_sigprint();

  if(is_Cube<T1>::value == true)
    {
    eop_cube_core<eop_cube_type>::apply_unwrap(out, x);
    }
  else
    {
    eop_cube_core<eop_cube_type>::apply_proxy(out, x);
    }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
inline
void
eop_cube_core<eop_cube_type>::apply_proxy(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x)
  {
  arma_extra_debug_sigprint();

  // eop_cube_type::apply() function is not allowed to unwrap things
  // (in order to get the input into a common format).
  // the proxy class is already providing objects with element access

  typedef typename T1::elem_type eT;

  const ProxyCube<T1>& P = x.P;

  out.set_size(P.n_rows, P.n_cols, P.n_slices);

        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;

  if(is_cube_generator<eop_cube_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] = eop_aux::generate<eT,eop_cube_type>();
      }
    }
  else
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] = eop_cube_core<eop_cube_type>::process(x, P[i]);
      }
    }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
inline
void
eop_cube_core<eop_cube_type>::apply_unwrap(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const ProxyCube<T1>& P = x.P;

  out.set_size(P.n_rows, P.n_cols, P.n_slices);

        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;

  if(is_cube_generator<eop_cube_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] = eop_aux::generate<eT,eop_cube_type>();
      }
    }
  else
    {
    const unwrap_cube<typename ProxyCube<T1>::stored_type> tmp(P.Q);

    const Cube<eT>& A     = tmp.M;
    const eT*       A_mem = A.memptr();

    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] = eop_cube_core<eop_cube_type>::process(x, A_mem[i]);
      }
    }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
inline
void
eop_cube_core<eop_cube_type>::apply_inplace_plus(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const ProxyCube<T1>& P = x.P;

  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, P.n_rows, P.n_cols, P.n_slices, "cube addition");

        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;

  if(is_cube_generator<eop_cube_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] += eop_aux::generate<eT,eop_cube_type>();
      }
    }
  else
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] += eop_cube_core<eop_cube_type>::process(x, P[i]);
      }
    }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
inline
void
eop_cube_core<eop_cube_type>::apply_inplace_minus(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const ProxyCube<T1>& P = x.P;

  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, P.n_rows, P.n_cols, P.n_slices, "cube subtraction");

        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;

  if(is_cube_generator<eop_cube_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] -= eop_aux::generate<eT,eop_cube_type>();
      }
    }
  else
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] -= eop_cube_core<eop_cube_type>::process(x, P[i]);
      }
    }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
inline
void
eop_cube_core<eop_cube_type>::apply_inplace_schur(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const ProxyCube<T1>& P = x.P;

  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, P.n_rows, P.n_cols, P.n_slices, "element-wise cube multiplication");

        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;

  if(is_cube_generator<eop_cube_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] *= eop_aux::generate<eT,eop_cube_type>();
      }
    }
  else
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] *= eop_cube_core<eop_cube_type>::process(x, P[i]);
      }
    }
  }



template<typename eop_cube_type>
template<typename T1>
arma_hot
inline
void
eop_cube_core<eop_cube_type>::apply_inplace_div(Cube<typename T1::elem_type>& out, const eOpCube<T1, eop_cube_type>& x)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const ProxyCube<T1>& P = x.P;

  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, P.n_rows, P.n_cols, P.n_slices, "element-wise cube division");

        eT* out_mem = out.memptr();
  const u32 n_elem  = P.n_elem;

  if(is_cube_generator<eop_cube_type>::value == true)
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] /= eop_aux::generate<eT,eop_cube_type>();
      }
    }
  else
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_mem[i] /= eop_cube_core<eop_cube_type>::process(x, P[i]);
      }
    }
  }



//! @}
