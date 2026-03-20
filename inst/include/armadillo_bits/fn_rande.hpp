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


//! \addtogroup fn_rande
//! @{



template<typename obj_type>
arma_warn_unused
inline
obj_type
rande(const uword n_rows, const uword n_cols, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename obj_type::elem_type eT;
  
  if(is_Col<obj_type>::value)
    {
    arma_conform_check( (n_cols != 1), "rande(): incompatible size" );
    }
  else
  if(is_Row<obj_type>::value)
    {
    arma_conform_check( (n_rows != 1), "rande(): incompatible size" );
    }
  
  double lambda = double(1);
  double unused = double(0);
  
  param.get_double_vals(lambda, unused);
  
  arma_conform_check( (lambda <= double(0)), "rande(): incorrect distribution parameters; lambda must be greater than zero" );
  
  obj_type out(n_rows, n_cols, arma_nozeros_indicator());
  
  arma_rng::rande<eT>::fill(out.memptr(), out.n_elem, lambda);
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
rande(const SizeMat& s, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  return rande<obj_type>(s.n_rows, s.n_cols, param);
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
rande(const uword n_elem, const distr_param& param = distr_param(), const arma_empty_class junk1 = arma_empty_class(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk2 = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  const uword n_rows = (is_Row<obj_type>::value) ? uword(1) : n_elem;
  const uword n_cols = (is_Row<obj_type>::value) ? n_elem   : uword(1);
  
  return rande<obj_type>(n_rows, n_cols, param);
  }



arma_warn_unused
inline
mat
rande(const uword n_rows, const uword n_cols, const distr_param& param = distr_param())
  {
  arma_debug_sigprint();
  
  return rande<mat>(n_rows, n_cols, param);
  }



arma_warn_unused
inline
mat
rande(const SizeMat& s, const distr_param& param = distr_param())
  {
  arma_debug_sigprint();
  
  return rande<mat>(s.n_rows, s.n_cols, param);
  }



arma_warn_unused
inline
vec
rande(const uword n_elem, const distr_param& param = distr_param())
  {
  arma_debug_sigprint();
  
  return rande<vec>(n_elem, uword(1), param);
  }



arma_warn_unused
inline
double
rande(const distr_param& param = distr_param())
  {
  arma_debug_sigprint();
  
  double lambda = double(1);
  double unused = double(0);
  
  param.get_double_vals(lambda, unused);
  
  arma_conform_check( (lambda <= double(0)), "rande(): incorrect distribution parameters; lambda must be greater than zero" );
  
  double out_val = double(0);
  
  arma_rng::rande<double>::fill(&out_val, uword(1), lambda);
  
  return out_val;
  }



template<typename eT>
arma_warn_unused
inline
typename arma_real_or_cx_only<eT>::result
rande(const distr_param& param = distr_param())
  {
  arma_debug_sigprint();
  
  double lambda = double(1);
  double unused = double(0);
  
  param.get_double_vals(lambda, unused);
  
  arma_conform_check( (lambda <= double(0)), "rande(): incorrect distribution parameters; lambda must be greater than zero" );
  
  eT out_val = eT(0);
  
  arma_rng::rande<eT>::fill(&out_val, uword(1), lambda);
  
  return out_val;
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
rande(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename cube_type::elem_type eT;
  
  double lambda = double(1);
  double unused = double(0);
  
  param.get_double_vals(lambda, unused);
  
  arma_conform_check( (lambda <= double(0)), "rande(): incorrect distribution parameters; lambda must be greater than zero" );
  
  cube_type out(n_rows, n_cols, n_slices, arma_nozeros_indicator());
  
  arma_rng::rande<eT>::fill(out.memptr(), out.n_elem, lambda);
  
  return out;
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
rande(const SizeCube& s, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_debug_sigprint();
  arma_ignore(junk);
  
  return rande<cube_type>(s.n_rows, s.n_cols, s.n_slices, param);
  }



arma_warn_unused
inline
cube
rande(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param())
  {
  arma_debug_sigprint();
  
  return rande<cube>(n_rows, n_cols, n_slices, param);
  }



arma_warn_unused
inline
cube
rande(const SizeCube& s, const distr_param& param = distr_param())
  {
  arma_debug_sigprint();
  
  return rande<cube>(s.n_rows, s.n_cols, s.n_slices, param);
  }



//! @}
