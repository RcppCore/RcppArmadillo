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


#ifndef ARMA_INCLUDES
#define ARMA_INCLUDES


#include <cstdlib>
#include <cstring>
#include <climits>
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <complex>
#include <vector>


#include "armadillo_bits/config.hpp"
#include "armadillo_bits/compiler_setup.hpp"
#include "armadillo_bits/undefine_conflicts.hpp"


#if defined(ARMA_USE_CXX11)
  #include <initializer_list>
#endif


#if defined(ARMA_HAVE_GETTIMEOFDAY)
  #include <sys/time.h>
  #undef ARMA_USE_BOOST_DATE
#endif


#if defined(ARMA_USE_BOOST_DATE)
  #include <boost/date_time/posix_time/posix_time.hpp>
#endif


#if defined(ARMA_HAVE_STD_TR1)
  // TODO: add handling of this functionality when use of C++11 is enabled
  #include <tr1/cmath>
  #include <tr1/complex>
#elif defined(ARMA_USE_BOOST)
  #include <boost/math/complex.hpp>
  #include <boost/math/special_functions/acosh.hpp>
  #include <boost/math/special_functions/asinh.hpp>
  #include <boost/math/special_functions/atanh.hpp>
#endif


#if defined(ARMA_USE_BOOST)
  #if defined(ARMA_EXTRA_DEBUG)
    #include <boost/format.hpp>
    #include <boost/current_function.hpp>
    #define ARMA_USE_BOOST_FORMAT
  #endif
#endif


#include "armadillo_bits/include_atlas.hpp"
#include "armadillo_bits/itpp_wrap.hpp"


//! \namespace arma namespace for Armadillo classes and functions
namespace arma
  {
  
  // preliminaries
  
  #include "armadillo_bits/forward_bones.hpp"
  #include "armadillo_bits/arma_static_check.hpp"
  #include "armadillo_bits/typedef.hpp"
  #include "armadillo_bits/typedef_blas_int.hpp"
  #include "armadillo_bits/format_wrap.hpp"
  #include "armadillo_bits/arma_version.hpp"
  #include "armadillo_bits/arma_config.hpp"
  #include "armadillo_bits/traits.hpp"
  #include "armadillo_bits/promote_type.hpp"
  #include "armadillo_bits/upgrade_val.hpp"
  #include "armadillo_bits/restrictors.hpp"
  #include "armadillo_bits/access.hpp"
  #include "armadillo_bits/span.hpp"
  #include "armadillo_bits/constants.hpp"
  
  
  //
  // class prototypes
  
  #include "armadillo_bits/Base_bones.hpp"
  #include "armadillo_bits/BaseCube_bones.hpp"
  
  #include "armadillo_bits/blas_bones.hpp"
  #include "armadillo_bits/lapack_bones.hpp"
  #include "armadillo_bits/atlas_bones.hpp"
  
  #include "armadillo_bits/blas_wrapper.hpp"
  #include "armadillo_bits/lapack_wrapper.hpp"
  #include "armadillo_bits/atlas_wrapper.hpp"
  
  #include "armadillo_bits/arrayops_bones.hpp"
  #include "armadillo_bits/podarray_bones.hpp"
  #include "armadillo_bits/auxlib_bones.hpp"
  
  #include "armadillo_bits/injector_bones.hpp"
  
  #include "armadillo_bits/Mat_bones.hpp"
  #include "armadillo_bits/Col_bones.hpp"
  #include "armadillo_bits/Row_bones.hpp"
  #include "armadillo_bits/Cube_bones.hpp"
  
  #include "armadillo_bits/typedef_fixed.hpp"
  
  #include "armadillo_bits/field_bones.hpp"
  #include "armadillo_bits/subview_bones.hpp"
  #include "armadillo_bits/subview_elem1_bones.hpp"
  #include "armadillo_bits/subview_field_bones.hpp"
  #include "armadillo_bits/subview_cube_bones.hpp"
  #include "armadillo_bits/diagview_bones.hpp"
  
  
  #include "armadillo_bits/diskio_bones.hpp"
  #include "armadillo_bits/wall_clock_bones.hpp"
  #include "armadillo_bits/running_stat_bones.hpp"
  #include "armadillo_bits/running_stat_vec_bones.hpp"
  
  #include "armadillo_bits/Op_bones.hpp"
  #include "armadillo_bits/OpCube_bones.hpp"
  
  #include "armadillo_bits/eOp_bones.hpp"
  #include "armadillo_bits/eOpCube_bones.hpp"
  
  #include "armadillo_bits/mtOp_bones.hpp"
  #include "armadillo_bits/mtOpCube_bones.hpp"
  
  #include "armadillo_bits/Glue_bones.hpp"
  #include "armadillo_bits/eGlue_bones.hpp"
  #include "armadillo_bits/mtGlue_bones.hpp"
  
  #include "armadillo_bits/GlueCube_bones.hpp"
  #include "armadillo_bits/eGlueCube_bones.hpp"
  #include "armadillo_bits/mtGlueCube_bones.hpp"
  
  #include "armadillo_bits/eop_core_bones.hpp"
  #include "armadillo_bits/eglue_core_bones.hpp"
  
  #include "armadillo_bits/Gen_bones.hpp"
  #include "armadillo_bits/GenCube_bones.hpp"
  
  #include "armadillo_bits/op_diagmat_bones.hpp"
  #include "armadillo_bits/op_diagvec_bones.hpp"
  #include "armadillo_bits/op_dot_bones.hpp"
  #include "armadillo_bits/op_inv_bones.hpp"
  #include "armadillo_bits/op_htrans_bones.hpp"
  #include "armadillo_bits/op_max_bones.hpp"
  #include "armadillo_bits/op_min_bones.hpp"
  #include "armadillo_bits/op_mean_bones.hpp"
  #include "armadillo_bits/op_median_bones.hpp"
  #include "armadillo_bits/op_sort_bones.hpp"
  #include "armadillo_bits/op_sum_bones.hpp"
  #include "armadillo_bits/op_stddev_bones.hpp"
  #include "armadillo_bits/op_strans_bones.hpp"
  #include "armadillo_bits/op_var_bones.hpp"
  #include "armadillo_bits/op_repmat_bones.hpp"
  #include "armadillo_bits/op_reshape_bones.hpp"
  #include "armadillo_bits/op_resize_bones.hpp"
  #include "armadillo_bits/op_cov_bones.hpp"
  #include "armadillo_bits/op_cor_bones.hpp"
  #include "armadillo_bits/op_shuffle_bones.hpp"
  #include "armadillo_bits/op_prod_bones.hpp"
  #include "armadillo_bits/op_pinv_bones.hpp"
  #include "armadillo_bits/op_dotext_bones.hpp"
  #include "armadillo_bits/op_flip_bones.hpp"
  #include "armadillo_bits/op_princomp_bones.hpp"
  #include "armadillo_bits/op_misc_bones.hpp"
  #include "armadillo_bits/op_relational_bones.hpp"
  #include "armadillo_bits/op_find_bones.hpp"
  #include "armadillo_bits/op_chol_bones.hpp"
  #include "armadillo_bits/op_cx_scalar_bones.hpp"
  #include "armadillo_bits/op_trimat_bones.hpp"
  #include "armadillo_bits/op_cumsum_bones.hpp"
  #include "armadillo_bits/op_symmat_bones.hpp"
  
  #include "armadillo_bits/glue_times_bones.hpp"
  #include "armadillo_bits/glue_mixed_bones.hpp"
  #include "armadillo_bits/glue_cov_bones.hpp"
  #include "armadillo_bits/glue_cor_bones.hpp"
  #include "armadillo_bits/glue_kron_bones.hpp"
  #include "armadillo_bits/glue_cross_bones.hpp"
  #include "armadillo_bits/glue_join_bones.hpp"
  #include "armadillo_bits/glue_relational_bones.hpp"
  #include "armadillo_bits/glue_solve_bones.hpp"
  #include "armadillo_bits/glue_conv_bones.hpp"
  #include "armadillo_bits/glue_toeplitz_bones.hpp"
  
  //
  // debugging functions
  
  #include "armadillo_bits/debug.hpp"
  
  //
  //
  
  #include "armadillo_bits/cmath_wrap.hpp"
  
  //
  // classes that underlay metaprogramming 
  
  #include "armadillo_bits/Proxy.hpp"
  #include "armadillo_bits/ProxyCube.hpp"
  
  #include "armadillo_bits/diagmat_proxy.hpp"

  #include "armadillo_bits/unwrap.hpp"
  #include "armadillo_bits/unwrap_cube.hpp"

  #include "armadillo_bits/strip.hpp"
  
  #include "armadillo_bits/Op_meat.hpp"
  #include "armadillo_bits/OpCube_meat.hpp"
  
  #include "armadillo_bits/mtOp_meat.hpp"
  #include "armadillo_bits/mtOpCube_meat.hpp"
  
  #include "armadillo_bits/Glue_meat.hpp"
  #include "armadillo_bits/GlueCube_meat.hpp"
  
  #include "armadillo_bits/eop_aux.hpp"
  
  #include "armadillo_bits/eOp_meat.hpp"
  #include "armadillo_bits/eOpCube_meat.hpp"
  
  #include "armadillo_bits/eGlue_meat.hpp"
  #include "armadillo_bits/eGlueCube_meat.hpp"

  #include "armadillo_bits/mtGlue_meat.hpp"
  #include "armadillo_bits/mtGlueCube_meat.hpp"
  
  #include "armadillo_bits/Base_meat.hpp"
  #include "armadillo_bits/BaseCube_meat.hpp"
  
  #include "armadillo_bits/Gen_meat.hpp"
  #include "armadillo_bits/GenCube_meat.hpp"
  
  
  //
  // ostream
  
  #include "armadillo_bits/arma_ostream_bones.hpp"
  #include "armadillo_bits/arma_ostream_meat.hpp"
  
  
  //
  // operators
  
  #include "armadillo_bits/operator_plus.hpp"
  #include "armadillo_bits/operator_minus.hpp"
  #include "armadillo_bits/operator_times.hpp"
  #include "armadillo_bits/operator_schur.hpp"
  #include "armadillo_bits/operator_div.hpp"
  #include "armadillo_bits/operator_relational.hpp"
  
  #include "armadillo_bits/operator_cube_plus.hpp"
  #include "armadillo_bits/operator_cube_minus.hpp"
  #include "armadillo_bits/operator_cube_times.hpp"
  #include "armadillo_bits/operator_cube_schur.hpp"
  #include "armadillo_bits/operator_cube_div.hpp"
  #include "armadillo_bits/operator_cube_relational.hpp"
  
  #include "armadillo_bits/operator_ostream.hpp"
  
  
  //
  // user accessible functions
  
  // the order of the fn_*.hpp include files matters,
  // as some files require functionality given in preceding files
  
  #include "armadillo_bits/fn_conv_to.hpp"
  #include "armadillo_bits/fn_min.hpp"
  #include "armadillo_bits/fn_max.hpp"
  #include "armadillo_bits/fn_accu.hpp"
  #include "armadillo_bits/fn_sum.hpp"
  #include "armadillo_bits/fn_diagmat.hpp"
  #include "armadillo_bits/fn_diagvec.hpp"
  #include "armadillo_bits/fn_inv.hpp"
  #include "armadillo_bits/fn_trace.hpp"
  #include "armadillo_bits/fn_trans.hpp"
  #include "armadillo_bits/fn_det.hpp"
  #include "armadillo_bits/fn_log_det.hpp"
  #include "armadillo_bits/fn_eig.hpp"
  #include "armadillo_bits/fn_lu.hpp"
  #include "armadillo_bits/fn_zeros.hpp"
  #include "armadillo_bits/fn_ones.hpp"
  #include "armadillo_bits/fn_eye.hpp"
  #include "armadillo_bits/fn_misc.hpp"
  #include "armadillo_bits/fn_elem.hpp"
  #include "armadillo_bits/fn_norm.hpp"
  #include "armadillo_bits/fn_dot.hpp"
  #include "armadillo_bits/fn_randu.hpp"
  #include "armadillo_bits/fn_randn.hpp"
  #include "armadillo_bits/fn_trig.hpp"
  #include "armadillo_bits/fn_mean.hpp"
  #include "armadillo_bits/fn_median.hpp"
  #include "armadillo_bits/fn_stddev.hpp"
  #include "armadillo_bits/fn_var.hpp"
  #include "armadillo_bits/fn_sort.hpp"
  #include "armadillo_bits/fn_sort_index.hpp"
  #include "armadillo_bits/fn_strans.hpp"
  #include "armadillo_bits/fn_chol.hpp"
  #include "armadillo_bits/fn_qr.hpp"
  #include "armadillo_bits/fn_svd.hpp"
  #include "armadillo_bits/fn_solve.hpp"
  #include "armadillo_bits/fn_repmat.hpp"
  #include "armadillo_bits/fn_reshape.hpp"
  #include "armadillo_bits/fn_resize.hpp"
  #include "armadillo_bits/fn_cov.hpp"
  #include "armadillo_bits/fn_cor.hpp"
  #include "armadillo_bits/fn_shuffle.hpp"
  #include "armadillo_bits/fn_prod.hpp"
  #include "armadillo_bits/fn_eps.hpp"
  #include "armadillo_bits/fn_pinv.hpp"
  #include "armadillo_bits/fn_rank.hpp"
  #include "armadillo_bits/fn_kron.hpp"
  #include "armadillo_bits/fn_flip.hpp"
  #include "armadillo_bits/fn_as_scalar.hpp"
  #include "armadillo_bits/fn_princomp.hpp"
  #include "armadillo_bits/fn_cross.hpp"
  #include "armadillo_bits/fn_join.hpp"
  #include "armadillo_bits/fn_conv.hpp"
  #include "armadillo_bits/fn_trunc_exp.hpp"
  #include "armadillo_bits/fn_trunc_log.hpp"
  #include "armadillo_bits/fn_toeplitz.hpp"
  #include "armadillo_bits/fn_trimat.hpp"
  #include "armadillo_bits/fn_cumsum.hpp"
  #include "armadillo_bits/fn_symmat.hpp"
  #include "armadillo_bits/fn_syl_lyap.hpp"
  
  //
  // class meat
  
  #include "armadillo_bits/gemv.hpp"
  #include "armadillo_bits/gemm.hpp"
  #include "armadillo_bits/gemm_mixed.hpp"
  
  #include "armadillo_bits/eop_core_meat.hpp"
  #include "armadillo_bits/eglue_core_meat.hpp"
  
  #include "armadillo_bits/arrayops_meat.hpp"
  #include "armadillo_bits/podarray_meat.hpp"
  #include "armadillo_bits/auxlib_meat.hpp"
  
  #include "armadillo_bits/injector_meat.hpp"
  
  #include "armadillo_bits/Mat_meat.hpp"
  #include "armadillo_bits/Col_meat.hpp"
  #include "armadillo_bits/Row_meat.hpp"
  #include "armadillo_bits/Cube_meat.hpp"
  #include "armadillo_bits/field_meat.hpp"
  #include "armadillo_bits/subview_meat.hpp"
  #include "armadillo_bits/subview_elem1_meat.hpp"
  #include "armadillo_bits/subview_field_meat.hpp"
  #include "armadillo_bits/subview_cube_meat.hpp"
  #include "armadillo_bits/diagview_meat.hpp"
  
  #include "armadillo_bits/diskio_meat.hpp"
  #include "armadillo_bits/wall_clock_meat.hpp"
  #include "armadillo_bits/running_stat_meat.hpp"
  #include "armadillo_bits/running_stat_vec_meat.hpp"
  
  #include "armadillo_bits/op_diagmat_meat.hpp"
  #include "armadillo_bits/op_diagvec_meat.hpp"
  #include "armadillo_bits/op_dot_meat.hpp"
  #include "armadillo_bits/op_inv_meat.hpp"
  #include "armadillo_bits/op_htrans_meat.hpp"
  #include "armadillo_bits/op_max_meat.hpp"
  #include "armadillo_bits/op_min_meat.hpp"
  #include "armadillo_bits/op_mean_meat.hpp"
  #include "armadillo_bits/op_median_meat.hpp"
  #include "armadillo_bits/op_sort_meat.hpp"
  #include "armadillo_bits/op_sum_meat.hpp"
  #include "armadillo_bits/op_stddev_meat.hpp"
  #include "armadillo_bits/op_strans_meat.hpp"
  #include "armadillo_bits/op_var_meat.hpp"
  #include "armadillo_bits/op_repmat_meat.hpp"
  #include "armadillo_bits/op_reshape_meat.hpp"
  #include "armadillo_bits/op_resize_meat.hpp"
  #include "armadillo_bits/op_cov_meat.hpp"
  #include "armadillo_bits/op_cor_meat.hpp"
  #include "armadillo_bits/op_shuffle_meat.hpp"
  #include "armadillo_bits/op_prod_meat.hpp"
  #include "armadillo_bits/op_pinv_meat.hpp"
  #include "armadillo_bits/op_dotext_meat.hpp"
  #include "armadillo_bits/op_flip_meat.hpp"
  #include "armadillo_bits/op_princomp_meat.hpp"
  #include "armadillo_bits/op_misc_meat.hpp"
  #include "armadillo_bits/op_relational_meat.hpp"
  #include "armadillo_bits/op_find_meat.hpp"
  #include "armadillo_bits/op_chol_meat.hpp"
  #include "armadillo_bits/op_cx_scalar_meat.hpp"
  #include "armadillo_bits/op_trimat_meat.hpp"
  #include "armadillo_bits/op_cumsum_meat.hpp"
  #include "armadillo_bits/op_symmat_meat.hpp"
  
  #include "armadillo_bits/glue_times_meat.hpp"
  #include "armadillo_bits/glue_mixed_meat.hpp"
  #include "armadillo_bits/glue_cov_meat.hpp"
  #include "armadillo_bits/glue_cor_meat.hpp"
  #include "armadillo_bits/glue_kron_meat.hpp"
  #include "armadillo_bits/glue_cross_meat.hpp"
  #include "armadillo_bits/glue_join_meat.hpp"
  #include "armadillo_bits/glue_relational_meat.hpp"
  #include "armadillo_bits/glue_solve_meat.hpp"
  #include "armadillo_bits/glue_conv_meat.hpp"
  #include "armadillo_bits/glue_toeplitz_meat.hpp"
  }
  
#endif

