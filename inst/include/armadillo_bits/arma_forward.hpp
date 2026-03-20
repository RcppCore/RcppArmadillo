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


using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::size_t;

template<typename elem_type, typename derived> struct Base;
template<typename elem_type, typename derived> struct BaseCube;

template<typename eT> class Mat;
template<typename eT> class Col;
template<typename eT> class Row;
template<typename eT> class Cube;
template<typename oT> class field;

template<typename eT>               struct xvec_htrans;
template<typename eT, bool do_conj> struct xtrans_mat;

template<typename eT> class subview;
template<typename eT> class subview_col;
template<typename eT> class subview_cols;
template<typename eT> class subview_row;
template<typename eT> class subview_row_strans;
template<typename eT> class subview_row_htrans;
template<typename eT> class subview_cube;
template<typename oT> class subview_field;

template<typename eT> class SpValProxy;
template<typename eT> class SpMat;
template<typename eT> class SpCol;
template<typename eT> class SpRow;
template<typename eT> class SpSubview;
template<typename eT> class SpSubview_col;
template<typename eT> class SpSubview_row;

template<typename eT> class diagview;
template<typename eT> class spdiagview;

template<typename eT> class MapMat;
template<typename eT> class MapMat_val;
template<typename eT> class SpMat_MapMat_val;
template<typename eT> class SpSubview_MapMat_val;

template<typename eT, typename T1>              class subview_elem1;
template<typename eT, typename T1, typename T2> class subview_elem2;

template<typename parent, unsigned int mode>              class subview_each1;
template<typename parent, unsigned int mode, typename TB> class subview_each2;

template<typename eT>              class subview_cube_each1;
template<typename eT, typename TB> class subview_cube_each2;
template<typename eT, typename T1> class subview_cube_slices;

template<typename eT, typename T1> class SpSubview_col_list;


template<typename eT> class   Mat_noalias;
template<typename eT> class SpMat_noalias;


struct SizeMat;
struct SizeCube;

struct arma_empty_class {};

struct diskio;

struct op_strans;
struct op_htrans;
struct op_htrans2;
struct op_inv_gen_default;
struct op_inv_spd_default;
struct op_inv_gen_full;
struct op_inv_spd_full;
struct op_diagmat;
struct op_trimat;
struct op_vectorise_row;
struct op_vectorise_col;
struct op_symmatu;
struct op_symmatl;
struct op_omit;

struct op_row_as_mat;
struct op_col_as_mat;

struct glue_times;
struct glue_times_diag;

struct glue_rel_lt;
struct glue_rel_gt;
struct glue_rel_lteq;
struct glue_rel_gteq;
struct glue_rel_eq;
struct glue_rel_noteq;
struct glue_rel_and;
struct glue_rel_or;

struct op_rel_lt_pre;
struct op_rel_lt_post;
struct op_rel_gt_pre;
struct op_rel_gt_post;
struct op_rel_lteq_pre;
struct op_rel_lteq_post;
struct op_rel_gteq_pre;
struct op_rel_gteq_post;
struct op_rel_eq;
struct op_rel_noteq;

struct gen_eye;
struct gen_ones;
struct gen_zeros;



struct spop_strans;
struct spop_htrans;
struct spop_vectorise_row;
struct spop_vectorise_col;
struct spop_square;

struct spop_rel_lt_pre;
struct spop_rel_lt_post;
struct spop_rel_gt_pre;
struct spop_rel_gt_post;
struct spop_rel_lteq_pre;
struct spop_rel_lteq_post;
struct spop_rel_gteq_pre;
struct spop_rel_gteq_post;
struct spop_rel_eq;
struct spop_rel_noteq;

struct spglue_plus;
struct spglue_minus;
struct spglue_schur;
struct spglue_times;
struct spglue_max;
struct spglue_min;
struct spglue_rel_lt;
struct spglue_rel_gt;

struct op_sp_as_dense;

struct op_internal_equ;
struct op_internal_plus;
struct op_internal_minus;
struct op_internal_schur;
struct op_internal_div;



struct traits_op_default
  {
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = false;
    static constexpr bool is_col  = false;
    static constexpr bool is_xvec = false;
    };
  };


struct traits_op_xvec
  {
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = false;
    static constexpr bool is_col  = false;
    static constexpr bool is_xvec = true;
    };
  };


struct traits_op_col
  {
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = false;
    static constexpr bool is_col  = true;
    static constexpr bool is_xvec = false;
    };
  };


struct traits_op_row
  {
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = true;
    static constexpr bool is_col  = false;
    static constexpr bool is_xvec = false;
    };
  };


struct traits_op_passthru
  {
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = T1::is_row;
    static constexpr bool is_col  = T1::is_col;
    static constexpr bool is_xvec = T1::is_xvec;
    };
  };


struct traits_glue_default
  {
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = false;
    static constexpr bool is_col  = false;
    static constexpr bool is_xvec = false;
    };
  };


struct traits_glue_or
  {
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = (T1::is_row  || T2::is_row );
    static constexpr bool is_col  = (T1::is_col  || T2::is_col );
    static constexpr bool is_xvec = (T1::is_xvec || T2::is_xvec);
    };
  };



template<const bool, const bool, const bool, const bool> struct gemm;
template<const bool, const bool, const bool>             struct gemv;


template<                 typename eT, typename gen_type> struct         Gen; 

template<                 typename T1, typename  op_type> struct          Op; 
template<                 typename T1, typename eop_type> struct         eOp;
template<                 typename T1, typename  op_type> struct     SpToDOp; 
template<                 typename T1, typename  op_type> struct CubeToMatOp;
template<typename out_eT, typename T1, typename  op_type> struct        mtOp;

template<                 typename T1, typename T2, typename  glue_type> struct      Glue;
template<                 typename T1, typename T2, typename eglue_type> struct     eGlue;
template<                 typename T1, typename T2, typename  glue_type> struct SpToDGlue;
template<typename out_eT, typename T1, typename T2, typename  glue_type> struct    mtGlue;



template<                 typename eT, typename gen_type> struct  GenCube; 

template<                 typename T1, typename  op_type> struct   OpCube; 
template<                 typename T1, typename eop_type> struct  eOpCube; 
template<typename out_eT, typename T1, typename  op_type> struct mtOpCube;

template<                 typename T1, typename T2, typename  glue_type> struct   GlueCube;
template<                 typename T1, typename T2, typename eglue_type> struct  eGlueCube;
template<typename out_eT, typename T1, typename T2, typename  glue_type> struct mtGlueCube;


template<typename T1> struct Proxy;
template<typename T1> struct ProxyCube;

template<typename T1> struct diagmat_proxy;

template<typename T1> struct unwrap;
template<typename T1> struct quasi_unwrap;
template<typename T1> struct unwrap_cube;
template<typename T1> struct unwrap_spmat;




struct state_type
  {
  #if   defined(ARMA_USE_OPENMP)
                int  state;
  #elif defined(ARMA_USE_STD_MUTEX)
    std::atomic<int> state;
  #else
                int  state;
  #endif
  
  arma_inline state_type() : state(int(0)) {}
  
  // openmp: "omp atomic" does an implicit flush on the affected variable
  // C++11:  std::atomic<>::load() and std::atomic<>::store() use std::memory_order_seq_cst by default, which has an implied fence
  
  arma_inline
  operator int () const
    {
    int out;
    
    #if   defined(ARMA_USE_OPENMP)
      #pragma omp atomic read seq_cst
      out = state;
    #elif defined(ARMA_USE_STD_MUTEX)
      out = state.load();
    #else
      out = state;
    #endif
    
    return out;
    }
  
  arma_inline
  void
  operator= (const int in_state)
    {
    #if   defined(ARMA_USE_OPENMP)
      #pragma omp atomic write seq_cst
      state = in_state;
    #elif defined(ARMA_USE_STD_MUTEX)
      state.store(in_state);
    #else
      state = in_state;
    #endif
    }
  };


template<                 typename T1, typename spop_type> struct   SpOp;
template<typename out_eT, typename T1, typename spop_type> struct mtSpOp;
template<typename out_eT, typename T1, typename   op_type> struct mtSpReduceOp;

template<                 typename T1, typename T2, typename spglue_type> struct   SpGlue;
template<typename out_eT, typename T1, typename T2, typename spglue_type> struct mtSpGlue;


template<typename T1> struct SpProxy;



struct arma_vec_indicator     {};
struct arma_fixed_indicator   {};
struct arma_reserve_indicator {};
struct arma_layout_indicator  {};

template<bool do_zeros> struct arma_initmode_indicator {};

struct arma_zeros_indicator   : public arma_initmode_indicator<true > {};
struct arma_nozeros_indicator : public arma_initmode_indicator<false> {};



//! \addtogroup injector
//! @{

template<typename Dummy = int> struct injector_end_of_row {};

// DEPRECATED: DO NOT USE IN NEW CODE
static const injector_end_of_row<> endr = injector_end_of_row<>();
//!< endr indicates "end of row" when using the << operator;
//!< similar conceptual meaning to std::endl

//! @}



//! \addtogroup diskio
//! @{


enum struct file_type : unsigned int
  {
  file_type_unknown,
  auto_detect,        //!< attempt to automatically detect the file type
  raw_ascii,          //!< raw text (ASCII), without a header
  arma_ascii,         //!< Armadillo text format, with a header specifying matrix type and size
  csv_ascii,          //!< comma separated values (CSV), without a header
  raw_binary,         //!< raw binary format (machine dependent), without a header
  arma_binary,        //!< Armadillo binary format (machine dependent), with a header specifying matrix type and size
  pgm_binary,         //!< Portable Grey Map (greyscale image)
  ppm_binary,         //!< Portable Pixel Map (colour image), used by the field and cube classes
  hdf5_binary,        //!< HDF5: open binary format, not specific to Armadillo, which can store arbitrary data
  hdf5_binary_trans,  //!< [NOTE: DO NOT USE - deprecated] as per hdf5_binary, but save/load the data with columns transposed to rows
  coord_ascii,        //!< simple co-ordinate format for sparse matrices (indices start at zero)
  ssv_ascii,          //!< similar to csv_ascii; uses semicolon (;) instead of comma (,) as the separator
  };


static constexpr file_type file_type_unknown  = file_type::file_type_unknown;
static constexpr file_type auto_detect        = file_type::auto_detect;
static constexpr file_type raw_ascii          = file_type::raw_ascii;
static constexpr file_type arma_ascii         = file_type::arma_ascii;
static constexpr file_type csv_ascii          = file_type::csv_ascii;
static constexpr file_type raw_binary         = file_type::raw_binary;
static constexpr file_type arma_binary        = file_type::arma_binary;
static constexpr file_type pgm_binary         = file_type::pgm_binary;
static constexpr file_type ppm_binary         = file_type::ppm_binary;
static constexpr file_type hdf5_binary        = file_type::hdf5_binary;
static constexpr file_type hdf5_binary_trans  = file_type::hdf5_binary_trans;
static constexpr file_type coord_ascii        = file_type::coord_ascii;
static constexpr file_type ssv_ascii          = file_type::ssv_ascii;


struct hdf5_name;
struct  csv_name;


//! @}



//! \addtogroup fn_spsolve
//! @{


struct spsolve_opts_base
  {
  const unsigned int id;
  
  inline spsolve_opts_base(const unsigned int in_id) : id(in_id) {}
  };


struct spsolve_opts_none : public spsolve_opts_base
  {
  inline spsolve_opts_none() : spsolve_opts_base(0) {}
  };


struct superlu_opts : public spsolve_opts_base
  {
  typedef enum {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD} permutation_type;
  
  typedef enum {REF_NONE, REF_SINGLE, REF_DOUBLE, REF_EXTRA} refine_type;
  
  bool             allow_ugly;
  bool             equilibrate;
  bool             symmetric;
  double           pivot_thresh;
  permutation_type permutation;
  refine_type      refine;
  
  inline superlu_opts()
    : spsolve_opts_base(1)
    {
    allow_ugly   = false;
    equilibrate  = false;
    symmetric    = false;
    pivot_thresh = 1.0;
    permutation  = COLAMD;
    refine       = REF_NONE;
    }
  };


//! @}



//! \ingroup fn_eigs_sym fs_eigs_gen
//! @{


struct eigs_opts
  {
  double       tol;     // tolerance
  unsigned int maxiter; // max iterations
  unsigned int subdim;  // subspace dimension
  
  inline eigs_opts()
    {
    tol     = 0.0;
    maxiter = 1000;
    subdim  = 0;
    }
  };


//! @}
