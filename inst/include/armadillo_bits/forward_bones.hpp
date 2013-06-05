// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2013 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


using std::cout;
using std::cerr;
using std::endl;
using std::ios;

template<typename elem_type, typename derived> struct Base;
template<typename elem_type, typename derived> struct BaseCube;

template<typename eT> class Mat;
template<typename eT> class Col;
template<typename eT> class Row;
template<typename eT> class Cube;
template<typename eT> class xvec_htrans;
template<typename oT> class field;

template<typename eT> class subview;
template<typename eT> class subview_col;
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

template<typename eT> class diagview;

template<typename eT, typename T1>              class subview_elem1;
template<typename eT, typename T1, typename T2> class subview_elem2;

template<typename parent, unsigned int mode>              class subview_each1;
template<typename parent, unsigned int mode, typename TB> class subview_each2;


class arma_empty_class {};

class diskio;

class op_min;
class op_max;

class op_strans;
class op_htrans;
class op_htrans2;
class op_inv;
class op_sum;
class op_abs;
class op_diagmat;
class op_trimat;
class op_diagvec;

class eop_conj;

class glue_times;
class glue_times_diag;

class glue_rel_lt;
class glue_rel_gt;
class glue_rel_lteq;
class glue_rel_gteq;
class glue_rel_eq;
class glue_rel_noteq;

class op_rel_lt_pre;
class op_rel_lt_post;
class op_rel_gt_pre;
class op_rel_gt_post;
class op_rel_lteq_pre;
class op_rel_lteq_post;
class op_rel_gteq_pre;
class op_rel_gteq_post;
class op_rel_eq;
class op_rel_noteq;

class gen_ones_diag;
class gen_ones_full;
class gen_zeros;
class gen_randu;
class gen_randn;

class glue_mixed_plus;
class glue_mixed_minus;
class glue_mixed_div;
class glue_mixed_schur;
class glue_mixed_times;

class op_cx_scalar_times;
class op_cx_scalar_plus;
class op_cx_scalar_minus_pre;
class op_cx_scalar_minus_post;
class op_cx_scalar_div_pre;
class op_cx_scalar_div_post;



class op_subview_elem_equ;
class op_subview_elem_inplace_plus;
class op_subview_elem_inplace_minus;
class op_subview_elem_inplace_schur;
class op_subview_elem_inplace_div;



template<const bool, const bool, const bool, const bool> class gemm;
template<const bool, const bool, const bool>             class gemv;


template<                 typename eT, typename gen_type> class  Gen; 

template<                 typename T1, typename  op_type> class   Op; 
template<                 typename T1, typename eop_type> class  eOp;
template<typename out_eT, typename T1, typename  op_type> class mtOp;

template<                 typename T1, typename T2, typename  glue_type> class   Glue;
template<                 typename T1, typename T2, typename eglue_type> class  eGlue;
template<typename out_eT, typename T1, typename T2, typename  glue_type> class mtGlue;



template<                 typename eT, typename gen_type> class  GenCube; 

template<                 typename T1, typename  op_type> class   OpCube; 
template<                 typename T1, typename eop_type> class  eOpCube; 
template<typename out_eT, typename T1, typename  op_type> class mtOpCube;

template<                 typename T1, typename T2, typename  glue_type> class   GlueCube;
template<                 typename T1, typename T2, typename eglue_type> class  eGlueCube;
template<typename out_eT, typename T1, typename T2, typename  glue_type> class mtGlueCube;


template<typename T1> class Proxy;
template<typename T1> class ProxyCube;



class spop_strans;
class spop_htrans;
class spop_scalar_times;

class spglue_plus;
class spglue_plus2;

class spglue_minus;
class spglue_minus2;

class spglue_times;
class spglue_times2;


template<                 typename T1, typename spop_type> class   SpOp;
template<typename out_eT, typename T1, typename spop_type> class mtSpOp;

template<typename T1, typename T2, typename spglue_type> class SpGlue;


template<typename T1> class SpProxy;



struct arma_vec_indicator   {};
struct arma_fixed_indicator {};


//! \addtogroup injector
//! @{

template<typename Dummy = int> struct injector_end_of_row {};

static const injector_end_of_row<> endr = injector_end_of_row<>();
//!< endr indicates "end of row" when using the << operator;
//!< similar conceptual meaning to std::endl

//! @}



//! \addtogroup diskio
//! @{


enum file_type
  {
  file_type_unknown,
  auto_detect,  //!< Automatically detect the file type
  raw_ascii,    //!< ASCII format (text), without any other information.
  arma_ascii,   //!< Armadillo ASCII format (text), with information about matrix type and size
  csv_ascii,    //!< comma separated values (CSV), without any other information
  raw_binary,   //!< raw binary format, without any other information.
  arma_binary,  //!< Armadillo binary format, with information about matrix type and size
  pgm_binary,   //!< Portable Grey Map (greyscale image)
  ppm_binary,   //!< Portable Pixel Map (colour image), used by the field and cube classes
  hdf5_binary,  //!< Open binary format, not specific to Armadillo, which can store arbitrary data
  coord_ascii   //!< simple co-ordinate format for sparse matrices
  };


//! @}


