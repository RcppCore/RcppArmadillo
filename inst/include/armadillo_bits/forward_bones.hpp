// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


using std::cout;
using std::cerr;
using std::endl;
using std::ios;

template<typename eT> class Mat;
template<typename eT> class Col;
template<typename eT> class Row;
template<typename eT> class Cube;
template<typename oT> class field;

template<typename eT> class subview;
template<typename eT> class subview_col;
template<typename eT> class subview_row;
template<typename eT> class subview_cube;
template<typename oT> class subview_field;

template<typename eT> class diagview;

template<typename eT, typename T1>              class subview_elem1;
template<typename eT, typename T1, typename T2> class subview_elem2;


class arma_empty_class {};

class diskio;

class op_min;
class op_max;

class op_strans;
class op_htrans;
class op_inv;
class op_sum;
class op_abs;
class op_diagmat;
class op_trimat;

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



//! \addtogroup injector
//! @{


struct injector_end_of_row {};

static const injector_end_of_row endr = injector_end_of_row();
//!< endr indicates "end of row" when using the << operator;
//!< similar conceptual meaning to std::endl

//! @}



//! \addtogroup diskio
//! @{


enum file_type
  {
  file_type_unknown,
  auto_detect,  //!< Automatically detect the file type (file must be one of the following types)
  raw_ascii,    //!< ASCII format (text), without any other information.
  arma_ascii,   //!< Armadillo ASCII format (text), with information about matrix type and size
  csv_ascii,    //!< comma separated values (CSV), without any other information
  raw_binary,   //!< raw binary format, without any other information.
  arma_binary,  //!< Armadillo binary format, with information about matrix type and size
  pgm_binary,   //!< Portable Grey Map (greyscale image)
  ppm_binary    //!< Portable Pixel Map (colour image), used by the field and cube classes
  };


//! @}


