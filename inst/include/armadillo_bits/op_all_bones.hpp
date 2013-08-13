// Copyright (C) 2013 Conrad Sanderson
// Copyright (C) 2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_all
//! @{



class op_all
  {
  public:
  
  
  template<typename T1>
  static inline bool
  all_vec_helper(const Base<typename T1::elem_type, T1>& X);
  
  
  template<typename T1, typename op_type>
  static inline bool
  all_vec_helper
    (
    const mtOp<uword, T1, op_type>& X,
    const typename arma_op_rel_only<op_type>::result           junk1 = 0,
    const typename arma_not_cx<typename T1::elem_type>::result junk2 = 0
    );
  
  
  template<typename T1, typename T2, typename glue_type>
  static inline bool
  all_vec_helper
    (
    const mtGlue<uword, T1, T2, glue_type>& X,
    const typename arma_glue_rel_only<glue_type>::result       junk1 = 0,
    const typename arma_not_cx<typename T1::elem_type>::result junk2 = 0,
    const typename arma_not_cx<typename T2::elem_type>::result junk3 = 0
    );
  
  
  template<typename T1>
  static inline bool all_vec(T1& X);
  
  
  template<typename T1>
  static inline void apply_helper(Mat<uword>& out, const Proxy<T1>& P, const uword dim);
  
  
  template<typename T1>
  static inline void apply(Mat<uword>& out, const mtOp<uword, T1, op_all>& X);
  };



//! @}
