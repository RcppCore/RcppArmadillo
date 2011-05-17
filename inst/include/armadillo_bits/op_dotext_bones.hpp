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


//! \addtogroup op_dotext
//! @{



class op_dotext
  {
  public:
  
  
  template<typename eT>
  inline static eT direct_rowvec_mat_colvec       (const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_transmat_colvec  (const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_diagmat_colvec   (const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_invdiagmat_colvec(const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_invdiagvec_colvec(const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  };



//! @}

