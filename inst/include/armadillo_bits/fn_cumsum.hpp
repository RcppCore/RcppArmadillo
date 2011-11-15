// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_cumsum
//! @{



template<typename T1>
arma_inline
const Op<T1, op_cumsum_mat>
cumsum(const Base<typename T1::elem_type,T1>& X, const uword dim = 0)
  {
  arma_extra_debug_sigprint();

  return Op<T1, op_cumsum_mat>(X.get_ref(), dim, 0);
  }



template<typename eT>
arma_inline
const Op<Row<eT>, op_cumsum_vec>
cumsum(const Row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<Row<eT>, op_cumsum_vec>(A);
  }



template<typename eT>
arma_inline
const Op<Col<eT>, op_cumsum_vec>
cumsum(const Col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<Col<eT>, op_cumsum_vec>(A);
  }



template<typename eT>
arma_inline
const Op<subview_row<eT>, op_cumsum_vec>
cumsum(const subview_row<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<subview_row<eT>, op_cumsum_vec>(A);
  }



template<typename eT>
arma_inline
const Op<subview_col<eT>, op_cumsum_vec>
cumsum(const subview_col<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<subview_col<eT>, op_cumsum_vec>(A);
  }



template<typename eT>
arma_inline
const Op<diagview<eT>, op_cumsum_vec>
cumsum(const diagview<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<diagview<eT>, op_cumsum_vec>(A);
  }



template<typename eT, typename T1>
arma_inline
const Op<subview_elem1<eT,T1>, op_cumsum_vec>
cumsum(const subview_elem1<eT,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<subview_elem1<eT,T1>, op_cumsum_vec>(A);
  }



//! @}
