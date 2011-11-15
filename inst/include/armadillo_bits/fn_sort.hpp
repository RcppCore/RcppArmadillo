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


//! \addtogroup fn_sort
//! @{


template<typename T1>
arma_inline
const Op<T1, op_sort>
sort(const Base<typename T1::elem_type,T1>& X, const uword sort_type = 0, const uword dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_sort>(X.get_ref(), sort_type, dim);
  }



template<typename eT>
arma_inline
const Op<Col<eT>, op_sort>
sort(const Col<eT>& X, const uword sort_type = 0)
  {
  arma_extra_debug_sigprint();
  
  const uword dim = 0;
  
  return Op<Col<eT>, op_sort>(X, sort_type, dim);
  }



template<typename eT>
arma_inline
const Op<Row<eT>, op_sort>
sort(const Row<eT>& X, const uword sort_type = 0)
  {
  arma_extra_debug_sigprint();
  
  const uword dim = 1;
  
  return Op<Row<eT>, op_sort>(X, sort_type, dim);
  }



//! @}
