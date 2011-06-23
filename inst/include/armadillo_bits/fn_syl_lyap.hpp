// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_syl_lyap
//! @{


//! find the solution of the Sylvester equation AX + XB = C
template<typename T1, typename T2, typename T3>
inline
bool
syl
  (
        Mat <typename T1::elem_type>   & out,
  const Base<typename T1::elem_type,T1>& in_A,
  const Base<typename T1::elem_type,T2>& in_B,
  const Base<typename T1::elem_type,T3>& in_C,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp_A(in_A.get_ref(), out);
  const unwrap_check<T2> tmp_B(in_B.get_ref(), out);
  const unwrap_check<T3> tmp_C(in_C.get_ref(), out);
  
  const Mat<eT>& A = tmp_A.M;
  const Mat<eT>& B = tmp_B.M;
  const Mat<eT>& C = tmp_C.M;
  
  const bool status = auxlib::syl(out, A, B, C);
  
  if(status == false)
    {
    out.reset();
    arma_bad("syl(): equation appears to be singular", false);
    }
  
  return false;
  }



template<typename T1, typename T2, typename T3>
inline
Mat<typename T1::elem_type>
syl
  (
  const Base<typename T1::elem_type,T1>& in_A,
  const Base<typename T1::elem_type,T2>& in_B,
  const Base<typename T1::elem_type,T3>& in_C,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> out;
  const bool status = syl(out, in_A, in_B, in_C);
  
  if(status == false)
    {
    out.reset();
    arma_bad("syl(): equation appears to be singular");
    }
  
  return out;
  }



//! @}
