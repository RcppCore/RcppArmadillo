// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_cov
//! @{



template<typename eT>
inline
void
op_cov::direct_cov(Mat<eT>& out, const Mat<eT>& A, const uword norm_type)
  {
  arma_extra_debug_sigprint();
  
  if(A.is_vec())
    {
    if(A.n_rows == 1)
      {
      out = var(trans(A), norm_type);      
      }
    else
      {
      out = var(A, norm_type);
      }
    }
  else
    {
    const uword N = A.n_rows;
    const eT norm_val = (norm_type == 0) ? ( (N > 1) ? eT(N-1) : eT(1) ) : eT(N);

    const Row<eT> acc = sum(A);

    out = trans(A) * A;
    out -= (trans(acc) * acc)/eT(N);
    out /= norm_val;
    }
  }



template<typename T>
inline
void
op_cov::direct_cov(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A, const uword norm_type)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  if(A.is_vec())
    {
    if(A.n_rows == 1)
      {
      const Mat<T> tmp_mat = var(trans(A), norm_type);
      out.set_size(1,1);
      out[0] = tmp_mat[0];
      }
    else
      {
      const Mat<T> tmp_mat = var(A, norm_type);
      out.set_size(1,1);
      out[0] = tmp_mat[0];
      }
    }
  else
    {
    const uword N = A.n_rows;
    const eT norm_val = (norm_type == 0) ? ( (N > 1) ? eT(N-1) : eT(1) ) : eT(N);
    
    const Row<eT> acc = sum(A);
    
    out = trans(A) * A;               // out = strans(conj(A)) * A;
    out -= (trans(acc) * acc)/eT(N);  // out -= (strans(conj(acc)) * acc)/eT(N);
    out /= norm_val;
    }
  }



template<typename T1>
inline
void
op_cov::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cov>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& A     = tmp.M;
  
  const uword norm_type = in.aux_uword_a;
  
  op_cov::direct_cov(out, A, norm_type);
  }



//! @}
