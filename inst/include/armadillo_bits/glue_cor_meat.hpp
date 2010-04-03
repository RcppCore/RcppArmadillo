// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup glue_cor
//! @{



template<typename eT>
inline
void
glue_cor::direct_cor(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const u32 norm_type)
  {
  arma_extra_debug_sigprint();

  if(A.is_vec() && B.is_vec())
    {
    arma_debug_check( (A.n_elem != B.n_elem), "cor(): the number of elements in A and B must match" );

    const eT* A_ptr = A.memptr();
    const eT* B_ptr = B.memptr();

    eT A_acc   = eT(0);
    eT B_acc   = eT(0);
    eT out_acc = eT(0);

    const u32 N = A.n_elem;

    for(u32 i=0; i<N; ++i)
      {
      const eT A_tmp = A_ptr[i];
      const eT B_tmp = B_ptr[i];

      A_acc += A_tmp;
      B_acc += B_tmp;

      out_acc += A_tmp * B_tmp;
      }

    out_acc -= (A_acc * B_acc)/eT(N);

    const eT norm_val = (norm_type == 0) ? ( (N > 1) ? eT(N-1) : eT(1) ) : eT(N);    

    out.set_size(1,1);
    out[0] = out_acc/norm_val;

    const Mat<eT> stddev_A = (A.n_rows == 1) ? stddev(trans(A)) : stddev(A);
    const Mat<eT> stddev_B = (B.n_rows == 1) ? stddev(trans(B)) : stddev(B);

    out /= stddev_A * stddev_B;
    }
  else
    {
    arma_debug_assert_same_size(A, B, "cor()");

    const u32 N = A.n_rows;
    const eT norm_val = (norm_type == 0) ? ( (N > 1) ? eT(N-1) : eT(1) ) : eT(N);

    out = trans(A) * B;
    out -= (trans(sum(A)) * sum(B))/eT(N);
    out /= norm_val;
    out /= trans(stddev(A)) * stddev(B);
    }
  }



template<typename T>
inline
void
glue_cor::direct_cor(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A, const Mat< std::complex<T> >& B, const u32 norm_type)
  {
  arma_extra_debug_sigprint();

  typedef typename std::complex<T> eT;

  if(A.is_vec() && B.is_vec())
    { 
    arma_debug_check( (A.n_elem != B.n_elem), "cor(): the number of elements in A and B must match" );

    const eT* A_ptr = A.memptr();
    const eT* B_ptr = B.memptr();        

    eT A_acc   = eT(0);
    eT B_acc   = eT(0);
    eT out_acc = eT(0);

    const u32 N = A.n_elem;

    for(u32 i=0; i<N; ++i)
      {
      const eT A_tmp = A_ptr[i];
      const eT B_tmp = B_ptr[i];

      A_acc += A_tmp;
      B_acc += B_tmp;

      out_acc += std::conj(A_tmp) * B_tmp;
      }

    out_acc -= (std::conj(A_acc) * B_acc)/eT(N);

    const eT norm_val = (norm_type == 0) ? ( (N > 1) ? eT(N-1) : eT(1) ) : eT(N);

    out.set_size(1,1);
    out[0] = out_acc/norm_val;

    const Mat<T> stddev_A = (A.n_rows == 1) ? stddev(trans(A)) : stddev(A);
    const Mat<T> stddev_B = (B.n_rows == 1) ? stddev(trans(B)) : stddev(B);

    out /= conv_to< Mat<eT> >::from( stddev_A * stddev_B );
    }
  else
    {
    arma_debug_assert_same_size(A, B, "cor()");

    const u32 N = A.n_rows;
    const eT norm_val = (norm_type == 0) ? ( (N > 1) ? eT(N-1) : eT(1) ) : eT(N);

    out = trans(conj(A)) * B;
    out -= (trans(conj(sum(A))) * sum(B))/eT(N);
    out /= norm_val;
    out /= conv_to< Mat<eT> >::from( trans(stddev(A)) * stddev(B) );
    }
  }



template<typename T1, typename T2>
inline
void
glue_cor::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_cor>& X)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const unwrap_check<T1> A_tmp(X.A, out);
  const unwrap_check<T2> B_tmp(X.B, out);

  const Mat<eT>& A = A_tmp.M;
  const Mat<eT>& B = B_tmp.M;
  
  const u32 norm_type = X.aux_u32;
  
  if(&A != &B)
    {
    glue_cor::direct_cor(out, A, B, norm_type);
    }
  else
    {
    op_cor::direct_cor(out, A, norm_type);
    }
  
  }



//! @}
