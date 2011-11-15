// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup glue_mixed
//! @{



//! matrix multiplication with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_times::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  // TODO: extend the unwrap_check framework to handle mixed matrix types
  
  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B);
  
  const Mat<eT1>& A = tmp1.M;
  const Mat<eT2>& B = tmp2.M;
  
  const bool A_is_alias = ( ((void *)&out) == ((void *)&A) );
  const bool B_is_alias = ( ((void *)&out) == ((void *)&B) );
  
  const Mat<eT1>* AA_ptr = A_is_alias ? new Mat<eT1>(A) : 0;
  const Mat<eT2>* BB_ptr = B_is_alias ? new Mat<eT2>(B) : 0;
  
  const Mat<eT1>& AA = A_is_alias ? *AA_ptr : A;
  const Mat<eT2>& BB = B_is_alias ? *BB_ptr : B;
  
  arma_debug_assert_mul_size(AA, BB, "multiplication");
  
  out.set_size(AA.n_rows, BB.n_cols);
  
  gemm_mixed<>::apply(out, AA, BB);
  
  if(A_is_alias == true)
    {
    delete AA_ptr;
    }
  
  if(B_is_alias == true)
    {
    delete BB_ptr;
    }
  }



//! matrix addition with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_plus::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "addition");
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
        out_eT* out_mem = out.memptr();
  const uword     n_elem  = out.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) + upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//! matrix subtraction with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_minus::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_minus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "subtraction");
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
        out_eT* out_mem = out.memptr();
  const uword     n_elem  = out.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) - upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//! element-wise matrix division with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_div::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_div>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "element-wise division");
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
        out_eT* out_mem = out.memptr();
  const uword     n_elem  = out.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) / upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//! element-wise matrix multiplication with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_schur::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "element-wise multiplication");
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
        out_eT* out_mem = out.memptr();
  const uword     n_elem  = out.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) * upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//
//
//



//! cube addition with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_plus::apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.A);
  const ProxyCube<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "addition");
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
        out_eT* out_mem = out.memptr();
  const uword     n_elem  = out.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) + upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//! cube subtraction with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_minus::apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_minus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.A);
  const ProxyCube<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "subtraction");
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
        out_eT* out_mem = out.memptr();
  const uword     n_elem  = out.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) - upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//! element-wise cube division with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_div::apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_div>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.A);
  const ProxyCube<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "element-wise division");
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
        out_eT* out_mem = out.memptr();
  const uword     n_elem  = out.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) / upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//! element-wise cube multiplication with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_schur::apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.A);
  const ProxyCube<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "element-wise multiplication");
  
  out.set_size(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  
        out_eT* out_mem = out.memptr();
  const uword     n_elem  = out.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) * upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//! @}
