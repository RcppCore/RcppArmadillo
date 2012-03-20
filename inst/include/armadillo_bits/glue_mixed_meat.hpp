// Copyright (C) 2009-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2012 Conrad Sanderson
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
  
  const uword n_rows = A.get_n_rows();
  const uword n_cols = A.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
        out_eT* out_mem = out.memptr();
  const uword   n_elem  = out.n_elem;
    
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type AA = A.get_ea();
    typename Proxy<T2>::ea_type BB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) + upgrade_val<eT1,eT2>::apply(BB[i]);
      }
    }
  else
    {
    uword i = 0;
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(A.at(row,col)) + upgrade_val<eT1,eT2>::apply(B.at(row,col));
      ++i;
      }
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
  
  const uword n_rows = A.get_n_rows();
  const uword n_cols = A.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
        out_eT* out_mem = out.memptr();
  const uword   n_elem  = out.n_elem;
    
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type AA = A.get_ea();
    typename Proxy<T2>::ea_type BB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) - upgrade_val<eT1,eT2>::apply(BB[i]);
      }
    }
  else
    {
    uword i = 0;
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(A.at(row,col)) - upgrade_val<eT1,eT2>::apply(B.at(row,col));
      ++i;
      }
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
  
  const uword n_rows = A.get_n_rows();
  const uword n_cols = A.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
        out_eT* out_mem = out.memptr();
  const uword   n_elem  = out.n_elem;
    
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type AA = A.get_ea();
    typename Proxy<T2>::ea_type BB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) / upgrade_val<eT1,eT2>::apply(BB[i]);
      }
    }
  else
    {
    uword i = 0;
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(A.at(row,col)) / upgrade_val<eT1,eT2>::apply(B.at(row,col));
      ++i;
      }
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
  
  const uword n_rows = A.get_n_rows();
  const uword n_cols = A.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
        out_eT* out_mem = out.memptr();
  const uword   n_elem  = out.n_elem;
    
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type AA = A.get_ea();
    typename Proxy<T2>::ea_type BB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) * upgrade_val<eT1,eT2>::apply(BB[i]);
      }
    }
  else
    {
    uword i = 0;
    
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(A.at(row,col)) * upgrade_val<eT1,eT2>::apply(B.at(row,col));
      ++i;
      }
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
  
  // TODO: add faster handling of subviews
  
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
  
  // TODO: add faster handling of subviews
  
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
  
  // TODO: add faster handling of subviews
  
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
  
  // TODO: add faster handling of subviews
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) * upgrade_val<eT1,eT2>::apply(B[i]);
    }
  }



//! @}
