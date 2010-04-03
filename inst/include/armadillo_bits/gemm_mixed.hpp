// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup gemm_mixed
//! @{



//! \brief
//! Matrix multplication where the matrices have different element types.
//! Uses caching for speedup.
//! Matrix 'C' is assumed to have been set to the correct size (i.e. taking into account transposes)

template<const bool do_trans_A=false, const bool do_trans_B=false, const bool use_alpha=false, const bool use_beta=false>
class gemm_mixed_cache
  {
  public:
  
  template<typename out_eT, typename in_eT1, typename in_eT2>
  arma_hot
  inline
  static
  void
  apply
    (
          Mat<out_eT>& C,
    const Mat<in_eT1>& A,
    const Mat<in_eT2>& B,
    const out_eT alpha = out_eT(1),
    const out_eT beta  = out_eT(0)
    )
    {
    arma_extra_debug_sigprint();
    
    const u32 A_n_rows = A.n_rows;
    const u32 A_n_cols = A.n_cols;
    
    const u32 B_n_rows = B.n_rows;
    const u32 B_n_cols = B.n_cols;
    
    if( (do_trans_A == false) && (do_trans_B == false) )
      {
      podarray<in_eT1> tmp(A_n_cols);
      in_eT1* A_rowdata = tmp.memptr();
      
      for(u32 row_A=0; row_A < A_n_rows; ++row_A)
        {
        
        for(u32 col_A=0; col_A < A_n_cols; ++col_A)
          {
          A_rowdata[col_A] = A.at(row_A,col_A);
          }
        
        for(u32 col_B=0; col_B < B_n_cols; ++col_B)
          {
          const in_eT2* B_coldata = B.colptr(col_B);
          
          out_eT acc = out_eT(0);
          for(u32 i=0; i < B_n_rows; ++i)
            {
            acc += upgrade_val<in_eT1,in_eT2>::apply(A_rowdata[i]) * upgrade_val<in_eT1,in_eT2>::apply(B_coldata[i]);
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(row_A,col_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(row_A,col_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(row_A,col_B) = acc + beta*C.at(row_A,col_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(row_A,col_B) = alpha*acc + beta*C.at(row_A,col_B);
            }
          
          }
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == false) )
      {
      for(u32 col_A=0; col_A < A_n_cols; ++col_A)
        {
        // col_A is interpreted as row_A when storing the results in matrix C
        
        const in_eT1* A_coldata = A.colptr(col_A);
        
        for(u32 col_B=0; col_B < B_n_cols; ++col_B)
          {
          const in_eT2* B_coldata = B.colptr(col_B);
          
          out_eT acc = out_eT(0);
          for(u32 i=0; i < B_n_rows; ++i)
            {
            acc += upgrade_val<in_eT1,in_eT2>::apply(A_coldata[i]) * upgrade_val<in_eT1,in_eT2>::apply(B_coldata[i]);
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(col_A,col_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(col_A,col_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(col_A,col_B) = acc + beta*C.at(col_A,col_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(col_A,col_B) = alpha*acc + beta*C.at(col_A,col_B);
            }
          
          }
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == true) )
      {
      Mat<in_eT2> B_tmp = trans(B);
      gemm_mixed_cache<false, false, use_alpha, use_beta>::apply(C, A, B_tmp, alpha, beta);
      }
    else
    if( (do_trans_A == true) && (do_trans_B == true) )
      {
      // mat B_tmp = trans(B);
      // dgemm_arma<true, false,  use_alpha, use_beta>::apply(C, A, B_tmp, alpha, beta);
      
      
      // By using the trans(A)*trans(B) = trans(B*A) equivalency,
      // transpose operations are not needed
      
      podarray<in_eT2> tmp(B.n_cols);
      in_eT2* B_rowdata = tmp.memptr();
      
      for(u32 row_B=0; row_B < B_n_rows; ++row_B)
        {
        
        for(u32 col_B=0; col_B < B_n_cols; ++col_B)
          {
          B_rowdata[col_B] = B.at(row_B,col_B);
          }
        
        for(u32 col_A=0; col_A < A_n_cols; ++col_A)
          {
          const in_eT1* A_coldata = A.colptr(col_A);
          
          out_eT acc = out_eT(0);
          for(u32 i=0; i < A_n_rows; ++i)
            {
            acc += upgrade_val<in_eT1,in_eT2>::apply(B_rowdata[i]) * upgrade_val<in_eT1,in_eT2>::apply(A_coldata[i]);
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(col_A,row_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(col_A,row_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(col_A,row_B) = acc + beta*C.at(col_A,row_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(col_A,row_B) = alpha*acc + beta*C.at(col_A,row_B);
            }
          
          }
        }
      
      }
    }
    
  };



//! Matrix multplication where the matrices have different element types.
//! Simple version (no caching).
//! Matrix 'C' is assumed to have been set to the correct size (i.e. taking into account transposes)
template<const bool do_trans_A=false, const bool do_trans_B=false, const bool use_alpha=false, const bool use_beta=false>
class gemm_mixed_simple
  {
  public:
  
  template<typename out_eT, typename in_eT1, typename in_eT2>
  arma_hot
  inline
  static
  void
  apply
    (
          Mat<out_eT>& C,
    const Mat<in_eT1>& A,
    const Mat<in_eT2>& B,
    const out_eT alpha = out_eT(1),
    const out_eT beta  = out_eT(0)
    )
    {
    arma_extra_debug_sigprint();
    
    const u32 A_n_rows = A.n_rows;
    const u32 A_n_cols = A.n_cols;
    
    const u32 B_n_rows = B.n_rows;
    const u32 B_n_cols = B.n_cols;
    
    if( (do_trans_A == false) && (do_trans_B == false) )
      {
      for(u32 row_A = 0; row_A < A_n_rows; ++row_A)
        {
        for(u32 col_B = 0; col_B < B_n_cols; ++col_B)
          {
          const in_eT2* B_coldata = B.colptr(col_B);
          
          out_eT acc = out_eT(0);
          for(u32 i = 0; i < B_n_rows; ++i)
            {
            const out_eT val1 = upgrade_val<in_eT1,in_eT2>::apply(A.at(row_A,i));
            const out_eT val2 = upgrade_val<in_eT1,in_eT2>::apply(B_coldata[i]);
            acc += val1 * val2;
            //acc += upgrade_val<in_eT1,in_eT2>::apply(A.at(row_A,i)) * upgrade_val<in_eT1,in_eT2>::apply(B_coldata[i]);
            }
          
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(row_A,col_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(row_A,col_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(row_A,col_B) = acc + beta*C.at(row_A,col_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(row_A,col_B) = alpha*acc + beta*C.at(row_A,col_B);
            }
          }
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == false) )
      {
      for(u32 col_A=0; col_A < A_n_cols; ++col_A)
        {
        // col_A is interpreted as row_A when storing the results in matrix C
        
        const in_eT1* A_coldata = A.colptr(col_A);
        
        for(u32 col_B=0; col_B < B_n_cols; ++col_B)
          {
          const in_eT2* B_coldata = B.colptr(col_B);
          
          out_eT acc = out_eT(0);
          for(u32 i=0; i < B_n_rows; ++i)
            {
            acc += upgrade_val<in_eT1,in_eT2>::apply(A_coldata[i]) * upgrade_val<in_eT1,in_eT2>::apply(B_coldata[i]);
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(col_A,col_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(col_A,col_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(col_A,col_B) = acc + beta*C.at(col_A,col_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(col_A,col_B) = alpha*acc + beta*C.at(col_A,col_B);
            }
          
          }
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == true) )
      {
      for(u32 row_A = 0; row_A < A_n_rows; ++row_A)
        {
        for(u32 row_B = 0; row_B < B_n_rows; ++row_B)
          {
          out_eT acc = out_eT(0);
          for(u32 i = 0; i < B_n_cols; ++i)
            {
            acc += upgrade_val<in_eT1,in_eT2>::apply(A.at(row_A,i)) * upgrade_val<in_eT1,in_eT2>::apply(B.at(row_B,i));
            }
          
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(row_A,row_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(row_A,row_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(row_A,row_B) = acc + beta*C.at(row_A,row_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(row_A,row_B) = alpha*acc + beta*C.at(row_A,row_B);
            }
          }
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == true) )
      {
      for(u32 row_B=0; row_B < B_n_rows; ++row_B)
        {
        
        for(u32 col_A=0; col_A < A_n_cols; ++col_A)
          {
          const in_eT1* A_coldata = A.colptr(col_A);
          
          out_eT acc = out_eT(0);
          for(u32 i=0; i < A_n_rows; ++i)
            {
            acc += upgrade_val<in_eT1,in_eT2>::apply(B.at(row_B,i)) * upgrade_val<in_eT1,in_eT2>::apply(A_coldata[i]);
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(col_A,row_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(col_A,row_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(col_A,row_B) = acc + beta*C.at(col_A,row_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(col_A,row_B) = alpha*acc + beta*C.at(col_A,row_B);
            }
          
          }
        }
      
      }
    }
    
  };





//! \brief
//! Matrix multplication where the matrices have different element types.

template<const bool do_trans_A=false, const bool do_trans_B=false, const bool use_alpha=false, const bool use_beta=false>
class gemm_mixed
  {
  public:
  
  //! immediate multiplication of matrices A and B, storing the result in C
  template<typename out_eT, typename in_eT1, typename in_eT2>
  inline
  static
  void
  apply
    (
          Mat<out_eT>& C,
    const Mat<in_eT1>& A,
    const Mat<in_eT2>& B,
    const out_eT alpha = out_eT(1),
    const out_eT beta  = out_eT(0)
    )
    {
    arma_extra_debug_sigprint();
    
    if( (A.n_elem <= 64u) && (B.n_elem <= 64u) )
      {
      gemm_mixed_simple<do_trans_A, do_trans_B, use_alpha, use_beta>::apply(C,A,B,alpha,beta);
      }
    else
      {
      gemm_mixed_cache<do_trans_A, do_trans_B, use_alpha, use_beta>::apply(C,A,B,alpha,beta);
      }
    }
  
  };



//! @}
