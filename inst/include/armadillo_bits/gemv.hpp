// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup gemv
//! @{



//! for small square matrices with n_rows <= 4
template<const bool do_trans_A=false, const bool use_alpha=false, const bool use_beta=false>
class gemv_emul_tiny
  {
  public:
  
  
  template<const bool do_flip, const u32 row, const u32 col>
  struct pos
    {
    static const u32 n2 = (do_flip == false) ? (row + col*2) : (col + row*2);
    static const u32 n3 = (do_flip == false) ? (row + col*3) : (col + row*3);
    static const u32 n4 = (do_flip == false) ? (row + col*4) : (col + row*4);
    };
  
  
  
  template<typename eT, const bool do_alpha_mul, const bool do_beta_mul, const u32 i>
  arma_hot
  arma_inline
  static
  void
  assign(eT* y, const eT acc, const eT alpha, const eT beta)
    {
    if(do_beta_mul == false)
      {
      y[i] = (do_alpha_mul == false) ? acc : alpha*acc;
      }
    else
      {
      const eT tmp = y[i];
      
      y[i] = beta*tmp + ( (do_alpha_mul == false) ? acc : alpha*acc );
      }
    }
  
  

  template<typename eT>
  arma_hot
  inline
  static
  void
  apply( eT* y, const Mat<eT>& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    arma_extra_debug_sigprint();
    
    const eT*  Am = A.memptr();
    const bool tA = do_trans_A;
    
    switch(A.n_rows)
      {
      case 1:
        {
        y[0] = Am[0] * x[0];
        }
        break;
      
      
      case 2:
        {
        const eT x0 = x[0];
        const eT x1 = x[1];
        
        const eT acc0 = Am[pos<tA,0,0>::n2]*x0 + Am[pos<tA,0,1>::n2]*x1;
        const eT acc1 = Am[pos<tA,1,0>::n2]*x0 + Am[pos<tA,1,1>::n2]*x1;
        
        assign<eT, use_alpha, use_beta, 0>(y, acc0, alpha, beta);
        assign<eT, use_alpha, use_beta, 1>(y, acc1, alpha, beta);
        }
        break;
      
        
      case 3:
        {
        const eT x0 = x[0];
        const eT x1 = x[1];
        const eT x2 = x[2];
        
        const eT acc0 = Am[pos<tA,0,0>::n3]*x0 + Am[pos<tA,0,1>::n3]*x1 + Am[pos<tA,0,2>::n3]*x2;
        const eT acc1 = Am[pos<tA,1,0>::n3]*x0 + Am[pos<tA,1,1>::n3]*x1 + Am[pos<tA,1,2>::n3]*x2;
        const eT acc2 = Am[pos<tA,2,0>::n3]*x0 + Am[pos<tA,2,1>::n3]*x1 + Am[pos<tA,2,2>::n3]*x2;
        
        assign<eT, use_alpha, use_beta, 0>(y, acc0, alpha, beta);
        assign<eT, use_alpha, use_beta, 1>(y, acc1, alpha, beta);
        assign<eT, use_alpha, use_beta, 2>(y, acc2, alpha, beta);
        }
        break;
      
      
      case 4:
        {
        const eT x0 = x[0];
        const eT x1 = x[1];
        const eT x2 = x[2];
        const eT x3 = x[3];
        
        const eT acc0 = Am[pos<tA,0,0>::n4]*x0 + Am[pos<tA,0,1>::n4]*x1 + Am[pos<tA,0,2>::n4]*x2 + Am[pos<tA,0,3>::n4]*x3;
        const eT acc1 = Am[pos<tA,1,0>::n4]*x0 + Am[pos<tA,1,1>::n4]*x1 + Am[pos<tA,1,2>::n4]*x2 + Am[pos<tA,1,3>::n4]*x3;
        const eT acc2 = Am[pos<tA,2,0>::n4]*x0 + Am[pos<tA,2,1>::n4]*x1 + Am[pos<tA,2,2>::n4]*x2 + Am[pos<tA,2,3>::n4]*x3;
        const eT acc3 = Am[pos<tA,3,0>::n4]*x0 + Am[pos<tA,3,1>::n4]*x1 + Am[pos<tA,3,2>::n4]*x2 + Am[pos<tA,3,3>::n4]*x3;
        
        assign<eT, use_alpha, use_beta, 0>(y, acc0, alpha, beta);
        assign<eT, use_alpha, use_beta, 1>(y, acc1, alpha, beta);
        assign<eT, use_alpha, use_beta, 2>(y, acc2, alpha, beta);
        assign<eT, use_alpha, use_beta, 3>(y, acc3, alpha, beta);
        }
        break;
      
      
      default:
        ;
      }
    }
    
  };



//! \brief
//! Partial emulation of ATLAS/BLAS gemv().
//! 'y' is assumed to have been set to the correct size (i.e. taking into account the transpose)

template<const bool do_trans_A=false, const bool use_alpha=false, const bool use_beta=false>
class gemv_emul_large
  {
  public:
  
  template<typename eT>
  arma_hot
  inline
  static
  void
  apply( eT* y, const Mat<eT>& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    arma_extra_debug_sigprint();
    
    const u32 A_n_rows = A.n_rows;
    const u32 A_n_cols = A.n_cols;
    
    if(do_trans_A == false)
      {
      for(u32 row=0; row < A_n_rows; ++row)
        {
        eT acc = eT(0);
        
        for(u32 i=0; i < A_n_cols; ++i)
          {
          acc += A.at(row,i) * x[i];
          }
        
        if( (use_alpha == false) && (use_beta == false) )
          {
          y[row] = acc;
          }
        else
        if( (use_alpha == true) && (use_beta == false) )
          {
          y[row] = alpha * acc;
          }
        else
        if( (use_alpha == false) && (use_beta == true) )
          {
          y[row] = acc + beta*y[row];
          }
        else
        if( (use_alpha == true) && (use_beta == true) )
          {
          y[row] = alpha*acc + beta*y[row];
          }
        }
      }
    else
    if(do_trans_A == true)
      {
      for(u32 col=0; col < A_n_cols; ++col)
        {
        // col is interpreted as row when storing the results in 'y'
        
        
        // const eT* A_coldata = A.colptr(col);
        // 
        // eT acc = eT(0);
        // for(u32 row=0; row < A_n_rows; ++row)
        //   {
        //   acc += A_coldata[row] * x[row];
        //   }
        
        const eT acc = op_dot::direct_dot_arma(A_n_rows, A.colptr(col), x);
        
        if( (use_alpha == false) && (use_beta == false) )
          {
          y[col] = acc;
          }
        else
        if( (use_alpha == true) && (use_beta == false) )
          {
          y[col] = alpha * acc;
          }
        else
        if( (use_alpha == false) && (use_beta == true) )
          {
          y[col] = acc + beta*y[col];
          }
        else
        if( (use_alpha == true) && (use_beta == true) )
          {
          y[col] = alpha*acc + beta*y[col];
          }
        
        }
      }
    }
  
  };



template<const bool do_trans_A=false, const bool use_alpha=false, const bool use_beta=false>
class gemv_emul
  {
  public:
  
  template<typename eT>
  arma_hot
  inline
  static
  void
  apply( eT* y, const Mat<eT>& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0), const typename arma_not_cx<eT>::result* junk = 0 )
    {
    arma_extra_debug_sigprint();
    arma_ignore(junk);
    
    if( A.n_rows > 4 )
      {
      gemv_emul_large<do_trans_A, use_alpha, use_beta>::apply(y, A, x, alpha, beta);
      }
    else
      {
      gemv_emul_tiny<do_trans_A, use_alpha, use_beta>::apply(y, A, x, alpha, beta);
      }
    }
  
  
  
  template<typename eT>
  arma_hot
  inline
  static
  void
  apply( eT* y, const Mat<eT>& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0), const typename arma_cx_only<eT>::result* junk = 0 )
    {
    arma_extra_debug_sigprint();
    
    Mat<eT> tmp_A;
    
    if(do_trans_A)
      {
      op_htrans::apply_noalias(tmp_A, A);
      }
    
    const Mat<eT>& AA = (do_trans_A == false) ? A : tmp_A;
    
    if( AA.n_rows > 4 )
      {
      gemv_emul_large<false, use_alpha, use_beta>::apply(y, AA, x, alpha, beta);
      }
    else
      {
      gemv_emul_tiny<false, use_alpha, use_beta>::apply(y, AA, x, alpha, beta);
      }
    }
  };



//! \brief
//! Wrapper for ATLAS/BLAS gemv function, using template arguments to control the arguments passed to gemv.
//! 'y' is assumed to have been set to the correct size (i.e. taking into account the transpose)

template<const bool do_trans_A=false, const bool use_alpha=false, const bool use_beta=false>
class gemv
  {
  public:
  
  template<typename eT>
  inline
  static
  void
  apply_blas_type( eT* y, const Mat<eT>& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    arma_extra_debug_sigprint();
    
    if(A.n_elem <= 64u)
      {
      gemv_emul<do_trans_A, use_alpha, use_beta>::apply(y,A,x,alpha,beta);
      }
    else
      {
      #if defined(ARMA_USE_ATLAS)
        {
        arma_extra_debug_print("atlas::cblas_gemv()");
        
        atlas::cblas_gemv<eT>
          (
          atlas::CblasColMajor,
          (do_trans_A) ? ( is_complex<eT>::value ? CblasConjTrans : atlas::CblasTrans ) : atlas::CblasNoTrans,
          A.n_rows,
          A.n_cols,
          (use_alpha) ? alpha : eT(1),
          A.mem,
          A.n_rows,
          x,
          1,
          (use_beta) ? beta : eT(0),
          y,
          1
          );
        }
      #elif defined(ARMA_USE_BLAS)
        {
        arma_extra_debug_print("blas::gemv()");
        
        const char      trans_A     = (do_trans_A) ? ( is_complex<eT>::value ? 'C' : 'T' ) : 'N';
        const blas_int  m           = A.n_rows;
        const blas_int  n           = A.n_cols;
        const eT        local_alpha = (use_alpha) ? alpha : eT(1);
        //const blas_int  lda         = A.n_rows;
        const blas_int  inc         = 1;
        const eT        local_beta  = (use_beta) ? beta : eT(0);
        
        arma_extra_debug_print( arma_boost::format("blas::gemv(): trans_A = %c") % trans_A );
        
        blas::gemv<eT>
          (
          &trans_A,
          &m,
          &n,
          &local_alpha,
          A.mem,
          &m,  // lda
          x,
          &inc,
          &local_beta,
          y,
          &inc
          );
        }
      #else
        {
        gemv_emul<do_trans_A, use_alpha, use_beta>::apply(y,A,x,alpha,beta);
        }
      #endif
      }
    
    }
  
  
  
  template<typename eT>
  arma_inline
  static
  void
  apply( eT* y, const Mat<eT>& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    gemv_emul<do_trans_A, use_alpha, use_beta>::apply(y,A,x,alpha,beta);
    }
  
  
  
  arma_inline
  static
  void
  apply
    (
          float*      y,
    const Mat<float>& A,
    const float*      x,
    const float       alpha = float(1),
    const float       beta  = float(0)
    )
    {
    gemv<do_trans_A, use_alpha, use_beta>::apply_blas_type(y,A,x,alpha,beta);
    }


  
  arma_inline
  static
  void
  apply
    (
          double*      y,
    const Mat<double>& A,
    const double*      x,
    const double       alpha = double(1),
    const double       beta  = double(0)
    )
    {
    gemv<do_trans_A, use_alpha, use_beta>::apply_blas_type(y,A,x,alpha,beta);
    }


  
  arma_inline
  static
  void
  apply
    (
          std::complex<float>*         y,
    const Mat< std::complex<float > >& A,
    const std::complex<float>*         x,
    const std::complex<float>          alpha = std::complex<float>(1),
    const std::complex<float>          beta  = std::complex<float>(0)
    )
    {
    gemv<do_trans_A, use_alpha, use_beta>::apply_blas_type(y,A,x,alpha,beta);
    }


  
  arma_inline
  static
  void
  apply
    (
          std::complex<double>*        y,
    const Mat< std::complex<double> >& A,
    const std::complex<double>*        x,
    const std::complex<double>         alpha = std::complex<double>(1),
    const std::complex<double>         beta  = std::complex<double>(0)
    )
    {
    gemv<do_trans_A, use_alpha, use_beta>::apply_blas_type(y,A,x,alpha,beta);
    }


  
  };


//! @}
