// Copyright (C) 2008-2013 Conrad Sanderson
// Copyright (C) 2008-2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup gemv
//! @{



//! for tiny square matrices, size <= 4x4
template<const bool do_trans_A=false, const bool use_alpha=false, const bool use_beta=false>
class gemv_emul_tinysq
  {
  public:
  
  
  template<const uword row, const uword col>
  struct pos
    {
    static const uword n2 = (do_trans_A == false) ? (row + col*2) : (col + row*2);
    static const uword n3 = (do_trans_A == false) ? (row + col*3) : (col + row*3);
    static const uword n4 = (do_trans_A == false) ? (row + col*4) : (col + row*4);
    };
  
  
  
  template<typename eT, const uword i>
  arma_hot
  arma_inline
  static
  void
  assign(eT* y, const eT acc, const eT alpha, const eT beta)
    {
    if(use_beta == false)
      {
      y[i] = (use_alpha == false) ? acc : alpha*acc;
      }
    else
      {
      const eT tmp = y[i];
      
      y[i] = beta*tmp + ( (use_alpha == false) ? acc : alpha*acc );
      }
    }
  
  

  template<typename eT, typename TA>
  arma_hot
  inline
  static
  void
  apply( eT* y, const TA& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    arma_extra_debug_sigprint();
    
    const eT*  Am = A.memptr();
    
    switch(A.n_rows)
      {
      case 1:
        {
        const eT acc = Am[0] * x[0];
        
        assign<eT, 0>(y, acc, alpha, beta);
        }
        break;
      
      
      case 2:
        {
        const eT x0 = x[0];
        const eT x1 = x[1];
        
        const eT acc0 = Am[pos<0,0>::n2]*x0 + Am[pos<0,1>::n2]*x1;
        const eT acc1 = Am[pos<1,0>::n2]*x0 + Am[pos<1,1>::n2]*x1;
        
        assign<eT, 0>(y, acc0, alpha, beta);
        assign<eT, 1>(y, acc1, alpha, beta);
        }
        break;
      
        
      case 3:
        {
        const eT x0 = x[0];
        const eT x1 = x[1];
        const eT x2 = x[2];
        
        const eT acc0 = Am[pos<0,0>::n3]*x0 + Am[pos<0,1>::n3]*x1 + Am[pos<0,2>::n3]*x2;
        const eT acc1 = Am[pos<1,0>::n3]*x0 + Am[pos<1,1>::n3]*x1 + Am[pos<1,2>::n3]*x2;
        const eT acc2 = Am[pos<2,0>::n3]*x0 + Am[pos<2,1>::n3]*x1 + Am[pos<2,2>::n3]*x2;
        
        assign<eT, 0>(y, acc0, alpha, beta);
        assign<eT, 1>(y, acc1, alpha, beta);
        assign<eT, 2>(y, acc2, alpha, beta);
        }
        break;
      
      
      case 4:
        {
        const eT x0 = x[0];
        const eT x1 = x[1];
        const eT x2 = x[2];
        const eT x3 = x[3];
        
        const eT acc0 = Am[pos<0,0>::n4]*x0 + Am[pos<0,1>::n4]*x1 + Am[pos<0,2>::n4]*x2 + Am[pos<0,3>::n4]*x3;
        const eT acc1 = Am[pos<1,0>::n4]*x0 + Am[pos<1,1>::n4]*x1 + Am[pos<1,2>::n4]*x2 + Am[pos<1,3>::n4]*x3;
        const eT acc2 = Am[pos<2,0>::n4]*x0 + Am[pos<2,1>::n4]*x1 + Am[pos<2,2>::n4]*x2 + Am[pos<2,3>::n4]*x3;
        const eT acc3 = Am[pos<3,0>::n4]*x0 + Am[pos<3,1>::n4]*x1 + Am[pos<3,2>::n4]*x2 + Am[pos<3,3>::n4]*x3;
        
        assign<eT, 0>(y, acc0, alpha, beta);
        assign<eT, 1>(y, acc1, alpha, beta);
        assign<eT, 2>(y, acc2, alpha, beta);
        assign<eT, 3>(y, acc3, alpha, beta);
        }
        break;
      
      
      default:
        ;
      }
    }
    
  };



class gemv_emul_large_helper
  {
  public:
  
  template<typename eT, typename TA>
  arma_hot
  inline
  static
  typename arma_not_cx<eT>::result
  dot_row_col( const TA& A, const eT* x, const uword row, const uword N )
    {
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    
    uword i,j;
    for(i=0, j=1; j < N; i+=2, j+=2)
      {
      const eT xi = x[i];
      const eT xj = x[j];
      
      acc1 += A.at(row,i) * xi;
      acc2 += A.at(row,j) * xj;
      }
    
    if(i < N)
      {
      acc1 += A.at(row,i) * x[i];
      }
    
    return (acc1 + acc2);
    }
  
  
  
  template<typename eT, typename TA>
  arma_hot
  inline
  static
  typename arma_cx_only<eT>::result
  dot_row_col( const TA& A, const eT* x, const uword row, const uword N )
    {
    typedef typename get_pod_type<eT>::result T;
    
    T val_real = T(0);
    T val_imag = T(0);
    
    for(uword i=0; i<N; ++i)
      {
      const std::complex<T>& Ai = A.at(row,i);
      const std::complex<T>& xi = x[i];
      
      const T a = Ai.real();
      const T b = Ai.imag();
      
      const T c = xi.real();
      const T d = xi.imag();
      
      val_real += (a*c) - (b*d);
      val_imag += (a*d) + (b*c);
      }
    
    return std::complex<T>(val_real, val_imag);
    }
  
  };



//! \brief
//! Partial emulation of ATLAS/BLAS gemv().
//! 'y' is assumed to have been set to the correct size (i.e. taking into account the transpose)

template<const bool do_trans_A=false, const bool use_alpha=false, const bool use_beta=false>
class gemv_emul_large
  {
  public:
  
  template<typename eT, typename TA>
  arma_hot
  inline
  static
  void
  apply( eT* y, const TA& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    arma_extra_debug_sigprint();
    
    const uword A_n_rows = A.n_rows;
    const uword A_n_cols = A.n_cols;
    
    if(do_trans_A == false)
      {
      if(A_n_rows == 1)
        {
        const eT acc = op_dot::direct_dot_arma(A_n_cols, A.memptr(), x);
        
             if( (use_alpha == false) && (use_beta == false) )  { y[0] =       acc;             }
        else if( (use_alpha == true ) && (use_beta == false) )  { y[0] = alpha*acc;             }
        else if( (use_alpha == false) && (use_beta == true ) )  { y[0] =       acc + beta*y[0]; }
        else if( (use_alpha == true ) && (use_beta == true ) )  { y[0] = alpha*acc + beta*y[0]; }
        }
      else
      for(uword row=0; row < A_n_rows; ++row)
        {
        const eT acc = gemv_emul_large_helper::dot_row_col(A, x, row, A_n_cols);
        
             if( (use_alpha == false) && (use_beta == false) )  { y[row] =       acc;               }
        else if( (use_alpha == true ) && (use_beta == false) )  { y[row] = alpha*acc;               }
        else if( (use_alpha == false) && (use_beta == true ) )  { y[row] =       acc + beta*y[row]; }
        else if( (use_alpha == true ) && (use_beta == true ) )  { y[row] = alpha*acc + beta*y[row]; }
        }
      }
    else
    if(do_trans_A == true)
      {
      for(uword col=0; col < A_n_cols; ++col)
        {
        // col is interpreted as row when storing the results in 'y'
        
        
        // const eT* A_coldata = A.colptr(col);
        // 
        // eT acc = eT(0);
        // for(uword row=0; row < A_n_rows; ++row)
        //   {
        //   acc += A_coldata[row] * x[row];
        //   }
        
        const eT acc = op_dot::direct_dot_arma(A_n_rows, A.colptr(col), x);
        
             if( (use_alpha == false) && (use_beta == false) )  { y[col] =       acc;               }
        else if( (use_alpha == true ) && (use_beta == false) )  { y[col] = alpha*acc;               }
        else if( (use_alpha == false) && (use_beta == true ) )  { y[col] =       acc + beta*y[col]; }
        else if( (use_alpha == true ) && (use_beta == true ) )  { y[col] = alpha*acc + beta*y[col]; }
        
        }
      }
    }
  
  };



template<const bool do_trans_A=false, const bool use_alpha=false, const bool use_beta=false>
class gemv_emul
  {
  public:
  
  template<typename eT, typename TA>
  arma_hot
  inline
  static
  void
  apply( eT* y, const TA& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0), const typename arma_not_cx<eT>::result* junk = 0 )
    {
    arma_extra_debug_sigprint();
    arma_ignore(junk);
    
    const uword A_n_rows = A.n_rows;
    const uword A_n_cols = A.n_cols;
    
    if( (A_n_rows <= 4) && (A_n_rows == A_n_cols) )
      {
      gemv_emul_tinysq<do_trans_A, use_alpha, use_beta>::apply(y, A, x, alpha, beta);
      }
    else
      {
      gemv_emul_large<do_trans_A, use_alpha, use_beta>::apply(y, A, x, alpha, beta);
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
    arma_ignore(junk);
    
    Mat<eT> tmp_A;
    
    if(do_trans_A)
      {
      op_htrans::apply_noalias(tmp_A, A);
      }
    
    const Mat<eT>& AA = (do_trans_A == false) ? A : tmp_A;
    
    const uword AA_n_rows = AA.n_rows;
    const uword AA_n_cols = AA.n_cols;
    
    if( (AA_n_rows <= 4) && (AA_n_rows == AA_n_cols) )
      {
      gemv_emul_tinysq<false, use_alpha, use_beta>::apply(y, AA, x, alpha, beta);
      }
    else
      {
      gemv_emul_large<false, use_alpha, use_beta>::apply(y, AA, x, alpha, beta);
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
  
  template<typename eT, typename TA>
  inline
  static
  void
  apply_blas_type( eT* y, const TA& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    arma_extra_debug_sigprint();
    
    //const uword threshold = (is_cx<eT>::yes) ? 16u : 64u;
    const uword threshold = (is_cx<eT>::yes) ? 64u : 100u;
    
    if(A.n_elem <= threshold)
      {
      gemv_emul<do_trans_A, use_alpha, use_beta>::apply(y,A,x,alpha,beta);
      }
    else
      {
      #if defined(ARMA_USE_ATLAS)
        {
        if(is_cx<eT>::no)
          {
          // use gemm() instead of gemv() to work around a speed issue in Atlas 3.8.4
          
          arma_extra_debug_print("atlas::cblas_gemm()");
          
          atlas::cblas_gemm<eT>
            (
            atlas::CblasColMajor,
            (do_trans_A) ? ( is_cx<eT>::yes ? CblasConjTrans : atlas::CblasTrans ) : atlas::CblasNoTrans,
            atlas::CblasNoTrans,
            (do_trans_A) ? A.n_cols : A.n_rows,
            1,
            (do_trans_A) ? A.n_rows : A.n_cols,
            (use_alpha) ? alpha : eT(1),
            A.mem,
            A.n_rows,
            x,
            (do_trans_A) ? A.n_rows : A.n_cols,
            (use_beta) ? beta : eT(0),
            y,
            (do_trans_A) ? A.n_cols : A.n_rows
            );
          }
        else
          {
          arma_extra_debug_print("atlas::cblas_gemv()");
          
          atlas::cblas_gemv<eT>
            (
            atlas::CblasColMajor,
            (do_trans_A) ? ( is_cx<eT>::yes ? CblasConjTrans : atlas::CblasTrans ) : atlas::CblasNoTrans,
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
        }
      #elif defined(ARMA_USE_BLAS)
        {
        arma_extra_debug_print("blas::gemv()");
        
        const char      trans_A     = (do_trans_A) ? ( is_cx<eT>::yes ? 'C' : 'T' ) : 'N';
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
  
  
  
  template<typename eT, typename TA>
  arma_inline
  static
  void
  apply( eT* y, const TA& A, const eT* x, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    gemv_emul<do_trans_A, use_alpha, use_beta>::apply(y,A,x,alpha,beta);
    }
  
  
  
  template<typename TA>
  arma_inline
  static
  void
  apply
    (
          float* y,
    const TA&    A,
    const float* x,
    const float  alpha = float(1),
    const float  beta  = float(0)
    )
    {
    gemv<do_trans_A, use_alpha, use_beta>::apply_blas_type(y,A,x,alpha,beta);
    }
  
  
  
  template<typename TA>
  arma_inline
  static
  void
  apply
    (
          double* y,
    const TA&     A,
    const double* x,
    const double  alpha = double(1),
    const double  beta  = double(0)
    )
    {
    gemv<do_trans_A, use_alpha, use_beta>::apply_blas_type(y,A,x,alpha,beta);
    }
  
  
  
  template<typename TA>
  arma_inline
  static
  void
  apply
    (
          std::complex<float>* y,
    const TA&                  A,
    const std::complex<float>* x,
    const std::complex<float>  alpha = std::complex<float>(1),
    const std::complex<float>  beta  = std::complex<float>(0)
    )
    {
    gemv<do_trans_A, use_alpha, use_beta>::apply_blas_type(y,A,x,alpha,beta);
    }


  
  template<typename TA>
  arma_inline
  static
  void
  apply
    (
          std::complex<double>* y,
    const TA&                   A,
    const std::complex<double>* x,
    const std::complex<double>  alpha = std::complex<double>(1),
    const std::complex<double>  beta  = std::complex<double>(0)
    )
    {
    gemv<do_trans_A, use_alpha, use_beta>::apply_blas_type(y,A,x,alpha,beta);
    }


  
  };


//! @}
