// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Edmund Highcock (edmund dot highcock at merton dot ox dot ac dot uk)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup auxlib
//! @{


//! immediate matrix inverse
template<typename eT>
inline
bool
auxlib::inv_noalias(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  switch(X.n_rows)
    {
    case 1:
      {
      out.set_size(1,1);
      out[0] = eT(1) / X[0];
      };
      break;
      
    case 2:
      {
      out.set_size(2,2);
      
      const eT a = X.at(0,0);
      const eT b = X.at(0,1);
      const eT c = X.at(1,0);
      const eT d = X.at(1,1);
      
      const eT k = eT(1) / (a*d - b*c);
      
      out.at(0,0) =  d*k;
      out.at(0,1) = -b*k;
      out.at(1,0) = -c*k;
      out.at(1,1) =  a*k;
      };
      break;
    
    case 3:
      {
      out.set_size(3,3);
      
      const eT* X_col0 = X.colptr(0);
      const eT a11 = X_col0[0];
      const eT a21 = X_col0[1];
      const eT a31 = X_col0[2];
      
      const eT* X_col1 = X.colptr(1);
      const eT a12 = X_col1[0];
      const eT a22 = X_col1[1];
      const eT a32 = X_col1[2];
      
      const eT* X_col2 = X.colptr(2);
      const eT a13 = X_col2[0];
      const eT a23 = X_col2[1];
      const eT a33 = X_col2[2];
      
      const eT k = eT(1) / ( a11*(a33*a22 - a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12 - a22*a13) );
      
      
      eT* out_col0 = out.colptr(0);
      out_col0[0] =  (a33*a22 - a32*a23) * k;
      out_col0[1] = -(a33*a21 - a31*a23) * k;
      out_col0[2] =  (a32*a21 - a31*a22) * k;
      
      eT* out_col1 = out.colptr(1);
      out_col1[0] = -(a33*a12 - a32*a13) * k;
      out_col1[1] =  (a33*a11 - a31*a13) * k;
      out_col1[2] = -(a32*a11 - a31*a12) * k;
      
      eT* out_col2 = out.colptr(2);
      out_col2[0] =  (a23*a12 - a22*a13) * k;
      out_col2[1] = -(a23*a11 - a21*a13) * k;
      out_col2[2] =  (a22*a11 - a21*a12) * k;
      };
      break;
      
    case 4:
      {
      out.set_size(4,4);
      
      out.at(0,0) = X.at(1,2)*X.at(2,3)*X.at(3,1) - X.at(1,3)*X.at(2,2)*X.at(3,1) + X.at(1,3)*X.at(2,1)*X.at(3,2) - X.at(1,1)*X.at(2,3)*X.at(3,2) - X.at(1,2)*X.at(2,1)*X.at(3,3) + X.at(1,1)*X.at(2,2)*X.at(3,3);
      out.at(1,0) = X.at(1,3)*X.at(2,2)*X.at(3,0) - X.at(1,2)*X.at(2,3)*X.at(3,0) - X.at(1,3)*X.at(2,0)*X.at(3,2) + X.at(1,0)*X.at(2,3)*X.at(3,2) + X.at(1,2)*X.at(2,0)*X.at(3,3) - X.at(1,0)*X.at(2,2)*X.at(3,3);
      out.at(2,0) = X.at(1,1)*X.at(2,3)*X.at(3,0) - X.at(1,3)*X.at(2,1)*X.at(3,0) + X.at(1,3)*X.at(2,0)*X.at(3,1) - X.at(1,0)*X.at(2,3)*X.at(3,1) - X.at(1,1)*X.at(2,0)*X.at(3,3) + X.at(1,0)*X.at(2,1)*X.at(3,3);
      out.at(3,0) = X.at(1,2)*X.at(2,1)*X.at(3,0) - X.at(1,1)*X.at(2,2)*X.at(3,0) - X.at(1,2)*X.at(2,0)*X.at(3,1) + X.at(1,0)*X.at(2,2)*X.at(3,1) + X.at(1,1)*X.at(2,0)*X.at(3,2) - X.at(1,0)*X.at(2,1)*X.at(3,2);
      
      out.at(0,1) = X.at(0,3)*X.at(2,2)*X.at(3,1) - X.at(0,2)*X.at(2,3)*X.at(3,1) - X.at(0,3)*X.at(2,1)*X.at(3,2) + X.at(0,1)*X.at(2,3)*X.at(3,2) + X.at(0,2)*X.at(2,1)*X.at(3,3) - X.at(0,1)*X.at(2,2)*X.at(3,3);
      out.at(1,1) = X.at(0,2)*X.at(2,3)*X.at(3,0) - X.at(0,3)*X.at(2,2)*X.at(3,0) + X.at(0,3)*X.at(2,0)*X.at(3,2) - X.at(0,0)*X.at(2,3)*X.at(3,2) - X.at(0,2)*X.at(2,0)*X.at(3,3) + X.at(0,0)*X.at(2,2)*X.at(3,3);
      out.at(2,1) = X.at(0,3)*X.at(2,1)*X.at(3,0) - X.at(0,1)*X.at(2,3)*X.at(3,0) - X.at(0,3)*X.at(2,0)*X.at(3,1) + X.at(0,0)*X.at(2,3)*X.at(3,1) + X.at(0,1)*X.at(2,0)*X.at(3,3) - X.at(0,0)*X.at(2,1)*X.at(3,3);
      out.at(3,1) = X.at(0,1)*X.at(2,2)*X.at(3,0) - X.at(0,2)*X.at(2,1)*X.at(3,0) + X.at(0,2)*X.at(2,0)*X.at(3,1) - X.at(0,0)*X.at(2,2)*X.at(3,1) - X.at(0,1)*X.at(2,0)*X.at(3,2) + X.at(0,0)*X.at(2,1)*X.at(3,2);
      
      out.at(0,2) = X.at(0,2)*X.at(1,3)*X.at(3,1) - X.at(0,3)*X.at(1,2)*X.at(3,1) + X.at(0,3)*X.at(1,1)*X.at(3,2) - X.at(0,1)*X.at(1,3)*X.at(3,2) - X.at(0,2)*X.at(1,1)*X.at(3,3) + X.at(0,1)*X.at(1,2)*X.at(3,3);
      out.at(1,2) = X.at(0,3)*X.at(1,2)*X.at(3,0) - X.at(0,2)*X.at(1,3)*X.at(3,0) - X.at(0,3)*X.at(1,0)*X.at(3,2) + X.at(0,0)*X.at(1,3)*X.at(3,2) + X.at(0,2)*X.at(1,0)*X.at(3,3) - X.at(0,0)*X.at(1,2)*X.at(3,3);
      out.at(2,2) = X.at(0,1)*X.at(1,3)*X.at(3,0) - X.at(0,3)*X.at(1,1)*X.at(3,0) + X.at(0,3)*X.at(1,0)*X.at(3,1) - X.at(0,0)*X.at(1,3)*X.at(3,1) - X.at(0,1)*X.at(1,0)*X.at(3,3) + X.at(0,0)*X.at(1,1)*X.at(3,3);
      out.at(3,2) = X.at(0,2)*X.at(1,1)*X.at(3,0) - X.at(0,1)*X.at(1,2)*X.at(3,0) - X.at(0,2)*X.at(1,0)*X.at(3,1) + X.at(0,0)*X.at(1,2)*X.at(3,1) + X.at(0,1)*X.at(1,0)*X.at(3,2) - X.at(0,0)*X.at(1,1)*X.at(3,2);
      
      out.at(0,3) = X.at(0,3)*X.at(1,2)*X.at(2,1) - X.at(0,2)*X.at(1,3)*X.at(2,1) - X.at(0,3)*X.at(1,1)*X.at(2,2) + X.at(0,1)*X.at(1,3)*X.at(2,2) + X.at(0,2)*X.at(1,1)*X.at(2,3) - X.at(0,1)*X.at(1,2)*X.at(2,3);
      out.at(1,3) = X.at(0,2)*X.at(1,3)*X.at(2,0) - X.at(0,3)*X.at(1,2)*X.at(2,0) + X.at(0,3)*X.at(1,0)*X.at(2,2) - X.at(0,0)*X.at(1,3)*X.at(2,2) - X.at(0,2)*X.at(1,0)*X.at(2,3) + X.at(0,0)*X.at(1,2)*X.at(2,3);
      out.at(2,3) = X.at(0,3)*X.at(1,1)*X.at(2,0) - X.at(0,1)*X.at(1,3)*X.at(2,0) - X.at(0,3)*X.at(1,0)*X.at(2,1) + X.at(0,0)*X.at(1,3)*X.at(2,1) + X.at(0,1)*X.at(1,0)*X.at(2,3) - X.at(0,0)*X.at(1,1)*X.at(2,3);
      out.at(3,3) = X.at(0,1)*X.at(1,2)*X.at(2,0) - X.at(0,2)*X.at(1,1)*X.at(2,0) + X.at(0,2)*X.at(1,0)*X.at(2,1) - X.at(0,0)*X.at(1,2)*X.at(2,1) - X.at(0,1)*X.at(1,0)*X.at(2,2) + X.at(0,0)*X.at(1,1)*X.at(2,2);      
      
      out /= det(X);
      };
      break;
      
    default:
      {
      #if defined(ARMA_USE_ATLAS)
        {
        out = X;
        podarray<int> ipiv(out.n_rows);
        
        int info = atlas::clapack_getrf(atlas::CblasColMajor, out.n_rows, out.n_cols, out.memptr(), out.n_rows, ipiv.memptr());
        
        if(info == 0)
          {
          info = atlas::clapack_getri(atlas::CblasColMajor, out.n_rows, out.memptr(), out.n_rows, ipiv.memptr());
          }
        
        return (info == 0);
        }
      #elif defined(ARMA_USE_LAPACK)
        {
        out = X;
        
        int n_rows = out.n_rows;
        int n_cols = out.n_cols;
        int info   = 0;
        
        podarray<int> ipiv(out.n_rows);
        
        // 84 was empirically found -- it is the maximum value suggested by LAPACK (as provided by ATLAS v3.6)
        // based on tests with various matrix types on 32-bit and 64-bit machines
        //
        // the "work" array is deliberately long so that a secondary (time-consuming)
        // memory allocation is avoided, if possible
        
        int work_len = (std::max)(1, n_rows*84);
        podarray<eT> work(work_len);
        
        lapack::getrf_(&n_rows, &n_cols, out.memptr(), &n_rows, ipiv.memptr(), &info);
        
        if(info == 0)
          {
          // query for optimum size of work_len
          
          int work_len_tmp = -1;
          lapack::getri_(&n_rows, out.memptr(), &n_rows, ipiv.memptr(), work.memptr(), &work_len_tmp, &info);
          
          if(info == 0)
            {
            int proposed_work_len = static_cast<int>(access::tmp_real(work[0]));
            
            // if necessary, allocate more memory
            if(work_len < proposed_work_len)
              {
              work_len = proposed_work_len;
              work.set_size(work_len);
              }
            }
          
          lapack::getri_(&n_rows, out.memptr(), &n_rows, ipiv.memptr(), work.memptr(), &work_len, &info);
          }
        
        return (info == 0);
        }
      #else
        {
        arma_stop("inv(): need ATLAS or LAPACK");
        }
      #endif
      };
    }
    
  return true;
  }
  
  
  
//! immediate inplace matrix inverse
template<typename eT>
inline
bool
auxlib::inv_inplace(Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  // for more info, see:
  // http://www.dr-lex.34sp.com/random/matrix_inv.html
  // http://www.cvl.iis.u-tokyo.ac.jp/~miyazaki/tech/teche23.html
  // http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
  // http://www.geometrictools.com//LibFoundation/Mathematics/Wm4Matrix4.inl
  
  switch(X.n_rows)
    {
    case 1:
      {
      X[0] = eT(1) / X[0];
      };
      break;
      
    case 2:
      {
      const eT a = X.at(0,0);
      const eT b = X.at(0,1);
      const eT c = X.at(1,0);
      const eT d = X.at(1,1);
      
      const eT k = eT(1) / (a*d - b*c);
      
      X.at(0,0) =  d*k;
      X.at(0,1) = -b*k;
      X.at(1,0) = -c*k;
      X.at(1,1) =  a*k;
      };
      break;
    
    case 3:
      {
      eT* X_col0 = X.colptr(0);
      eT* X_col1 = X.colptr(1);
      eT* X_col2 = X.colptr(2);
      
      const eT a11 = X_col0[0];
      const eT a21 = X_col0[1];
      const eT a31 = X_col0[2];
      
      const eT a12 = X_col1[0];
      const eT a22 = X_col1[1];
      const eT a32 = X_col1[2];
      
      const eT a13 = X_col2[0];
      const eT a23 = X_col2[1];
      const eT a33 = X_col2[2];
      
      const eT k = eT(1) / ( a11*(a33*a22 - a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12 - a22*a13) );
      
      X_col0[0] =  (a33*a22 - a32*a23) * k;
      X_col0[1] = -(a33*a21 - a31*a23) * k;
      X_col0[2] =  (a32*a21 - a31*a22) * k;
      
      X_col1[0] = -(a33*a12 - a32*a13) * k;
      X_col1[1] =  (a33*a11 - a31*a13) * k;
      X_col1[2] = -(a32*a11 - a31*a12) * k;
      
      X_col2[0] =  (a23*a12 - a22*a13) * k;
      X_col2[1] = -(a23*a11 - a21*a13) * k;
      X_col2[2] =  (a22*a11 - a21*a12) * k;
      };
      break;
      
    case 4:
      {
      const Mat<eT> A(X);
      
      X.at(0,0) = A.at(1,2)*A.at(2,3)*A.at(3,1) - A.at(1,3)*A.at(2,2)*A.at(3,1) + A.at(1,3)*A.at(2,1)*A.at(3,2) - A.at(1,1)*A.at(2,3)*A.at(3,2) - A.at(1,2)*A.at(2,1)*A.at(3,3) + A.at(1,1)*A.at(2,2)*A.at(3,3);
      X.at(1,0) = A.at(1,3)*A.at(2,2)*A.at(3,0) - A.at(1,2)*A.at(2,3)*A.at(3,0) - A.at(1,3)*A.at(2,0)*A.at(3,2) + A.at(1,0)*A.at(2,3)*A.at(3,2) + A.at(1,2)*A.at(2,0)*A.at(3,3) - A.at(1,0)*A.at(2,2)*A.at(3,3);
      X.at(2,0) = A.at(1,1)*A.at(2,3)*A.at(3,0) - A.at(1,3)*A.at(2,1)*A.at(3,0) + A.at(1,3)*A.at(2,0)*A.at(3,1) - A.at(1,0)*A.at(2,3)*A.at(3,1) - A.at(1,1)*A.at(2,0)*A.at(3,3) + A.at(1,0)*A.at(2,1)*A.at(3,3);
      X.at(3,0) = A.at(1,2)*A.at(2,1)*A.at(3,0) - A.at(1,1)*A.at(2,2)*A.at(3,0) - A.at(1,2)*A.at(2,0)*A.at(3,1) + A.at(1,0)*A.at(2,2)*A.at(3,1) + A.at(1,1)*A.at(2,0)*A.at(3,2) - A.at(1,0)*A.at(2,1)*A.at(3,2);
      
      X.at(0,1) = A.at(0,3)*A.at(2,2)*A.at(3,1) - A.at(0,2)*A.at(2,3)*A.at(3,1) - A.at(0,3)*A.at(2,1)*A.at(3,2) + A.at(0,1)*A.at(2,3)*A.at(3,2) + A.at(0,2)*A.at(2,1)*A.at(3,3) - A.at(0,1)*A.at(2,2)*A.at(3,3);
      X.at(1,1) = A.at(0,2)*A.at(2,3)*A.at(3,0) - A.at(0,3)*A.at(2,2)*A.at(3,0) + A.at(0,3)*A.at(2,0)*A.at(3,2) - A.at(0,0)*A.at(2,3)*A.at(3,2) - A.at(0,2)*A.at(2,0)*A.at(3,3) + A.at(0,0)*A.at(2,2)*A.at(3,3);
      X.at(2,1) = A.at(0,3)*A.at(2,1)*A.at(3,0) - A.at(0,1)*A.at(2,3)*A.at(3,0) - A.at(0,3)*A.at(2,0)*A.at(3,1) + A.at(0,0)*A.at(2,3)*A.at(3,1) + A.at(0,1)*A.at(2,0)*A.at(3,3) - A.at(0,0)*A.at(2,1)*A.at(3,3);
      X.at(3,1) = A.at(0,1)*A.at(2,2)*A.at(3,0) - A.at(0,2)*A.at(2,1)*A.at(3,0) + A.at(0,2)*A.at(2,0)*A.at(3,1) - A.at(0,0)*A.at(2,2)*A.at(3,1) - A.at(0,1)*A.at(2,0)*A.at(3,2) + A.at(0,0)*A.at(2,1)*A.at(3,2);
      
      X.at(0,2) = A.at(0,2)*A.at(1,3)*A.at(3,1) - A.at(0,3)*A.at(1,2)*A.at(3,1) + A.at(0,3)*A.at(1,1)*A.at(3,2) - A.at(0,1)*A.at(1,3)*A.at(3,2) - A.at(0,2)*A.at(1,1)*A.at(3,3) + A.at(0,1)*A.at(1,2)*A.at(3,3);
      X.at(1,2) = A.at(0,3)*A.at(1,2)*A.at(3,0) - A.at(0,2)*A.at(1,3)*A.at(3,0) - A.at(0,3)*A.at(1,0)*A.at(3,2) + A.at(0,0)*A.at(1,3)*A.at(3,2) + A.at(0,2)*A.at(1,0)*A.at(3,3) - A.at(0,0)*A.at(1,2)*A.at(3,3);
      X.at(2,2) = A.at(0,1)*A.at(1,3)*A.at(3,0) - A.at(0,3)*A.at(1,1)*A.at(3,0) + A.at(0,3)*A.at(1,0)*A.at(3,1) - A.at(0,0)*A.at(1,3)*A.at(3,1) - A.at(0,1)*A.at(1,0)*A.at(3,3) + A.at(0,0)*A.at(1,1)*A.at(3,3);
      X.at(3,2) = A.at(0,2)*A.at(1,1)*A.at(3,0) - A.at(0,1)*A.at(1,2)*A.at(3,0) - A.at(0,2)*A.at(1,0)*A.at(3,1) + A.at(0,0)*A.at(1,2)*A.at(3,1) + A.at(0,1)*A.at(1,0)*A.at(3,2) - A.at(0,0)*A.at(1,1)*A.at(3,2);
      
      X.at(0,3) = A.at(0,3)*A.at(1,2)*A.at(2,1) - A.at(0,2)*A.at(1,3)*A.at(2,1) - A.at(0,3)*A.at(1,1)*A.at(2,2) + A.at(0,1)*A.at(1,3)*A.at(2,2) + A.at(0,2)*A.at(1,1)*A.at(2,3) - A.at(0,1)*A.at(1,2)*A.at(2,3);
      X.at(1,3) = A.at(0,2)*A.at(1,3)*A.at(2,0) - A.at(0,3)*A.at(1,2)*A.at(2,0) + A.at(0,3)*A.at(1,0)*A.at(2,2) - A.at(0,0)*A.at(1,3)*A.at(2,2) - A.at(0,2)*A.at(1,0)*A.at(2,3) + A.at(0,0)*A.at(1,2)*A.at(2,3);
      X.at(2,3) = A.at(0,3)*A.at(1,1)*A.at(2,0) - A.at(0,1)*A.at(1,3)*A.at(2,0) - A.at(0,3)*A.at(1,0)*A.at(2,1) + A.at(0,0)*A.at(1,3)*A.at(2,1) + A.at(0,1)*A.at(1,0)*A.at(2,3) - A.at(0,0)*A.at(1,1)*A.at(2,3);
      X.at(3,3) = A.at(0,1)*A.at(1,2)*A.at(2,0) - A.at(0,2)*A.at(1,1)*A.at(2,0) + A.at(0,2)*A.at(1,0)*A.at(2,1) - A.at(0,0)*A.at(1,2)*A.at(2,1) - A.at(0,1)*A.at(1,0)*A.at(2,2) + A.at(0,0)*A.at(1,1)*A.at(2,2);      
      
      X /= det(A);
      };
      break;
      
    default:
      {
      #if defined(ARMA_USE_ATLAS)
        {
        Mat<eT>& out = X;
        podarray<int> ipiv(out.n_rows);
        
        int info = atlas::clapack_getrf(atlas::CblasColMajor, out.n_rows, out.n_cols, out.memptr(), out.n_rows, ipiv.memptr());
        
        if(info == 0)
          {
          info = atlas::clapack_getri(atlas::CblasColMajor, out.n_rows, out.memptr(), out.n_rows, ipiv.memptr());
          }
        
        return (info == 0);
        }
      #elif defined(ARMA_USE_LAPACK)
        {
        Mat<eT>& out = X;
        
        int n_rows = out.n_rows;
        int n_cols = out.n_cols;
        int info   = 0;
        
        podarray<int> ipiv(out.n_rows);
        
        // 84 was empirically found -- it is the maximum value suggested by LAPACK (as provided by ATLAS v3.6)
        // based on tests with various matrix types on 32-bit and 64-bit machines
        //
        // the "work" array is deliberately long so that a secondary (time-consuming)
        // memory allocation is avoided, if possible
        
        int work_len = (std::max)(1, n_rows*84);
        podarray<eT> work(work_len);
        
        lapack::getrf_(&n_rows, &n_cols, out.memptr(), &n_rows, ipiv.memptr(), &info);
        
        if(info == 0)
          {
          // query for optimum size of work_len
          
          int work_len_tmp = -1;
          lapack::getri_(&n_rows, out.memptr(), &n_rows, ipiv.memptr(), work.memptr(), &work_len_tmp, &info);
          
          if(info == 0)
            {
            int proposed_work_len = static_cast<int>(access::tmp_real(work[0]));
            
            // if necessary, allocate more memory
            if(work_len < proposed_work_len)
              {
              work_len = proposed_work_len;
              work.set_size(work_len);
              }
            }
          
          lapack::getri_(&n_rows, out.memptr(), &n_rows, ipiv.memptr(), work.memptr(), &work_len, &info);
          }
        
        return (info == 0);
        }
      #else
        {
        arma_stop("inv(): need ATLAS or LAPACK");
        }
      #endif
      }
    
    }
  
  return true;
  }


//! immediate determinant of a matrix using ATLAS or LAPACK
template<typename eT>
inline
eT
auxlib::det(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  switch(X.n_rows)
    {
    case 0:
      return eT(0);
    
    case 1:
      return X[0];
    
    case 2:
      return (X.at(0,0)*X.at(1,1) - X.at(0,1)*X.at(1,0));
    
    case 3:
      {
      const eT* a_col0 = X.colptr(0);
      const eT  a11    = a_col0[0];
      const eT  a21    = a_col0[1];
      const eT  a31    = a_col0[2];
      
      const eT* a_col1 = X.colptr(1);
      const eT  a12    = a_col1[0];
      const eT  a22    = a_col1[1];
      const eT  a32    = a_col1[2];
      
      const eT* a_col2 = X.colptr(2);
      const eT  a13    = a_col2[0];
      const eT  a23    = a_col2[1];
      const eT  a33    = a_col2[2];
      
      return ( a11*(a33*a22 - a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12 - a22*a13) );
      
      // const double tmp1 = X.at(0,0) * X.at(1,1) * X.at(2,2);
      // const double tmp2 = X.at(0,1) * X.at(1,2) * X.at(2,0);
      // const double tmp3 = X.at(0,2) * X.at(1,0) * X.at(2,1);
      // const double tmp4 = X.at(2,0) * X.at(1,1) * X.at(0,2);
      // const double tmp5 = X.at(2,1) * X.at(1,2) * X.at(0,0);
      // const double tmp6 = X.at(2,2) * X.at(1,0) * X.at(0,1);
      // return (tmp1+tmp2+tmp3) - (tmp4+tmp5+tmp6);
      }
      
    case 4:
      {
      const eT val = \
          X.at(0,3) * X.at(1,2) * X.at(2,1) * X.at(3,0) \
        - X.at(0,2) * X.at(1,3) * X.at(2,1) * X.at(3,0) \
        - X.at(0,3) * X.at(1,1) * X.at(2,2) * X.at(3,0) \
        + X.at(0,1) * X.at(1,3) * X.at(2,2) * X.at(3,0) \
        + X.at(0,2) * X.at(1,1) * X.at(2,3) * X.at(3,0) \
        - X.at(0,1) * X.at(1,2) * X.at(2,3) * X.at(3,0) \
        - X.at(0,3) * X.at(1,2) * X.at(2,0) * X.at(3,1) \
        + X.at(0,2) * X.at(1,3) * X.at(2,0) * X.at(3,1) \
        + X.at(0,3) * X.at(1,0) * X.at(2,2) * X.at(3,1) \
        - X.at(0,0) * X.at(1,3) * X.at(2,2) * X.at(3,1) \
        - X.at(0,2) * X.at(1,0) * X.at(2,3) * X.at(3,1) \
        + X.at(0,0) * X.at(1,2) * X.at(2,3) * X.at(3,1) \
        + X.at(0,3) * X.at(1,1) * X.at(2,0) * X.at(3,2) \
        - X.at(0,1) * X.at(1,3) * X.at(2,0) * X.at(3,2) \
        - X.at(0,3) * X.at(1,0) * X.at(2,1) * X.at(3,2) \
        + X.at(0,0) * X.at(1,3) * X.at(2,1) * X.at(3,2) \
        + X.at(0,1) * X.at(1,0) * X.at(2,3) * X.at(3,2) \
        - X.at(0,0) * X.at(1,1) * X.at(2,3) * X.at(3,2) \
        - X.at(0,2) * X.at(1,1) * X.at(2,0) * X.at(3,3) \
        + X.at(0,1) * X.at(1,2) * X.at(2,0) * X.at(3,3) \
        + X.at(0,2) * X.at(1,0) * X.at(2,1) * X.at(3,3) \
        - X.at(0,0) * X.at(1,2) * X.at(2,1) * X.at(3,3) \
        - X.at(0,1) * X.at(1,0) * X.at(2,2) * X.at(3,3) \
        + X.at(0,0) * X.at(1,1) * X.at(2,2) * X.at(3,3) \
        ;
      
      return val;
      }
      
    default:  
      {
      #if defined(ARMA_USE_ATLAS)
        {
        Mat<eT> tmp = X;
        podarray<int> ipiv(tmp.n_rows);
        
        atlas::clapack_getrf(atlas::CblasColMajor, tmp.n_rows, tmp.n_cols, tmp.memptr(), tmp.n_rows, ipiv.memptr());
        
        // on output tmp appears to be L+U_alt, where U_alt is U with the main diagonal set to zero
        eT val = tmp.at(0,0);
        for(u32 i=1; i < tmp.n_rows; ++i)
          {
          val *= tmp.at(i,i);
          }
        
        int sign = +1;
        for(u32 i=0; i < tmp.n_rows; ++i)
          {
          if( int(i) != ipiv.mem[i] )  // NOTE: no adjustment required, as the clapack version of getrf() assumes counting from 0
            {
            sign *= -1;
            }
          }
        
        return val * eT(sign);
        }
      #elif defined(ARMA_USE_LAPACK)
        {
        Mat<eT> tmp = X;
        podarray<int> ipiv(tmp.n_rows);
        
        int info   = 0;
        int n_rows = int(tmp.n_rows);
        int n_cols = int(tmp.n_cols);
        
        lapack::getrf_(&n_rows, &n_cols, tmp.memptr(), &n_rows, ipiv.memptr(), &info);
        
        // on output tmp appears to be L+U_alt, where U_alt is U with the main diagonal set to zero
        eT val = tmp.at(0,0);
        for(u32 i=1; i < tmp.n_rows; ++i)
          {
          val *= tmp.at(i,i);
          }
        
        int sign = +1;
        for(u32 i=0; i < tmp.n_rows; ++i)
          {
          if( int(i) != (ipiv.mem[i] - 1) )  // NOTE: adjustment of -1 is required as Fortran counts from 1
            {
            sign *= -1;
            }
          }
        
        return val * eT(sign);
        }
      #else
        {
        arma_stop("det(): need ATLAS or LAPACK");
        return eT(0);
        }
      #endif
      }
    }
  }



//! immediate log determinant of a matrix using ATLAS or LAPACK
template<typename eT>
inline
void
auxlib::log_det(eT& out_val, typename get_pod_type<eT>::result& out_sign, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  #if defined(ARMA_USE_ATLAS)
    {
    Mat<eT> tmp = X;
    podarray<int> ipiv(tmp.n_rows);
    
    atlas::clapack_getrf(atlas::CblasColMajor, tmp.n_rows, tmp.n_cols, tmp.memptr(), tmp.n_rows, ipiv.memptr());
    
    // on output tmp appears to be L+U_alt, where U_alt is U with the main diagonal set to zero
    
    s32 sign = (is_complex<eT>::value == false) ? ( (access::tmp_real( tmp.at(0,0) ) < T(0)) ? -1 : +1 ) : +1;
    eT   val = (is_complex<eT>::value == false) ? std::log( (access::tmp_real( tmp.at(0,0) ) < T(0)) ? tmp.at(0,0)*T(-1) : tmp.at(0,0) ) : std::log( tmp.at(0,0) );
    
    for(u32 i=1; i < tmp.n_rows; ++i)
      {
      const eT x = tmp.at(i,i);
      
      sign *= (is_complex<eT>::value == false) ? ( (access::tmp_real(x) < T(0)) ? -1 : +1 ) : +1;
      val  += (is_complex<eT>::value == false) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x ) : std::log(x);
      }
    
    for(u32 i=0; i < tmp.n_rows; ++i)
      {
      if( int(i) != ipiv.mem[i] )  // NOTE: no adjustment required, as the clapack version of getrf() assumes counting from 0
        {
        sign *= -1;
        }
      }
    
    out_val  = val;
    out_sign = T(sign);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    Mat<eT> tmp = X;
    podarray<int> ipiv(tmp.n_rows);
    
    int info   = 0;
    int n_rows = int(tmp.n_rows);
    int n_cols = int(tmp.n_cols);
    
    lapack::getrf_(&n_rows, &n_cols, tmp.memptr(), &n_rows, ipiv.memptr(), &info);
    
    // on output tmp appears to be L+U_alt, where U_alt is U with the main diagonal set to zero
    
    s32 sign = (is_complex<eT>::value == false) ? ( (access::tmp_real( tmp.at(0,0) ) < T(0)) ? -1 : +1 ) : +1;
    eT   val = (is_complex<eT>::value == false) ? std::log( (access::tmp_real( tmp.at(0,0) ) < T(0)) ? tmp.at(0,0)*T(-1) : tmp.at(0,0) ) : std::log( tmp.at(0,0) );
    
    for(u32 i=1; i < tmp.n_rows; ++i)
      {
      const eT x = tmp.at(i,i);
      
      sign *= (is_complex<eT>::value == false) ? ( (access::tmp_real(x) < T(0)) ? -1 : +1 ) : +1;
      val  += (is_complex<eT>::value == false) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x ) : std::log(x);
      }
    
    for(u32 i=0; i < tmp.n_rows; ++i)
      {
      if( int(i) != (ipiv.mem[i] - 1) )  // NOTE: adjustment of -1 is required as Fortran counts from 1
        {
        sign *= -1;
        }
      }
    
    out_val  = val;
    out_sign = T(sign);
    }
  #else
    {
    arma_stop("log_det(): need ATLAS or LAPACK");
    
    out_val  = eT(0);
    out_sign =  T(0);
    }
  #endif
  }



//! immediate LU decomposition of a matrix using ATLAS or LAPACK
template<typename eT>
inline
void
auxlib::lu(Mat<eT>& L, Mat<eT>& U, podarray<int>& ipiv, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  U = X;
  
  #if defined(ARMA_USE_ATLAS) || defined(ARMA_USE_LAPACK)
    {
    
    #if defined(ARMA_USE_ATLAS)
      {
      ipiv.set_size(U.n_rows);
    
      //int info = 
      atlas::clapack_getrf(atlas::CblasColMajor, U.n_rows, U.n_cols, U.memptr(), U.n_rows, ipiv.memptr());
      }
    #elif defined(ARMA_USE_LAPACK)
      {
      ipiv.set_size(U.n_rows);
      int info = 0;
      
      int n_rows = U.n_rows;
      int n_cols = U.n_cols;
      
      lapack::getrf_(&n_rows, &n_cols, U.memptr(), &n_rows, ipiv.memptr(), &info);
      
      // take into account that Fortran counts from 1
      for(u32 i=0; i<U.n_rows; ++i)
        {
        ipiv[i] -= 1;
        }
      
      }
    #endif
    
    
    L.set_size(U.n_rows, U.n_rows);
    
    for(u32 col=0; col<L.n_cols; ++col)
      {
      
      for(u32 row=0; row<col; ++row)
        {
        L.at(row,col) = eT(0);
        }
      
      L.at(col,col) = eT(1);
      
      for(u32 row=col+1; row<L.n_rows; ++row)
        {
        L.at(row,col) = U.at(row,col);
        U.at(row,col) = eT(0);
        }
      
      }
    }
  #else
    {
    arma_stop("lu(): need ATLAS or LAPACK");
    }
  #endif
  
  }



template<typename eT>
inline
void
auxlib::lu(Mat<eT>& L, Mat<eT>& U, Mat<eT>& P, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  podarray<int> ipiv;
  auxlib::lu(L, U, ipiv, X);
  
  const u32 n = ipiv.n_elem;

  Mat<u32> P_tmp(n,n);
  Mat<u32> ident = eye< Mat<u32> >(n,n);
  
  for(u32 i=n; i>0; --i)
    {
    const u32 j = i-1;
    const u32 k = ipiv[j];
    
    ident.swap_rows(j,k);
    
    if(i == n)
      {
      P_tmp = ident;
      }
    else
      {
      P_tmp *= ident;
      }
    
    ident.swap_rows(j,k);
    }
  
  P = conv_to< Mat<eT> >::from(P_tmp);
  }



template<typename eT>
inline
void
auxlib::lu(Mat<eT>& L, Mat<eT>& U, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  podarray<int> ipiv;
  auxlib::lu(L, U, ipiv, X);
  }
  


//! immediate eigenvalues of a symmetric real matrix using LAPACK
template<typename eT>
inline
void
auxlib::eig_sym(Col<eT>& eigval, const Mat<eT>& A_orig)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    const unwrap_check<Mat<eT> > tmp(A_orig, eigval);
    const Mat<eT>& A = tmp.M;
    
    arma_debug_check( (A.n_rows != A.n_cols), "eig_sym(): given matrix is not square");
    
    // rudimentary "better-than-nothing" test for symmetry
    //arma_debug_check( (A.at(A.n_rows-1, 0) != A.at(0, A.n_cols-1)), "auxlib::eig(): given matrix is not symmetric" );

    char jobz  = 'N';
    char uplo  = 'U';
    
    int n_rows = A.n_rows;
    int lwork  = (std::max)(1,3*n_rows-1);
    
    eigval.set_size(n_rows);
    podarray<eT> work(lwork);
  
    Mat<eT> A_copy = A;
    int info;

    arma_extra_debug_print("lapack::syev_()");
    lapack::syev_(&jobz, &uplo, &n_rows, A_copy.memptr(), &n_rows, eigval.memptr(), work.memptr(), &lwork, &info);
    }
  #else
    {
    arma_stop("eig_sym(): need LAPACK");
    }
  #endif
  }



//! immediate eigenvalues of a hermitian complex matrix using LAPACK
template<typename T> 
inline
void
auxlib::eig_sym(Col<T>& eigval, const Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();

  typedef typename std::complex<T> eT;
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_check( (A.n_rows != A.n_cols), "eig_sym(): given matrix is not hermitian");

    char jobz  = 'N'; 
    char uplo  = 'U';

    int n_rows = A.n_rows;
    int lda    = A.n_rows;
    int lwork  = (std::max)(1, 2*n_rows - 1);  // TODO: automatically find best size of lwork

    eigval.set_size(n_rows);

    podarray<eT> work(lwork);
    podarray<T>  rwork( (std::max)(1, 3*n_rows - 2) );
  
    Mat<eT> A_copy = A;
    int info;
  
    arma_extra_debug_print("lapack::heev_()");
    lapack::heev_(&jobz, &uplo, &n_rows, A_copy.memptr(), &lda, eigval.memptr(), work.memptr(), &lwork, rwork.memptr(), &info);
    }
  #else
    {
    arma_stop("eig_sym(): need LAPACK");
    }
  #endif
  }



//! immediate eigenvalues and eigenvectors of a symmetric real matrix using LAPACK
template<typename eT>
inline
void
auxlib::eig_sym(Col<eT>& eigval, Mat<eT>& eigvec, const Mat<eT>& A_orig)
  {
  arma_extra_debug_sigprint();
  
  // TODO: check for aliasing
  
  #if defined(ARMA_USE_LAPACK)
    {
    const unwrap_check< Mat<eT> > tmp1(A_orig, eigval);
    const Mat<eT>& A_tmp = tmp1.M;
    
    const unwrap_check< Mat<eT> > tmp2(A_tmp, eigvec);
    const Mat<eT>& A = tmp2.M;
    
    arma_debug_check( (A.n_rows != A.n_cols), "eig_sym(): given matrix is not square" );
    
    // rudimentary "better-than-nothing" test for symmetry
    //arma_debug_check( (A.at(A.n_rows-1, 0) != A.at(0, A.n_cols-1)), "auxlib::eig(): given matrix is not symmetric" );
    
    
    char jobz  = 'V';
    char uplo  = 'U';
    
    int n_rows = A.n_rows;
    int lwork  = (std::max)(1, 3*n_rows-1);
    
    eigval.set_size(n_rows);
    podarray<eT> work(lwork);
  
    eigvec = A;
    int info;
    
    arma_extra_debug_print("lapack::syev_()");
    lapack::syev_(&jobz, &uplo, &n_rows, eigvec.memptr(), &n_rows, eigval.memptr(), work.memptr(), &lwork, &info);
    }
  #else
    {
    arma_stop("eig_sym(): need LAPACK");
    }
  #endif
  
  }



//! immediate eigenvalues and eigenvectors of a hermitian complex matrix using LAPACK
template<typename T>
inline
void
auxlib::eig_sym(Col<T>& eigval, Mat< std::complex<T> >& eigvec, const Mat< std::complex<T> >& A_orig) 
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;

  #if defined(ARMA_USE_LAPACK)
    {
    const unwrap_check< Mat<eT> > tmp(A_orig, eigvec);
    const Mat<eT>& A = tmp.M;
    
    arma_debug_check( (A.n_rows != A.n_cols), "eig_sym(): given matrix is not hermitian" );
    
    char jobz  = 'V';
    char uplo  = 'U';

    int n_rows = A.n_rows;
    int lda    = A.n_rows;
    int lwork  = (std::max)(1, 2*n_rows - 1);  // TODO: automatically find best size of lwork
    
    eigval.set_size(n_rows);

    podarray<eT> work(lwork);
    podarray<T>  rwork( (std::max)(1, 3*n_rows - 2) );
  
    eigvec = A;
    int info;
  
    arma_extra_debug_print("lapack::heev_()");
    lapack::heev_(&jobz, &uplo, &n_rows, eigvec.memptr(), &lda, eigval.memptr(), work.memptr(), &lwork, rwork.memptr(), &info);
    }
  #else
    {
    arma_stop("eig_sym(): need LAPACK");
    }
  #endif
  
  }



//! Eigenvalues and eigenvectors of a general square real matrix using LAPACK.
//! The argument 'side' specifies which eigenvectors should be calculated
//! (see code for mode details).
template<typename T>
inline
void
auxlib::eig_gen
  (
  Col< std::complex<T> >& eigval,
  Mat<T>& l_eigvec,
  Mat<T>& r_eigvec,
  const Mat<T>& A, 
  const char side
  )
  {
  arma_extra_debug_sigprint();

  // TODO: check for aliasing
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_check( (A.n_rows != A.n_cols), "eig_gen(): given matrix is not square" );
    
    char jobvl;
    char jobvr;

    switch(side)
      {
      case 'l':  // left
        jobvl = 'V';
        jobvr = 'N';
        break;
        
      case 'r':  // right
        jobvl = 'N';
        jobvr = 'V';
        break;
        
      case 'b':  // both
        jobvl = 'V';
        jobvr = 'V';
        break;
        
      case 'n':  // neither
        jobvl = 'N';
        jobvr = 'N';
        break;

      default:
        arma_stop("eig_gen(): parameter 'side' is invalid");
      }

       
    int n_rows = A.n_rows;
    int lda    = A.n_rows;
    int lwork  = (std::max)(1, 4*n_rows);  // TODO: automatically find best size of lwork
    
    eigval.set_size(n_rows);
    l_eigvec.set_size(n_rows, n_rows);
    r_eigvec.set_size(n_rows, n_rows);
    
    podarray<T> work(lwork);
    podarray<T> rwork( (std::max)(1, 3*n_rows) );
    
    podarray<T> wr(n_rows);
    podarray<T> wi(n_rows);
    
    Mat<T> A_copy = A;
    int info;
    
    arma_extra_debug_print("lapack::cx_geev_()");
    lapack::geev_(&jobvl, &jobvr, &n_rows, A_copy.memptr(), &lda, wr.memptr(), wi.memptr(), l_eigvec.memptr(), &n_rows, r_eigvec.memptr(), &n_rows, work.memptr(), &lwork, &info);
    
    
    eigval.set_size(n_rows);
    for(u32 i=0; i<u32(n_rows); ++i)
      {
      eigval[i] = std::complex<T>(wr[i], wi[i]);
      }
    
    }
  #else
    {
    arma_stop("eig_gen(): need LAPACK");
    }
  #endif
  
  }





//! Eigenvalues and eigenvectors of a general square complex matrix using LAPACK
//! The argument 'side' specifies which eigenvectors should be calculated
//! (see code for mode details).
template<typename T>
inline
void
auxlib::eig_gen
  (
  Col< std::complex<T> >& eigval,
  Mat< std::complex<T> >& l_eigvec, 
  Mat< std::complex<T> >& r_eigvec, 
  const Mat< std::complex<T> >& A, 
  const char side
  )
  {
  arma_extra_debug_sigprint();

  // TODO: check for aliasing
  
  typedef typename std::complex<T> eT;

  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_check( (A.n_rows != A.n_cols), "eig_gen(): given matrix is not square" );
    
    char jobvl;
    char jobvr;

    switch(side)
      {
      case 'l':  // left
        jobvl = 'V';
        jobvr = 'N';
        break;
        
      case 'r':  // right
        jobvl = 'N';
        jobvr = 'V';
        break;
        
      case 'b':  // both
        jobvl = 'V';
        jobvr = 'V';
        break;
        
      case 'n':  // neither
        jobvl = 'N';
        jobvr = 'N';
        break;

      default:
        arma_stop("eig_gen(): parameter 'side' is invalid");
      }
    
       
    int n_rows = A.n_rows;
    int lda    = A.n_rows;
    int lwork  = (std::max)(1, 4*n_rows);  // TODO: automatically find best size of lwork
    
    eigval.set_size(n_rows);
    l_eigvec.set_size(n_rows, n_rows);
    r_eigvec.set_size(n_rows, n_rows);
    
    podarray<eT> work(lwork);
    podarray<T>  rwork( (std::max)(1, 3*n_rows) );  // was 2,3
    
    Mat<eT> A_copy = A;
    int info;
    
    arma_extra_debug_print("lapack::cx_geev_()");
    lapack::cx_geev_(&jobvl, &jobvr, &n_rows, A_copy.memptr(), &lda, eigval.memptr(), l_eigvec.memptr(), &n_rows, r_eigvec.memptr(), &n_rows, work.memptr(), &lwork, rwork.memptr(), &info);
    }
  #else
    {
    arma_stop("eig_gen(): need LAPACK");
    }
  #endif
  
  }



template<typename eT> 
inline
bool
auxlib::chol(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    char uplo = 'U';
    int  n    = X.n_rows;
    int  lda  = X.n_rows;
    int  info;
    
    out = X;
    lapack::potrf_(&uplo, &n, out.memptr(), &lda, &info);
    
    for(u32 col=0; col<X.n_rows; ++col)
      {
      eT* colptr = out.colptr(col);
      for(u32 row=col+1; row<X.n_rows; ++row)
        {
        colptr[row] = eT(0);
        }
      }
    
    return (info == 0);
    }
  #else
    {
    arma_stop("chol(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool 
auxlib::qr(Mat<eT>& Q, Mat<eT>& R, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    int m            = static_cast<int>(X.n_rows);
    int n            = static_cast<int>(X.n_cols);
    int work_len     = (std::max)(1,n);
    int work_len_tmp;
    int k            = (std::min)(m,n);
    int info;
    
    podarray<eT> tau(k);
    podarray<eT> work(work_len);
    
    R = X;
    
    // query for the optimum value of work_len
    work_len_tmp = -1;
    lapack::geqrf_(&m, &n, R.memptr(), &m, tau.memptr(), work.memptr(), &work_len_tmp, &info);
    
    if(info == 0)
      {
      work_len = static_cast<int>(access::tmp_real(work[0]));
      work.set_size(work_len);
      }
    
    lapack::geqrf_(&m, &n, R.memptr(), &m, tau.memptr(), work.memptr(), &work_len, &info);
    
    Q.set_size(X.n_rows, X.n_rows);
    
          eT* Q_mem = Q.memptr();
    const eT* R_mem = R.mem;
    
    const u32 n_elem_copy = (std::min)(Q.n_elem, R.n_elem);
    for(u32 i=0; i < n_elem_copy; ++i)
      {
      Q_mem[i] = R_mem[i];
      }
    
    
    // construct R
    for(u32 row=0; row < R.n_rows; ++row)
      {
      const u32 n_elem_tmp = (std::min)(row, R.n_cols);
      for(u32 col=0; col < n_elem_tmp; ++col)
        {
        R.at(row,col) = eT(0);
        }
      }
      
    
    if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
      {
      // query for the optimum value of work_len
      work_len_tmp = -1;
      lapack::orgqr_(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &work_len_tmp, &info);
      
      if(info == 0)
        {
        work_len = static_cast<int>(access::tmp_real(work[0]));
        work.set_size(work_len);
        }
      
      lapack::orgqr_(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &work_len, &info);
      }
    else
    if( (is_supported_complex_float<eT>::value == true) || (is_supported_complex_double<eT>::value == true) )
      {
      // query for the optimum value of work_len
      work_len_tmp = -1;
      lapack::ungqr_(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &work_len_tmp, &info);
      
      if(info == 0)
        {
        work_len = static_cast<int>(access::tmp_real(work[0]));
        work.set_size(work_len);
        }
      
      lapack::ungqr_(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &work_len, &info);
      }
    
    return (info == 0);
    }
  #else
    {
    arma_stop("qr(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename eT> 
inline
bool
auxlib::svd(Col<eT>& S, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A = X;
    
    Mat<eT> U(1, 1);
    Mat<eT> V(1, A.n_cols);
    
    char jobu  = 'N';
    char jobvt = 'N';
    
    int  m     = A.n_rows;
    int  n     = A.n_cols;
    int  lda   = A.n_rows;
    int  ldu   = U.n_rows;
    int  ldvt  = V.n_rows;
    int  lwork = 2 * (std::max)(1, (std::max)( (3*(std::min)(m,n) + (std::max)(m,n)), 5*(std::min)(m,n) ) );
    int  info;
    
    S.set_size( (std::min)(m, n) );
    
    podarray<eT> work(lwork);
  
  
    // let gesvd_() calculate the optimum size of the workspace
    int lwork_tmp = -1;
    
    lapack::gesvd_<eT>
      (
      &jobu, &jobvt,
      &m,&n,
      A.memptr(), &lda,
      S.memptr(),
      U.memptr(), &ldu,
      V.memptr(), &ldvt,
      work.memptr(), &lwork_tmp,
      &info
      );
    
    if(info == 0)
      {
      int proposed_lwork = static_cast<int>(work[0]);
      
      if(proposed_lwork > lwork)
        {
        lwork = proposed_lwork;
        work.set_size(lwork);
        }
      
      lapack::gesvd_<eT>
        (
        &jobu, &jobvt,
        &m, &n,
        A.memptr(), &lda,
        S.memptr(),
        U.memptr(), &ldu,
        V.memptr(), &ldvt,
        work.memptr(), &lwork,
        &info
        );
      }
    
    return (info == 0);
    }
  #else
    {
    arma_stop("svd(): need LAPACK");
    return false;
    }
  #endif
  }


  
template<typename T>
inline
bool
auxlib::svd(Col<T>& S, const Mat< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A = X;
    
    Mat<eT> U(1, 1);
    Mat<eT> V(1, A.n_cols);
    
    char jobu  = 'N';
    char jobvt = 'N';
    
    int  m     = A.n_rows;
    int  n     = A.n_cols;
    int  lda   = A.n_rows;
    int  ldu   = U.n_rows;
    int  ldvt  = V.n_rows;
    int  lwork = 2 * (std::max)(1, 2*(std::min)(m,n)+(std::max)(m,n) );
    int  info;
    
    S.set_size( (std::min)(m,n) );
    
    podarray<eT> work(lwork);
    podarray<T>  rwork( 5*(std::min)(m,n) );
  
    // let gesvd_() calculate the optimum size of the workspace
    int lwork_tmp = -1;
    
    lapack::cx_gesvd_<T>
      (
      &jobu, &jobvt,
      &m, &n,
      A.memptr(), &lda,
      S.memptr(),
      U.memptr(), &ldu,
      V.memptr(), &ldvt,
      work.memptr(), &lwork_tmp,
      rwork.memptr(),
      &info
      );
    
    if(info == 0)
      {
      int proposed_lwork = static_cast<int>(real(work[0]));
      if(proposed_lwork > lwork)
        {
        lwork = proposed_lwork;
        work.set_size(lwork);
        }
      
      lapack::cx_gesvd_<T>
        (
        &jobu, &jobvt,
        &m, &n,
        A.memptr(), &lda,
        S.memptr(),
        U.memptr(), &ldu,
        V.memptr(), &ldvt,
        work.memptr(), &lwork,
        rwork.memptr(),
        &info
        );
      }
        
    return (info == 0);
    }
  #else
    {
    arma_stop("svd(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::svd(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A = X;
    
    U.set_size(A.n_rows, A.n_rows);
    V.set_size(A.n_cols, A.n_cols);
    
    char jobu  = 'A';
    char jobvt = 'A';
    
    int  m     = A.n_rows;
    int  n     = A.n_cols;
    int  lda   = A.n_rows;
    int  ldu   = U.n_rows;
    int  ldvt  = V.n_rows;
    int  lwork = 2 * (std::max)(1, (std::max)( (3*(std::min)(m,n) + (std::max)(m,n)), 5*(std::min)(m,n) ) );
    int  info;
    
    
    S.set_size( (std::min)(m,n) );
    podarray<eT> work(lwork);
  
    // let gesvd_() calculate the optimum size of the workspace
    int lwork_tmp = -1;
    
    lapack::gesvd_<eT>
      (
      &jobu, &jobvt,
      &m, &n,
      A.memptr(), &lda,
      S.memptr(),
      U.memptr(), &ldu,
      V.memptr(), &ldvt,
      work.memptr(), &lwork_tmp,
      &info
      );
    
    if(info == 0)
      {
      int proposed_lwork = static_cast<int>(work[0]);
      if(proposed_lwork > lwork)
        {
        lwork = proposed_lwork;
        work.set_size(lwork);
        }
      
      lapack::gesvd_<eT>
        (
        &jobu, &jobvt,
        &m, &n,
        A.memptr(), &lda,
        S.memptr(),
        U.memptr(), &ldu,
        V.memptr(), &ldvt,
        work.memptr(), &lwork,
        &info
        );
      
      op_trans::apply(V,V);  // op_trans will work out that an in-place transpose can be done
      }
    
    return (info == 0);
    }
  #else
    {
    arma_stop("svd(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
auxlib::svd(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, const Mat< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A = X;
    
    U.set_size(A.n_rows, A.n_rows);
    V.set_size(A.n_cols, A.n_cols);
    
    char jobu  = 'A';
    char jobvt = 'A';
    
    int  m     = A.n_rows;
    int  n     = A.n_cols;
    int  lda   = A.n_rows;
    int  ldu   = U.n_rows;
    int  ldvt  = V.n_rows;
    int  lwork = 2 * (std::max)(1, 2*(std::min)(m,n)+(std::max)(m,n) );
    int  info;
    
    S.set_size( (std::min)(m,n) );
    
    podarray<eT> work(lwork);
    podarray<T>  rwork( 5*(std::min)(m,n) );
  
    // let gesvd_() calculate the optimum size of the workspace
    int lwork_tmp = -1;
    lapack::cx_gesvd_<T>
     (
     &jobu, &jobvt,
     &m, &n,
     A.memptr(), &lda,
     S.memptr(),
     U.memptr(), &ldu,
     V.memptr(), &ldvt,
     work.memptr(), &lwork_tmp,
     rwork.memptr(),
     &info
     );
    
    if(info == 0)
      {
      int proposed_lwork = static_cast<int>(real(work[0]));
      if(proposed_lwork > lwork)
        {
        lwork = proposed_lwork;
        work.set_size(lwork);
        }
      
      lapack::cx_gesvd_<T>
        (
        &jobu, &jobvt,
        &m, &n,
        A.memptr(), &lda,
        S.memptr(),
        U.memptr(), &ldu,
        V.memptr(), &ldvt,
        work.memptr(), &lwork,
        rwork.memptr(),
        &info
        );
      
      op_htrans::apply(V,V);  // op_htrans will work out that an in-place transpose can be done
      }
    
    return (info == 0);
    }
  #else
    {
    arma_stop("svd(): need LAPACK");
    return false;
    }
  #endif
  
  }


//! Solve a system of linear equations
//! Assumes that A.n_rows = A.n_cols
//! and B.n_rows = A.n_rows
template<typename eT>
inline
bool
auxlib::solve(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    int n    = A.n_rows;
    int lda  = A.n_rows;
    int ldb  = A.n_rows;
    int nrhs = B.n_cols;
    int info;
    
    podarray<int> ipiv(n);
    
    out = B;
    Mat<eT> A_copy = A;
  
    lapack::gesv_<eT>(&n, &nrhs, A_copy.memptr(), &lda, ipiv.memptr(), out.memptr(), &ldb, &info);
  
    return (info == 0);
    }
  #else
    {
    arma_stop("solve(): need LAPACK");
    return false;
    }
  #endif
  }


               
//! Solve an over-determined system.
//! Assumes that A.n_rows > A.n_cols
//! and B.n_rows = A.n_rows
template<typename eT>
inline
bool
auxlib::solve_od(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    char trans = 'N';
    
    int  m     = A.n_rows;
    int  n     = A.n_cols;
    int  lda   = A.n_rows;
    int  ldb   = A.n_rows;
    int  nrhs  = B.n_cols;
    int  lwork = n + (std::max)(n, nrhs);
    int  info;
    
    Mat<eT> A_copy = A;
    Mat<eT> tmp    = B;
    
    
    podarray<eT> work(lwork);
    
    arma_extra_debug_print("lapack::gels_()");
    
    // NOTE:
    // the dgels() function in the lapack library supplied by ATLAS 3.6
    // seems to have problems
    
    lapack::gels_<eT>
      (
      &trans, &m, &n, &nrhs,
      A_copy.memptr(), &lda,
      tmp.memptr(), &ldb,
      work.memptr(), &lwork,
      &info
      );
    
    arma_extra_debug_print("lapack::gels_() -- finished");
    
    out.set_size(A.n_cols, B.n_cols);
    
    for(u32 col=0; col<B.n_cols; ++col)
      {
      syslib::copy_elem( out.colptr(col), tmp.colptr(col), A.n_cols );
      }
    
    return (info == 0);
    }
  #else
    {
    arma_stop("auxlib::solve_od(): need LAPACK");
    return false;
    }
  #endif
  }



//! Solve an under-determined system.
//! Assumes that A.n_rows < A.n_cols
//! and B.n_rows = A.n_rows
template<typename eT>
inline
bool
auxlib::solve_ud(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    char trans = 'N';
    
    int  m     = A.n_rows;
    int  n     = A.n_cols;
    int  lda   = A.n_rows;
    int  ldb   = A.n_cols;
    int  nrhs  = B.n_cols;
    int  lwork = m + (std::max)(m,nrhs);
    int  info;
    
    
    Mat<eT> A_copy = A;
    
    Mat<eT> tmp;
    tmp.zeros(A.n_cols, B.n_cols);
    
    for(u32 col=0; col<B.n_cols; ++col)
      {
      eT* tmp_colmem = tmp.colptr(col);
      
      syslib::copy_elem( tmp_colmem, B.colptr(col), B.n_rows );
      
      for(u32 row=B.n_rows; row<A.n_cols; ++row)
        {
        tmp_colmem[row] = eT(0);
        }
      }
    
    podarray<eT> work(lwork);
    
    arma_extra_debug_print("lapack::gels_()");
    
    // NOTE:
    // the dgels() function in the lapack library supplied by ATLAS 3.6
    // seems to have problems
    
    lapack::gels_<eT>
      (
      &trans, &m, &n, &nrhs,
      A_copy.memptr(), &lda,
      tmp.memptr(), &ldb,
      work.memptr(), &lwork,
      &info
      );
    
    arma_extra_debug_print("lapack::gels_() -- finished");
    
    out.set_size(A.n_cols, B.n_cols);
    
    for(u32 col=0; col<B.n_cols; ++col)
      {
      syslib::copy_elem( out.colptr(col), tmp.colptr(col), A.n_cols );
      }
  
    return (info == 0);
    }
  #else
    {
    arma_stop("auxlib::solve_ud(): need LAPACK");
    return false;
    }
  #endif
  }



//! @}
