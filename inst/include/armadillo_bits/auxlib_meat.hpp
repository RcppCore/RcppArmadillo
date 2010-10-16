// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// Copyright (C) 2009      Edmund Highcock
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
template<typename eT, typename T1>
inline
bool
auxlib::inv(Mat<eT>& out, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  bool status;
  
  out = X.get_ref();
  
  arma_debug_check( (out.is_square() == false), "inv(): given matrix is not square" );
  
  const u32 N = out.n_rows;
  
  if(N <= 4)
    {
    status = auxlib::inv_inplace_tinymat(out, N);
    }
    
  if( (N > 4) || (status == false) )
    {
    status = auxlib::inv_inplace_lapack(out);
    }
  
  if(status == false)
    {
    arma_print("inv(): matrix appears to be singular" );
    out.reset();
    }
  
  return status;
  }



template<typename eT>
inline
bool
auxlib::inv(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (X.is_square() == false), "inv(): given matrix is not square" );
  
  bool status;
  
  const u32 N = X.n_rows;
  
  if(N <= 4)
    {
    status = (&out != &X) ? auxlib::inv_noalias_tinymat(out, X, N) : auxlib::inv_inplace_tinymat(out, N);
    }
  
  if( (N > 4) || (status == false) )
    {
    out = X;
    status = auxlib::inv_inplace_lapack(out);
    }
  
  if(status == false)
    {
    arma_print("inv(): matrix appears to be singular" );
    out.reset();
    }
  
  return status;
  }



template<typename eT>
inline
bool
auxlib::inv_noalias_tinymat(Mat<eT>& out, const Mat<eT>& X, const u32 N)
  {
  arma_extra_debug_sigprint();
  
  bool det_ok = true;
  
  out.set_size(N,N);
  
  switch(N)
    {
    case 1:
      {
      out[0] = eT(1) / X[0];
      };
      break;
      
    case 2:
      {
      const eT a = X.at(0,0);
      const eT b = X.at(0,1);
      const eT c = X.at(1,0);
      const eT d = X.at(1,1);
      
      const eT tmp_det = (a*d - b*c);
      
      if(tmp_det != eT(0))
        {
        out.at(0,0) =  d / tmp_det;
        out.at(0,1) = -b / tmp_det;
        out.at(1,0) = -c / tmp_det;
        out.at(1,1) =  a / tmp_det;
        }
      else
        {
        det_ok = false;
        }
      };
      break;
    
    case 3:
      {
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
      
      const eT tmp_det = a11*(a33*a22 - a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12 - a22*a13);
      
      if(tmp_det != eT(0))
        {
        eT* out_col0 = out.colptr(0);
        out_col0[0] =  (a33*a22 - a32*a23) / tmp_det;
        out_col0[1] = -(a33*a21 - a31*a23) / tmp_det;
        out_col0[2] =  (a32*a21 - a31*a22) / tmp_det;
        
        eT* out_col1 = out.colptr(1);
        out_col1[0] = -(a33*a12 - a32*a13) / tmp_det;
        out_col1[1] =  (a33*a11 - a31*a13) / tmp_det;
        out_col1[2] = -(a32*a11 - a31*a12) / tmp_det;
        
        eT* out_col2 = out.colptr(2);
        out_col2[0] =  (a23*a12 - a22*a13) / tmp_det;
        out_col2[1] = -(a23*a11 - a21*a13) / tmp_det;
        out_col2[2] =  (a22*a11 - a21*a12) / tmp_det;
        }
      else
        {
        det_ok = false;
        }
      };
      break;
      
    case 4:
      {
      const eT tmp_det = det(X);
      
      if(tmp_det != eT(0))
        {
        out.at(0,0) = ( X.at(1,2)*X.at(2,3)*X.at(3,1) - X.at(1,3)*X.at(2,2)*X.at(3,1) + X.at(1,3)*X.at(2,1)*X.at(3,2) - X.at(1,1)*X.at(2,3)*X.at(3,2) - X.at(1,2)*X.at(2,1)*X.at(3,3) + X.at(1,1)*X.at(2,2)*X.at(3,3) ) / tmp_det;
        out.at(1,0) = ( X.at(1,3)*X.at(2,2)*X.at(3,0) - X.at(1,2)*X.at(2,3)*X.at(3,0) - X.at(1,3)*X.at(2,0)*X.at(3,2) + X.at(1,0)*X.at(2,3)*X.at(3,2) + X.at(1,2)*X.at(2,0)*X.at(3,3) - X.at(1,0)*X.at(2,2)*X.at(3,3) ) / tmp_det;
        out.at(2,0) = ( X.at(1,1)*X.at(2,3)*X.at(3,0) - X.at(1,3)*X.at(2,1)*X.at(3,0) + X.at(1,3)*X.at(2,0)*X.at(3,1) - X.at(1,0)*X.at(2,3)*X.at(3,1) - X.at(1,1)*X.at(2,0)*X.at(3,3) + X.at(1,0)*X.at(2,1)*X.at(3,3) ) / tmp_det;
        out.at(3,0) = ( X.at(1,2)*X.at(2,1)*X.at(3,0) - X.at(1,1)*X.at(2,2)*X.at(3,0) - X.at(1,2)*X.at(2,0)*X.at(3,1) + X.at(1,0)*X.at(2,2)*X.at(3,1) + X.at(1,1)*X.at(2,0)*X.at(3,2) - X.at(1,0)*X.at(2,1)*X.at(3,2) ) / tmp_det;
        
        out.at(0,1) = ( X.at(0,3)*X.at(2,2)*X.at(3,1) - X.at(0,2)*X.at(2,3)*X.at(3,1) - X.at(0,3)*X.at(2,1)*X.at(3,2) + X.at(0,1)*X.at(2,3)*X.at(3,2) + X.at(0,2)*X.at(2,1)*X.at(3,3) - X.at(0,1)*X.at(2,2)*X.at(3,3) ) / tmp_det;
        out.at(1,1) = ( X.at(0,2)*X.at(2,3)*X.at(3,0) - X.at(0,3)*X.at(2,2)*X.at(3,0) + X.at(0,3)*X.at(2,0)*X.at(3,2) - X.at(0,0)*X.at(2,3)*X.at(3,2) - X.at(0,2)*X.at(2,0)*X.at(3,3) + X.at(0,0)*X.at(2,2)*X.at(3,3) ) / tmp_det;
        out.at(2,1) = ( X.at(0,3)*X.at(2,1)*X.at(3,0) - X.at(0,1)*X.at(2,3)*X.at(3,0) - X.at(0,3)*X.at(2,0)*X.at(3,1) + X.at(0,0)*X.at(2,3)*X.at(3,1) + X.at(0,1)*X.at(2,0)*X.at(3,3) - X.at(0,0)*X.at(2,1)*X.at(3,3) ) / tmp_det;
        out.at(3,1) = ( X.at(0,1)*X.at(2,2)*X.at(3,0) - X.at(0,2)*X.at(2,1)*X.at(3,0) + X.at(0,2)*X.at(2,0)*X.at(3,1) - X.at(0,0)*X.at(2,2)*X.at(3,1) - X.at(0,1)*X.at(2,0)*X.at(3,2) + X.at(0,0)*X.at(2,1)*X.at(3,2) ) / tmp_det;
        
        out.at(0,2) = ( X.at(0,2)*X.at(1,3)*X.at(3,1) - X.at(0,3)*X.at(1,2)*X.at(3,1) + X.at(0,3)*X.at(1,1)*X.at(3,2) - X.at(0,1)*X.at(1,3)*X.at(3,2) - X.at(0,2)*X.at(1,1)*X.at(3,3) + X.at(0,1)*X.at(1,2)*X.at(3,3) ) / tmp_det;
        out.at(1,2) = ( X.at(0,3)*X.at(1,2)*X.at(3,0) - X.at(0,2)*X.at(1,3)*X.at(3,0) - X.at(0,3)*X.at(1,0)*X.at(3,2) + X.at(0,0)*X.at(1,3)*X.at(3,2) + X.at(0,2)*X.at(1,0)*X.at(3,3) - X.at(0,0)*X.at(1,2)*X.at(3,3) ) / tmp_det;
        out.at(2,2) = ( X.at(0,1)*X.at(1,3)*X.at(3,0) - X.at(0,3)*X.at(1,1)*X.at(3,0) + X.at(0,3)*X.at(1,0)*X.at(3,1) - X.at(0,0)*X.at(1,3)*X.at(3,1) - X.at(0,1)*X.at(1,0)*X.at(3,3) + X.at(0,0)*X.at(1,1)*X.at(3,3) ) / tmp_det;
        out.at(3,2) = ( X.at(0,2)*X.at(1,1)*X.at(3,0) - X.at(0,1)*X.at(1,2)*X.at(3,0) - X.at(0,2)*X.at(1,0)*X.at(3,1) + X.at(0,0)*X.at(1,2)*X.at(3,1) + X.at(0,1)*X.at(1,0)*X.at(3,2) - X.at(0,0)*X.at(1,1)*X.at(3,2) ) / tmp_det;
        
        out.at(0,3) = ( X.at(0,3)*X.at(1,2)*X.at(2,1) - X.at(0,2)*X.at(1,3)*X.at(2,1) - X.at(0,3)*X.at(1,1)*X.at(2,2) + X.at(0,1)*X.at(1,3)*X.at(2,2) + X.at(0,2)*X.at(1,1)*X.at(2,3) - X.at(0,1)*X.at(1,2)*X.at(2,3) ) / tmp_det;
        out.at(1,3) = ( X.at(0,2)*X.at(1,3)*X.at(2,0) - X.at(0,3)*X.at(1,2)*X.at(2,0) + X.at(0,3)*X.at(1,0)*X.at(2,2) - X.at(0,0)*X.at(1,3)*X.at(2,2) - X.at(0,2)*X.at(1,0)*X.at(2,3) + X.at(0,0)*X.at(1,2)*X.at(2,3) ) / tmp_det;
        out.at(2,3) = ( X.at(0,3)*X.at(1,1)*X.at(2,0) - X.at(0,1)*X.at(1,3)*X.at(2,0) - X.at(0,3)*X.at(1,0)*X.at(2,1) + X.at(0,0)*X.at(1,3)*X.at(2,1) + X.at(0,1)*X.at(1,0)*X.at(2,3) - X.at(0,0)*X.at(1,1)*X.at(2,3) ) / tmp_det;
        out.at(3,3) = ( X.at(0,1)*X.at(1,2)*X.at(2,0) - X.at(0,2)*X.at(1,1)*X.at(2,0) + X.at(0,2)*X.at(1,0)*X.at(2,1) - X.at(0,0)*X.at(1,2)*X.at(2,1) - X.at(0,1)*X.at(1,0)*X.at(2,2) + X.at(0,0)*X.at(1,1)*X.at(2,2) ) / tmp_det;
        }
      else
        {
        det_ok = false;
        }
      };
      break;
    
    default:
      ;
    }
  
  return det_ok;
  }



template<typename eT>
inline
bool
auxlib::inv_inplace_tinymat(Mat<eT>& X, const u32 N)
  {
  arma_extra_debug_sigprint();
  
  bool det_ok = true;
  
  // for more info, see:
  // http://www.dr-lex.34sp.com/random/matrix_inv.html
  // http://www.cvl.iis.u-tokyo.ac.jp/~miyazaki/tech/teche23.html
  // http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
  // http://www.geometrictools.com//LibFoundation/Mathematics/Wm4Matrix4.inl
  
  switch(N)
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
      
      const eT tmp_det = (a*d - b*c);
      
      if(tmp_det != eT(0))
        {
        X.at(0,0) =  d / tmp_det;
        X.at(0,1) = -b / tmp_det;
        X.at(1,0) = -c / tmp_det;
        X.at(1,1) =  a / tmp_det;
        }
      else
        {
        det_ok = false;
        }
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
      
      const eT tmp_det = a11*(a33*a22 - a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12 - a22*a13);
      
      if(tmp_det != eT(0))
        {
        X_col0[0] =  (a33*a22 - a32*a23) / tmp_det;
        X_col0[1] = -(a33*a21 - a31*a23) / tmp_det;
        X_col0[2] =  (a32*a21 - a31*a22) / tmp_det;
        
        X_col1[0] = -(a33*a12 - a32*a13) / tmp_det;
        X_col1[1] =  (a33*a11 - a31*a13) / tmp_det;
        X_col1[2] = -(a32*a11 - a31*a12) / tmp_det;
        
        X_col2[0] =  (a23*a12 - a22*a13) / tmp_det;
        X_col2[1] = -(a23*a11 - a21*a13) / tmp_det;
        X_col2[2] =  (a22*a11 - a21*a12) / tmp_det;
        }
      else
        {
        det_ok = false;
        }
      };
      break;
      
    case 4:
      {
      const eT tmp_det = det(X);
      
      if(tmp_det != eT(0))
        {
        const Mat<eT> A(X);
        
        X.at(0,0) = ( A.at(1,2)*A.at(2,3)*A.at(3,1) - A.at(1,3)*A.at(2,2)*A.at(3,1) + A.at(1,3)*A.at(2,1)*A.at(3,2) - A.at(1,1)*A.at(2,3)*A.at(3,2) - A.at(1,2)*A.at(2,1)*A.at(3,3) + A.at(1,1)*A.at(2,2)*A.at(3,3) ) / tmp_det;
        X.at(1,0) = ( A.at(1,3)*A.at(2,2)*A.at(3,0) - A.at(1,2)*A.at(2,3)*A.at(3,0) - A.at(1,3)*A.at(2,0)*A.at(3,2) + A.at(1,0)*A.at(2,3)*A.at(3,2) + A.at(1,2)*A.at(2,0)*A.at(3,3) - A.at(1,0)*A.at(2,2)*A.at(3,3) ) / tmp_det;
        X.at(2,0) = ( A.at(1,1)*A.at(2,3)*A.at(3,0) - A.at(1,3)*A.at(2,1)*A.at(3,0) + A.at(1,3)*A.at(2,0)*A.at(3,1) - A.at(1,0)*A.at(2,3)*A.at(3,1) - A.at(1,1)*A.at(2,0)*A.at(3,3) + A.at(1,0)*A.at(2,1)*A.at(3,3) ) / tmp_det;
        X.at(3,0) = ( A.at(1,2)*A.at(2,1)*A.at(3,0) - A.at(1,1)*A.at(2,2)*A.at(3,0) - A.at(1,2)*A.at(2,0)*A.at(3,1) + A.at(1,0)*A.at(2,2)*A.at(3,1) + A.at(1,1)*A.at(2,0)*A.at(3,2) - A.at(1,0)*A.at(2,1)*A.at(3,2) ) / tmp_det;
        
        X.at(0,1) = ( A.at(0,3)*A.at(2,2)*A.at(3,1) - A.at(0,2)*A.at(2,3)*A.at(3,1) - A.at(0,3)*A.at(2,1)*A.at(3,2) + A.at(0,1)*A.at(2,3)*A.at(3,2) + A.at(0,2)*A.at(2,1)*A.at(3,3) - A.at(0,1)*A.at(2,2)*A.at(3,3) ) / tmp_det;
        X.at(1,1) = ( A.at(0,2)*A.at(2,3)*A.at(3,0) - A.at(0,3)*A.at(2,2)*A.at(3,0) + A.at(0,3)*A.at(2,0)*A.at(3,2) - A.at(0,0)*A.at(2,3)*A.at(3,2) - A.at(0,2)*A.at(2,0)*A.at(3,3) + A.at(0,0)*A.at(2,2)*A.at(3,3) ) / tmp_det;
        X.at(2,1) = ( A.at(0,3)*A.at(2,1)*A.at(3,0) - A.at(0,1)*A.at(2,3)*A.at(3,0) - A.at(0,3)*A.at(2,0)*A.at(3,1) + A.at(0,0)*A.at(2,3)*A.at(3,1) + A.at(0,1)*A.at(2,0)*A.at(3,3) - A.at(0,0)*A.at(2,1)*A.at(3,3) ) / tmp_det;
        X.at(3,1) = ( A.at(0,1)*A.at(2,2)*A.at(3,0) - A.at(0,2)*A.at(2,1)*A.at(3,0) + A.at(0,2)*A.at(2,0)*A.at(3,1) - A.at(0,0)*A.at(2,2)*A.at(3,1) - A.at(0,1)*A.at(2,0)*A.at(3,2) + A.at(0,0)*A.at(2,1)*A.at(3,2) ) / tmp_det;
        
        X.at(0,2) = ( A.at(0,2)*A.at(1,3)*A.at(3,1) - A.at(0,3)*A.at(1,2)*A.at(3,1) + A.at(0,3)*A.at(1,1)*A.at(3,2) - A.at(0,1)*A.at(1,3)*A.at(3,2) - A.at(0,2)*A.at(1,1)*A.at(3,3) + A.at(0,1)*A.at(1,2)*A.at(3,3) ) / tmp_det;
        X.at(1,2) = ( A.at(0,3)*A.at(1,2)*A.at(3,0) - A.at(0,2)*A.at(1,3)*A.at(3,0) - A.at(0,3)*A.at(1,0)*A.at(3,2) + A.at(0,0)*A.at(1,3)*A.at(3,2) + A.at(0,2)*A.at(1,0)*A.at(3,3) - A.at(0,0)*A.at(1,2)*A.at(3,3) ) / tmp_det;
        X.at(2,2) = ( A.at(0,1)*A.at(1,3)*A.at(3,0) - A.at(0,3)*A.at(1,1)*A.at(3,0) + A.at(0,3)*A.at(1,0)*A.at(3,1) - A.at(0,0)*A.at(1,3)*A.at(3,1) - A.at(0,1)*A.at(1,0)*A.at(3,3) + A.at(0,0)*A.at(1,1)*A.at(3,3) ) / tmp_det;
        X.at(3,2) = ( A.at(0,2)*A.at(1,1)*A.at(3,0) - A.at(0,1)*A.at(1,2)*A.at(3,0) - A.at(0,2)*A.at(1,0)*A.at(3,1) + A.at(0,0)*A.at(1,2)*A.at(3,1) + A.at(0,1)*A.at(1,0)*A.at(3,2) - A.at(0,0)*A.at(1,1)*A.at(3,2) ) / tmp_det;
        
        X.at(0,3) = ( A.at(0,3)*A.at(1,2)*A.at(2,1) - A.at(0,2)*A.at(1,3)*A.at(2,1) - A.at(0,3)*A.at(1,1)*A.at(2,2) + A.at(0,1)*A.at(1,3)*A.at(2,2) + A.at(0,2)*A.at(1,1)*A.at(2,3) - A.at(0,1)*A.at(1,2)*A.at(2,3) ) / tmp_det;
        X.at(1,3) = ( A.at(0,2)*A.at(1,3)*A.at(2,0) - A.at(0,3)*A.at(1,2)*A.at(2,0) + A.at(0,3)*A.at(1,0)*A.at(2,2) - A.at(0,0)*A.at(1,3)*A.at(2,2) - A.at(0,2)*A.at(1,0)*A.at(2,3) + A.at(0,0)*A.at(1,2)*A.at(2,3) ) / tmp_det;
        X.at(2,3) = ( A.at(0,3)*A.at(1,1)*A.at(2,0) - A.at(0,1)*A.at(1,3)*A.at(2,0) - A.at(0,3)*A.at(1,0)*A.at(2,1) + A.at(0,0)*A.at(1,3)*A.at(2,1) + A.at(0,1)*A.at(1,0)*A.at(2,3) - A.at(0,0)*A.at(1,1)*A.at(2,3) ) / tmp_det;
        X.at(3,3) = ( A.at(0,1)*A.at(1,2)*A.at(2,0) - A.at(0,2)*A.at(1,1)*A.at(2,0) + A.at(0,2)*A.at(1,0)*A.at(2,1) - A.at(0,0)*A.at(1,2)*A.at(2,1) - A.at(0,1)*A.at(1,0)*A.at(2,2) + A.at(0,0)*A.at(1,1)*A.at(2,2) ) / tmp_det;
        }
      else
        {
        det_ok = false;
        }
      };
      break;
      
    default:
      ;
    }
  
  return det_ok;
  }



template<typename eT>
inline
bool
auxlib::inv_inplace_lapack(Mat<eT>& out)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_ATLAS)
    {
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
    blas_int n_rows = out.n_rows;
    blas_int n_cols = out.n_cols;
    blas_int info   = 0;
    
    podarray<blas_int> ipiv(out.n_rows);
    
    // 84 was empirically found -- it is the maximum value suggested by LAPACK (as provided by ATLAS v3.6)
    // based on tests with various matrix types on 32-bit and 64-bit machines
    //
    // the "work" array is deliberately long so that a secondary (time-consuming)
    // memory allocation is avoided, if possible
    
    blas_int work_len = (std::max)(blas_int(1), n_rows*84);
    podarray<eT> work(work_len);
    
    lapack::getrf(&n_rows, &n_cols, out.memptr(), &n_rows, ipiv.memptr(), &info);
    
    if(info == 0)
      {
      // query for optimum size of work_len
      
      blas_int work_len_tmp = -1;
      lapack::getri(&n_rows, out.memptr(), &n_rows, ipiv.memptr(), work.memptr(), &work_len_tmp, &info);
      
      if(info == 0)
        {
        blas_int proposed_work_len = static_cast<blas_int>(access::tmp_real(work[0]));
        
        // if necessary, allocate more memory
        if(work_len < proposed_work_len)
          {
          work_len = proposed_work_len;
          work.set_size(work_len);
          }
        }
      
      lapack::getri(&n_rows, out.memptr(), &n_rows, ipiv.memptr(), work.memptr(), &work_len, &info);
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_stop("inv(): need ATLAS or LAPACK");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
eT
auxlib::det(const Base<eT,T1>& X)
  {
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( !A.is_square(), "det(): matrix is not square" );
  
  const bool make_copy = (is_Mat<T1>::value == true) ? true : false;
  
  const u32 N = A.n_rows;
  
  switch(N)
    {
    case 0:
    case 1:
    case 2:
      return auxlib::det_tinymat(A, N);
      break;
    
    case 3:
    case 4:
      {
      const eT tmp_det = auxlib::det_tinymat(A, N);
      return (tmp_det != eT(0)) ? tmp_det : auxlib::det_lapack(A, make_copy);
      }
      break;
    
    default:
      return auxlib::det_lapack(A, make_copy);
    }
  
  return eT(0);  // prevent compiler warnings
  }



template<typename eT>
inline
eT
auxlib::det_tinymat(const Mat<eT>& X, const u32 N)
  {
  arma_extra_debug_sigprint();
  
  switch(N)
    {
    case 0:
      return eT(0);
    
    case 1:
      return X[0];
    
    case 2:
      return (X.at(0,0)*X.at(1,1) - X.at(0,1)*X.at(1,0));
    
    case 3:
      {
      // const double tmp1 = X.at(0,0) * X.at(1,1) * X.at(2,2);
      // const double tmp2 = X.at(0,1) * X.at(1,2) * X.at(2,0);
      // const double tmp3 = X.at(0,2) * X.at(1,0) * X.at(2,1);
      // const double tmp4 = X.at(2,0) * X.at(1,1) * X.at(0,2);
      // const double tmp5 = X.at(2,1) * X.at(1,2) * X.at(0,0);
      // const double tmp6 = X.at(2,2) * X.at(1,0) * X.at(0,1);
      // return (tmp1+tmp2+tmp3) - (tmp4+tmp5+tmp6);
      
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
      return eT(0);
      ;
    }
  }



//! immediate determinant of a matrix using ATLAS or LAPACK
template<typename eT>
inline
eT
auxlib::det_lapack(const Mat<eT>& X, const bool make_copy)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> X_copy;
  
  if(make_copy == true)
    {
    X_copy = X;
    }
  
  Mat<eT>& tmp = (make_copy == true) ? X_copy : const_cast< Mat<eT>& >(X);
  
  #if defined(ARMA_USE_ATLAS)
    {
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
    
    return ( (sign < 0) ? -val : val );
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    podarray<blas_int> ipiv(tmp.n_rows);
    
    blas_int info   = 0;
    blas_int n_rows = blas_int(tmp.n_rows);
    blas_int n_cols = blas_int(tmp.n_cols);
    
    lapack::getrf(&n_rows, &n_cols, tmp.memptr(), &n_rows, ipiv.memptr(), &info);
    
    // on output tmp appears to be L+U_alt, where U_alt is U with the main diagonal set to zero
    eT val = tmp.at(0,0);
    for(u32 i=1; i < tmp.n_rows; ++i)
      {
      val *= tmp.at(i,i);
      }
    
    blas_int sign = +1;
    for(u32 i=0; i < tmp.n_rows; ++i)
      {
      if( blas_int(i) != (ipiv.mem[i] - 1) )  // NOTE: adjustment of -1 is required as Fortran counts from 1
        {
        sign *= -1;
        }
      }
    
    return ( (sign < 0) ? -val : val );
    }
  #else
    {
    arma_ignore(X);
    arma_ignore(make_copy);
    arma_ignore(tmp);
    arma_stop("det(): need ATLAS or LAPACK");
    return eT(0);
    }
  #endif
  }



//! immediate log determinant of a matrix using ATLAS or LAPACK
template<typename eT, typename T1>
inline
void
auxlib::log_det(eT& out_val, typename get_pod_type<eT>::result& out_sign, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  #if defined(ARMA_USE_ATLAS)
    {
    Mat<eT> tmp(X.get_ref());
    arma_debug_check( (tmp.is_square() == false), "log_det(): given matrix is not square" );
    
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
    Mat<eT> tmp(X.get_ref());
    arma_debug_check( (tmp.is_square() == false), "log_det(): given matrix is not square" );
    
    podarray<blas_int> ipiv(tmp.n_rows);
    
    blas_int info   = 0;
    blas_int n_rows = blas_int(tmp.n_rows);
    blas_int n_cols = blas_int(tmp.n_cols);
    
    lapack::getrf(&n_rows, &n_cols, tmp.memptr(), &n_rows, ipiv.memptr(), &info);
    
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
      if( blas_int(i) != (ipiv.mem[i] - 1) )  // NOTE: adjustment of -1 is required as Fortran counts from 1
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
template<typename eT, typename T1>
inline
void
auxlib::lu(Mat<eT>& L, Mat<eT>& U, podarray<blas_int>& ipiv, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  U = X.get_ref();
  
  if(U.is_empty())
    {
    ipiv.reset();
    L.reset();
    return;
    }
  
  const u32 U_n_rows = U.n_rows;
  const u32 U_n_cols = U.n_cols;
  
  #if defined(ARMA_USE_ATLAS) || defined(ARMA_USE_LAPACK)
    {
    #if defined(ARMA_USE_ATLAS)
      {
      ipiv.set_size( (std::min)(U_n_rows, U_n_cols) );
      
      //int info = 
      atlas::clapack_getrf(atlas::CblasColMajor, U_n_rows, U_n_cols, U.memptr(), U_n_rows, ipiv.memptr());
      }
    #elif defined(ARMA_USE_LAPACK)
      {
      ipiv.set_size( (std::min)(U_n_rows, U_n_cols) );
      
      blas_int info = 0;
      
      blas_int n_rows = U_n_rows;
      blas_int n_cols = U_n_cols;
      
      
      lapack::getrf(&n_rows, &n_cols, U.memptr(), &n_rows, ipiv.memptr(), &info);
      
      // take into account that Fortran counts from 1
      arrayops::inplace_minus(ipiv.memptr(), blas_int(1), ipiv.n_elem);
      }
    #endif
    
    L.copy_size(U);
    
    for(u32 col=0; col < U_n_cols; ++col)
      {
      for(u32 row=0; (row < col) && (row < U_n_rows); ++row)
        {
        L.at(row,col) = eT(0);
        }
      
      if( L.in_range(col,col) == true )
        {
        L.at(col,col) = eT(1);
        }
      
      for(u32 row = (col+1); row < U_n_rows; ++row)
        {
        L.at(row,col) = U.at(row,col);
        U.at(row,col) = eT(0);
        }
      }
    }
  #else
    {
    arma_ignore(L);
    arma_ignore(ipiv);
    arma_ignore(U_n_rows);
    arma_ignore(U_n_cols);
    arma_stop("lu(): need ATLAS or LAPACK");
    }
  #endif
  }



template<typename eT, typename T1>
inline
void
auxlib::lu(Mat<eT>& L, Mat<eT>& U, Mat<eT>& P, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  podarray<blas_int> ipiv1;
  auxlib::lu(L, U, ipiv1, X);
  
  const u32 n      = ipiv1.n_elem;
  const u32 P_rows = U.n_rows;
  
  podarray<blas_int> ipiv2(P_rows);
  
  const blas_int* ipiv1_mem = ipiv1.memptr();
        blas_int* ipiv2_mem = ipiv2.memptr();
  
  for(u32 i=0; i<P_rows; ++i)
    {
    ipiv2_mem[i] = i;
    }
  
  for(u32 i=0; i<n; ++i)
    {
    const u32 k = ipiv1_mem[i];
    
    if( ipiv2_mem[i] != ipiv2_mem[k] )
      {
      std::swap( ipiv2_mem[i], ipiv2_mem[k] );
      }
    }
  
  P.zeros(P_rows, P_rows);
  
  for(u32 row=0; row<P_rows; ++row)
    {
    P.at(row, ipiv2_mem[row]) = eT(1);
    }
  
  if(L.n_cols > U.n_rows)
    {
    L.shed_cols(U.n_rows, L.n_cols-1);
    }
    
  if(U.n_rows > L.n_cols)
    {
    U.shed_rows(L.n_cols, U.n_rows-1);
    }
  }



template<typename eT, typename T1>
inline
void
auxlib::lu(Mat<eT>& L, Mat<eT>& U, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  podarray<blas_int> ipiv1;
  auxlib::lu(L, U, ipiv1, X);
  
  const u32 n      = ipiv1.n_elem;
  const u32 P_rows = U.n_rows;
  
  podarray<blas_int> ipiv2(P_rows);
  
  const blas_int* ipiv1_mem = ipiv1.memptr();
        blas_int* ipiv2_mem = ipiv2.memptr();
  
  for(u32 i=0; i<P_rows; ++i)
    {
    ipiv2_mem[i] = i;
    }
  
  for(u32 i=0; i<n; ++i)
    {
    const u32 k = ipiv1_mem[i];
    
    if( ipiv2_mem[i] != ipiv2_mem[k] )
      {
      std::swap( ipiv2_mem[i], ipiv2_mem[k] );
      L.swap_rows( u32(ipiv2_mem[i]), u32(ipiv2_mem[k]) );
      }
    }
  
  if(L.n_cols > U.n_rows)
    {
    L.shed_cols(U.n_rows, L.n_cols-1);
    }
    
  if(U.n_rows > L.n_cols)
    {
    U.shed_rows(L.n_cols, U.n_rows-1);
    }
  }



//! immediate eigenvalues of a symmetric real matrix using LAPACK
template<typename eT, typename T1>
inline
bool
auxlib::eig_sym(Col<eT>& eigval, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A(X.get_ref());
    
    arma_debug_check( (A.is_square() == false), "eig_sym(): given matrix is not square");
    
    // rudimentary "better-than-nothing" test for symmetry
    //arma_debug_check( (A.at(A.n_rows-1, 0) != A.at(0, A.n_cols-1)), "auxlib::eig(): given matrix is not symmetric" );
    
    char jobz  = 'N';
    char uplo  = 'U';
    
    blas_int n_rows = A.n_rows;
    blas_int lwork  = (std::max)(blas_int(1), 3*n_rows-1);
    
    eigval.set_size(n_rows);
    podarray<eT> work(lwork);
    
    blas_int info;
    
    arma_extra_debug_print("lapack::syev()");
    lapack::syev(&jobz, &uplo, &n_rows, A.memptr(), &n_rows, eigval.memptr(), work.memptr(), &lwork, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(X);
    arma_stop("eig_sym(): need LAPACK");
    return false;
    }
  #endif
  }



//! immediate eigenvalues of a hermitian complex matrix using LAPACK
template<typename T, typename T1>
inline
bool
auxlib::eig_sym(Col<T>& eigval, const Base<std::complex<T>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A(X.get_ref());
    arma_debug_check( (A.is_square() == false), "eig_sym(): given matrix is not hermitian");
    
    char jobz  = 'N'; 
    char uplo  = 'U';
    
    blas_int n_rows = A.n_rows;
    blas_int lda    = A.n_rows;
    blas_int lwork  = (std::max)(blas_int(1), 2*n_rows - 1);  // TODO: automatically find best size of lwork
    
    eigval.set_size(n_rows);
    
    podarray<eT> work(lwork);
    podarray<T>  rwork( (std::max)(blas_int(1), 3*n_rows - 2) );
    
    blas_int info;
    
    arma_extra_debug_print("lapack::heev()");
    lapack::heev(&jobz, &uplo, &n_rows, A.memptr(), &lda, eigval.memptr(), work.memptr(), &lwork, rwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(X);
    arma_stop("eig_sym(): need LAPACK");
    return false;
    }
  #endif
  }



//! immediate eigenvalues and eigenvectors of a symmetric real matrix using LAPACK
template<typename eT, typename T1>
inline
bool
auxlib::eig_sym(Col<eT>& eigval, Mat<eT>& eigvec, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    eigvec = X.get_ref();
    
    arma_debug_check( (eigvec.is_square() == false), "eig_sym(): given matrix is not square" );
    
    // rudimentary "better-than-nothing" test for symmetry
    //arma_debug_check( (A.at(A.n_rows-1, 0) != A.at(0, A.n_cols-1)), "auxlib::eig(): given matrix is not symmetric" );
    
    char jobz  = 'V';
    char uplo  = 'U';
    
    blas_int n_rows = eigvec.n_rows;
    blas_int lwork  = (std::max)(blas_int(1), 3*n_rows-1);
    
    eigval.set_size(n_rows);
    podarray<eT> work(lwork);
    
    blas_int info;
    
    arma_extra_debug_print("lapack::syev()");
    lapack::syev(&jobz, &uplo, &n_rows, eigvec.memptr(), &n_rows, eigval.memptr(), work.memptr(), &lwork, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_stop("eig_sym(): need LAPACK");
    return false;
    }
  #endif
  }



//! immediate eigenvalues and eigenvectors of a hermitian complex matrix using LAPACK
template<typename T, typename T1>
inline
bool
auxlib::eig_sym(Col<T>& eigval, Mat< std::complex<T> >& eigvec, const Base<std::complex<T>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  #if defined(ARMA_USE_LAPACK)
    {
    eigvec = X.get_ref();
    
    arma_debug_check( (eigvec.is_square() == false), "eig_sym(): given matrix is not hermitian" );
    
    char jobz  = 'V';
    char uplo  = 'U';
    
    blas_int n_rows = eigvec.n_rows;
    blas_int lda    = eigvec.n_rows;
    blas_int lwork  = (std::max)(blas_int(1), 2*n_rows - 1);  // TODO: automatically find best size of lwork
    
    eigval.set_size(n_rows);
    
    podarray<eT> work(lwork);
    podarray<T>  rwork( (std::max)(blas_int(1), 3*n_rows - 2) );
    
    blas_int info;
    
    arma_extra_debug_print("lapack::heev()");
    lapack::heev(&jobz, &uplo, &n_rows, eigvec.memptr(), &lda, eigval.memptr(), work.memptr(), &lwork, rwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_stop("eig_sym(): need LAPACK");
    return false;
    }
  #endif
  }



//! Eigenvalues and eigenvectors of a general square real matrix using LAPACK.
//! The argument 'side' specifies which eigenvectors should be calculated
//! (see code for mode details).
template<typename T, typename T1>
inline
bool
auxlib::eig_gen
  (
  Col< std::complex<T> >& eigval,
  Mat<T>& l_eigvec,
  Mat<T>& r_eigvec,
  const Base<T,T1>& X,
  const char side
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<T> A(X.get_ref());
    arma_debug_check( (A.is_square() == false), "eig_gen(): given matrix is not square" );
    
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
    
       
    blas_int n_rows = A.n_rows;
    blas_int lda    = A.n_rows;
    blas_int lwork  = (std::max)(blas_int(1), 4*n_rows);  // TODO: automatically find best size of lwork
    
    eigval.set_size(n_rows);
    l_eigvec.set_size(n_rows, n_rows);
    r_eigvec.set_size(n_rows, n_rows);
    
    podarray<T> work(lwork);
    podarray<T> rwork( (std::max)(blas_int(1), 3*n_rows) );
    
    podarray<T> wr(n_rows);
    podarray<T> wi(n_rows);
    
    Mat<T> A_copy = A;
    blas_int info;
    
    arma_extra_debug_print("lapack::geev()");
    lapack::geev(&jobvl, &jobvr, &n_rows, A_copy.memptr(), &lda, wr.memptr(), wi.memptr(), l_eigvec.memptr(), &n_rows, r_eigvec.memptr(), &n_rows, work.memptr(), &lwork, &info);
    
    
    eigval.set_size(n_rows);
    for(u32 i=0; i<u32(n_rows); ++i)
      {
      eigval[i] = std::complex<T>(wr[i], wi[i]);
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(l_eigvec);
    arma_ignore(r_eigvec);
    arma_ignore(X);
    arma_ignore(side);
    arma_stop("eig_gen(): need LAPACK");
    return false;
    }
  #endif
  }





//! Eigenvalues and eigenvectors of a general square complex matrix using LAPACK
//! The argument 'side' specifies which eigenvectors should be calculated
//! (see code for mode details).
template<typename T, typename T1>
inline
bool
auxlib::eig_gen
  (
  Col< std::complex<T> >& eigval,
  Mat< std::complex<T> >& l_eigvec, 
  Mat< std::complex<T> >& r_eigvec, 
  const Base< std::complex<T>, T1 >& X, 
  const char side
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A(X.get_ref());
    arma_debug_check( (A.is_square() == false), "eig_gen(): given matrix is not square" );
    
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
    
       
    blas_int n_rows = A.n_rows;
    blas_int lda    = A.n_rows;
    blas_int lwork  = (std::max)(blas_int(1), 4*n_rows);  // TODO: automatically find best size of lwork
    
    eigval.set_size(n_rows);
    l_eigvec.set_size(n_rows, n_rows);
    r_eigvec.set_size(n_rows, n_rows);
    
    podarray<eT> work(lwork);
    podarray<T>  rwork( (std::max)(blas_int(1), 3*n_rows) );  // was 2,3
    
    blas_int info;
    
    arma_extra_debug_print("lapack::cx_geev()");
    lapack::cx_geev(&jobvl, &jobvr, &n_rows, A.memptr(), &lda, eigval.memptr(), l_eigvec.memptr(), &n_rows, r_eigvec.memptr(), &n_rows, work.memptr(), &lwork, rwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(l_eigvec);
    arma_ignore(r_eigvec);
    arma_ignore(X);
    arma_ignore(side);
    arma_stop("eig_gen(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
bool
auxlib::chol(Mat<eT>& out, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    out = X.get_ref();
    
    arma_debug_check( (out.is_square() == false), "chol(): given matrix is not square" );
    
    const u32 out_n_rows = out.n_rows;
    
    char      uplo = 'U';
    blas_int  n    = out_n_rows;
    blas_int  info;
    
    lapack::potrf(&uplo, &n, out.memptr(), &n, &info);
    
    for(u32 col=0; col<out_n_rows; ++col)
      {
      eT* colptr = out.colptr(col);
      
      for(u32 row=(col+1); row < out_n_rows; ++row)
        {
        colptr[row] = eT(0);
        }
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(X);
    arma_stop("chol(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
bool 
auxlib::qr(Mat<eT>& Q, Mat<eT>& R, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    R = X.get_ref();
    
    if(R.is_empty())
      {
      Q.reset();
      return false;
      }
    
    const u32 R_n_rows = R.n_rows;
    const u32 R_n_cols = R.n_cols;
    
    blas_int m            = static_cast<blas_int>(R_n_rows);
    blas_int n            = static_cast<blas_int>(R_n_cols);
    blas_int work_len     = (std::max)(blas_int(1),n);
    blas_int work_len_tmp;
    blas_int k            = (std::min)(m,n);
    blas_int info;
    
    podarray<eT> tau(k);
    podarray<eT> work(work_len);
    
    // query for the optimum value of work_len
    work_len_tmp = -1;
    lapack::geqrf(&m, &n, R.memptr(), &m, tau.memptr(), work.memptr(), &work_len_tmp, &info);
    
    if(info == 0)
      {
      work_len = static_cast<blas_int>(access::tmp_real(work[0]));
      work.set_size(work_len);
      }
    
    lapack::geqrf(&m, &n, R.memptr(), &m, tau.memptr(), work.memptr(), &work_len, &info);
    
    Q.set_size(R_n_rows, R_n_rows);
    
    syslib::copy_elem( Q.memptr(), R.memptr(), (std::min)(Q.n_elem, R.n_elem) );
    
    //
    // construct R
    
    for(u32 col=0; col < R_n_cols; ++col)
      {
      for(u32 row=(col+1); row < R_n_rows; ++row)
        {
        R.at(row,col) = eT(0);
        }
      }
    
    
    if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
      {
      // query for the optimum value of work_len
      work_len_tmp = -1;
      lapack::orgqr(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &work_len_tmp, &info);
      
      if(info == 0)
        {
        work_len = static_cast<blas_int>(access::tmp_real(work[0]));
        work.set_size(work_len);
        }
      
      lapack::orgqr(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &work_len, &info);
      }
    else
    if( (is_supported_complex_float<eT>::value == true) || (is_supported_complex_double<eT>::value == true) )
      {
      // query for the optimum value of work_len
      work_len_tmp = -1;
      lapack::ungqr(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &work_len_tmp, &info);
      
      if(info == 0)
        {
        work_len = static_cast<blas_int>(access::tmp_real(work[0]));
        work.set_size(work_len);
        }
      
      lapack::ungqr(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &work_len, &info);
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(Q);
    arma_ignore(R);
    arma_ignore(X);
    arma_stop("qr(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
bool
auxlib::svd(Col<eT>& S, const Base<eT,T1>& X, u32& X_n_rows, u32& X_n_cols)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A(X.get_ref());
    
    X_n_rows = A.n_rows;
    X_n_cols = A.n_cols;
    
    if(A.is_empty())
      {
      S.reset();
      return false;
      }
    
    Mat<eT> U(1, 1);
    Mat<eT> V(1, A.n_cols);
    
    char jobu  = 'N';
    char jobvt = 'N';
    
    blas_int  m     = A.n_rows;
    blas_int  n     = A.n_cols;
    blas_int  lda   = A.n_rows;
    blas_int  ldu   = U.n_rows;
    blas_int  ldvt  = V.n_rows;
    blas_int  lwork = 2 * (std::max)(blas_int(1), (std::max)( (3*(std::min)(m,n) + (std::max)(m,n)), 5*(std::min)(m,n) ) );
    blas_int  info;
    
    S.set_size( (std::min)(m, n) );
    
    podarray<eT> work(lwork);
  
  
    // let gesvd_() calculate the optimum size of the workspace
    blas_int lwork_tmp = -1;
    
    lapack::gesvd<eT>
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
      blas_int proposed_lwork = static_cast<blas_int>(work[0]);
      
      if(proposed_lwork > lwork)
        {
        lwork = proposed_lwork;
        work.set_size(lwork);
        }
      
      lapack::gesvd<eT>
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
    arma_ignore(S);
    arma_ignore(X);
    arma_ignore(X_n_rows);
    arma_ignore(X_n_cols);
    arma_stop("svd(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename T, typename T1>
inline
bool
auxlib::svd(Col<T>& S, const Base<std::complex<T>, T1>& X, u32& X_n_rows, u32& X_n_cols)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A(X.get_ref());
    
    X_n_rows = A.n_rows;
    X_n_cols = A.n_cols;
    
    if(A.is_empty())
      {
      S.reset();
      return false;
      }
    
    Mat<eT> U(1, 1);
    Mat<eT> V(1, A.n_cols);
    
    char jobu  = 'N';
    char jobvt = 'N';
    
    blas_int  m     = A.n_rows;
    blas_int  n     = A.n_cols;
    blas_int  lda   = A.n_rows;
    blas_int  ldu   = U.n_rows;
    blas_int  ldvt  = V.n_rows;
    blas_int  lwork = 2 * (std::max)(blas_int(1), 2*(std::min)(m,n)+(std::max)(m,n) );
    blas_int  info;
    
    S.set_size( (std::min)(m,n) );
    
    podarray<eT> work(lwork);
    podarray<T>  rwork( 5*(std::min)(m,n) );
    
    // let gesvd_() calculate the optimum size of the workspace
    blas_int lwork_tmp = -1;
    
    lapack::cx_gesvd<T>
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
      blas_int proposed_lwork = static_cast<blas_int>(real(work[0]));
      if(proposed_lwork > lwork)
        {
        lwork = proposed_lwork;
        work.set_size(lwork);
        }
      
      lapack::cx_gesvd<T>
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
    arma_ignore(S);
    arma_ignore(X);
    arma_ignore(X_n_rows);
    arma_ignore(X_n_cols);

    arma_stop("svd(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
bool
auxlib::svd(Col<eT>& S, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  u32 junk;
  return auxlib::svd(S, X, junk, junk);
  }



template<typename T, typename T1>
inline
bool
auxlib::svd(Col<T>& S, const Base<std::complex<T>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  u32 junk;
  return auxlib::svd(S, X, junk, junk);
  }



template<typename eT, typename T1>
inline
bool
auxlib::svd(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A(X.get_ref());
    
    if(A.is_empty())
      {
      U.reset();
      S.reset();
      V.reset();
      return false;
      }
    
    U.set_size(A.n_rows, A.n_rows);
    V.set_size(A.n_cols, A.n_cols);
    
    char jobu  = 'A';
    char jobvt = 'A';
    
    blas_int  m     = A.n_rows;
    blas_int  n     = A.n_cols;
    blas_int  lda   = A.n_rows;
    blas_int  ldu   = U.n_rows;
    blas_int  ldvt  = V.n_rows;
    blas_int  lwork = 2 * (std::max)(blas_int(1), (std::max)( (3*(std::min)(m,n) + (std::max)(m,n)), 5*(std::min)(m,n) ) );
    blas_int  info;
    
    
    S.set_size( (std::min)(m,n) );
    podarray<eT> work(lwork);
  
    // let gesvd_() calculate the optimum size of the workspace
    blas_int lwork_tmp = -1;
    
    lapack::gesvd<eT>
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
      blas_int proposed_lwork = static_cast<blas_int>(work[0]);
      if(proposed_lwork > lwork)
        {
        lwork = proposed_lwork;
        work.set_size(lwork);
        }
      
      lapack::gesvd<eT>
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
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(X);
    arma_stop("svd(): need LAPACK");
    return false;
    }
  #endif
  }



template<typename T, typename T1>
inline
bool
auxlib::svd(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, const Base< std::complex<T>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  #if defined(ARMA_USE_LAPACK)
    {
    Mat<eT> A(X.get_ref());
    
    if(A.is_empty())
      {
      U.reset();
      S.reset();
      V.reset();
      return false;
      }
    
    U.set_size(A.n_rows, A.n_rows);
    V.set_size(A.n_cols, A.n_cols);
    
    char jobu  = 'A';
    char jobvt = 'A';
    
    blas_int  m     = A.n_rows;
    blas_int  n     = A.n_cols;
    blas_int  lda   = A.n_rows;
    blas_int  ldu   = U.n_rows;
    blas_int  ldvt  = V.n_rows;
    blas_int  lwork = 2 * (std::max)(blas_int(1), 2*(std::min)(m,n)+(std::max)(m,n) );
    blas_int  info;
    
    S.set_size( (std::min)(m,n) );
    
    podarray<eT> work(lwork);
    podarray<T>  rwork( 5*(std::min)(m,n) );
    
    // let gesvd_() calculate the optimum size of the workspace
    blas_int lwork_tmp = -1;
    lapack::cx_gesvd<T>
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
      blas_int proposed_lwork = static_cast<blas_int>(real(work[0]));
      if(proposed_lwork > lwork)
        {
        lwork = proposed_lwork;
        work.set_size(lwork);
        }
      
      lapack::cx_gesvd<T>
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
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(X);
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
  
  if(A.is_empty() || B.is_empty())
    {
    out.reset();
    return false;
    }
  
  #if defined(ARMA_USE_LAPACK)
    {
    blas_int n    = A.n_rows;
    blas_int lda  = A.n_rows;
    blas_int ldb  = A.n_rows;
    blas_int nrhs = B.n_cols;
    blas_int info;
    
    podarray<blas_int> ipiv(n);
    
    out = B;
    Mat<eT> A_copy = A;
  
    lapack::gesv<eT>(&n, &nrhs, A_copy.memptr(), &lda, ipiv.memptr(), out.memptr(), &ldb, &info);
  
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B);
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
  
  if(A.is_empty() || B.is_empty())
    {
    out.reset();
    return false;
    }
  
  #if defined(ARMA_USE_LAPACK)
    {
    char trans = 'N';
    
    blas_int  m     = A.n_rows;
    blas_int  n     = A.n_cols;
    blas_int  lda   = A.n_rows;
    blas_int  ldb   = A.n_rows;
    blas_int  nrhs  = B.n_cols;
    blas_int  lwork = n + (std::max)(n, nrhs);
    blas_int  info;
    
    Mat<eT> A_copy = A;
    Mat<eT> tmp    = B;
    
    
    podarray<eT> work(lwork);
    
    arma_extra_debug_print("lapack::gels()");
    
    // NOTE:
    // the dgels() function in the lapack library supplied by ATLAS 3.6
    // seems to have problems
    
    lapack::gels<eT>
      (
      &trans, &m, &n, &nrhs,
      A_copy.memptr(), &lda,
      tmp.memptr(), &ldb,
      work.memptr(), &lwork,
      &info
      );
    
    arma_extra_debug_print("lapack::gels() -- finished");
    
    out.set_size(A.n_cols, B.n_cols);
    
    for(u32 col=0; col<B.n_cols; ++col)
      {
      syslib::copy_elem( out.colptr(col), tmp.colptr(col), A.n_cols );
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B);
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
  
  if(A.is_empty() || B.is_empty())
    {
    out.reset();
    return false;
    }
  
  #if defined(ARMA_USE_LAPACK)
    {
    char trans = 'N';
    
    blas_int  m     = A.n_rows;
    blas_int  n     = A.n_cols;
    blas_int  lda   = A.n_rows;
    blas_int  ldb   = A.n_cols;
    blas_int  nrhs  = B.n_cols;
    blas_int  lwork = m + (std::max)(m,nrhs);
    blas_int  info;
    
    
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
    
    arma_extra_debug_print("lapack::gels()");
    
    // NOTE:
    // the dgels() function in the lapack library supplied by ATLAS 3.6
    // seems to have problems
    
    lapack::gels<eT>
      (
      &trans, &m, &n, &nrhs,
      A_copy.memptr(), &lda,
      tmp.memptr(), &ldb,
      work.memptr(), &lwork,
      &info
      );
    
    arma_extra_debug_print("lapack::gels() -- finished");
    
    out.set_size(A.n_cols, B.n_cols);
    
    for(u32 col=0; col<B.n_cols; ++col)
      {
      syslib::copy_elem( out.colptr(col), tmp.colptr(col), A.n_cols );
      }
  
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B);
    arma_stop("auxlib::solve_ud(): need LAPACK");
    return false;
    }
  #endif
  }



//! @}
