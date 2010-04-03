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


//! \addtogroup operator_relational
//! @{



template<typename T1, typename T2>
inline
umat
operator==
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;

  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
  
  arma_debug_assert_same_size(A, B, "operator==");
  
  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] == B[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
umat
operator!=
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;

  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
    
  arma_debug_assert_same_size(A, B, "operator!=");
  
  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] != B[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
umat
operator>=
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;

  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
    
  arma_debug_assert_same_size(A, B, "operator>=");
  
  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] >= B[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
umat
operator<=
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;

  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
    
  arma_debug_assert_same_size(A, B, "operator<=");
  
  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] <= B[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
umat
operator>
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;

  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
    
  arma_debug_assert_same_size(A, B, "operator>");
  
  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] > B[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
umat
operator<
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;

  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
    
  arma_debug_assert_same_size(A, B, "operator<");
  
  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] < B[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator==
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] == val)
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator!=
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;

  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] != val)
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator>=
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] >= val)
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator<=
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] <= val)
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator>
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] > val)
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator<
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A[i] < val)
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator==
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  return operator==(X,val);
  }



template<typename T1>
inline
umat
operator!=
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  return operator!=(X,val);
  }



template<typename T1>
inline
umat
operator>=
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(val >= A[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator<=
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(val <= A[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator>
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(val > A[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator<
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename umat::elem_type umat_eT;
  
  const Proxy<T1> A(X.get_ref());

  umat out(A.n_rows, A.n_cols);
  
  umat_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(val < A[i])
      {
      out_mem[i] = umat_eT(1);
      }
    else
      {
      out_mem[i] = umat_eT(0);
      }
    }
  
  return out;
  }



//! @}
