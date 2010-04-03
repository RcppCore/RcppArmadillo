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


//! \addtogroup fn_misc
//! @{



//! \brief
//! Generate a vector with 'num' elements.
//! The values of the elements linearly increase from 'start' upto (and including) 'end'.

template<typename vec_type>
inline
vec_type
linspace
  (
  const typename vec_type::pod_type start,
  const typename vec_type::pod_type end,
  const u32 num,
  const u32 dim = 0
  )
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< (is_Mat<vec_type>::value == false) >::apply();
  
  arma_debug_check( (num < 2), "linspace(): num must be >= 2");
  
  arma_warn( (dim != 0), "linspace(): the 'dim' argument is deprecated -- please use template based specification instead" );
  
  typedef typename vec_type::elem_type eT;
  typedef typename vec_type::pod_type   T;
  
  // // this will be the default in the future:
  // const u32 n_rows = (is_Row<vec_type>::value == true) ? 1   : num;
  // const u32 n_cols = (is_Row<vec_type>::value == true) ? num : 1;
  
  // for temporary compatibility with old user code:
  const u32 n_rows = (is_Row<vec_type>::value == true) ? 1   : ( (dim == 0) ? num : 1   );
  const u32 n_cols = (is_Row<vec_type>::value == true) ? num : ( (dim == 0) ? 1   : num );
  
  const eT delta = (end-start)/T(num-1);
  
  Mat<eT> x(n_rows, n_cols);
  eT* x_mem = x.memptr();
  
  x_mem[0] = start;
  
  for(u32 i=1; i<num; ++i)
    {
    x_mem[i] = x_mem[i-1] + delta;
    }
  
  return x;
  }



inline
mat
linspace(const double start, const double end, const u32 num, const u32 dim = 0)
  {
  arma_extra_debug_sigprint();
  return linspace<mat>(start, end, num, dim);
  }



//
// real

template<typename T1>
arma_inline
const T1&
real(const Base<typename T1::pod_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.get_ref();
  }



template<typename T1>
arma_inline
const T1&
real(const BaseCube<typename T1::pod_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.get_ref();
  }



template<typename T1>
inline
Mat<typename T1::pod_type>
real(const Base<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> A(X.get_ref());
  
  Mat<T> out(A.n_rows, A.n_cols);
  
  T* out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::real(A[i]);
    }
  
  return out;
  }



template<typename T1>
inline
Cube<typename T1::pod_type>
real(const BaseCube<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const ProxyCube<T1> A(X.get_ref());
  
  Cube<T> out(A.n_rows, A.n_cols, A.n_slices);
  
  T* out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::real(A[i]);
    }
  
  return out;
  }



//
// imag

template<typename T1>
inline
const eOp<Mat<typename T1::pod_type>, eop_zeros>
imag(const Base<typename T1::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.get_ref());
  
  return eOp<Mat<typename T1::pod_type>, eop_zeros>(A.n_rows, A.n_cols);
  }



template<typename T1>
inline
const eOpCube<Cube<typename T1::pod_type>, eop_cube_zeros>
imag(const BaseCube<typename T1::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<T1> A(X.get_ref());
  
  return eOpCube<Cube<typename T1::pod_type>, eop_cube_zeros>(A.n_rows, A.n_cols, A.n_slices);
  }



template<typename T1>
inline
Mat<typename T1::pod_type>
imag(const Base<std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> A(X.get_ref());
  
  Mat<T> out(A.n_rows, A.n_cols);
  
  T* out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::imag(A[i]);
    }
  
  return out;
  }



template<typename T1>
inline
Cube<typename T1::pod_type>
imag(const BaseCube<std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const ProxyCube<T1> A(X.get_ref());
  
  Cube<T> out(A.n_rows, A.n_cols, A.n_slices);
  
  T* out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::imag(A[i]);
    }
  
  return out;
  }



//
// log_add

template<typename eT>
inline
eT
log_add(eT log_a, eT log_b)
  {
  if(log_a < log_b)
    {
    std::swap(log_a, log_b);
    }
  
  const eT negdelta = log_b - log_a;
  
  if( (negdelta < Math<eT>::log_min()) || (arma_isfinite(negdelta) == false) )
    {
    return log_a;
    }
  else
    {
    #if defined(ARMA_HAVE_LOG1P)
      return (log_a + log1p(std::exp(negdelta)));
    #else
      return (log_a + std::log(1.0 + std::exp(negdelta)));
    #endif
    }
  }



//
// log

template<typename T1>
arma_inline
const eOp<T1, eop_log>
log(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_log>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_log>
log(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_log>(A.get_ref());
  }



//
// trunc_log

template<typename T1>
arma_inline
const eOp<T1, eop_trunc_log>
trunc_log(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_trunc_log>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_trunc_log>
trunc_log(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_trunc_log>(A.get_ref());
  }



//
// log10

template<typename T1>
arma_inline
const eOp<T1, eop_log10>
log10(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_log10>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_log10>
log10(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_log10>(A.get_ref());
  }



//
// exp

template<typename T1>
arma_inline
const eOp<T1, eop_exp>
exp(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_exp>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_exp>
exp(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_exp>(A.get_ref());
  }



//
// trunc_exp

template<typename T1>
arma_inline
const eOp<T1, eop_trunc_exp>
trunc_exp(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_trunc_exp>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_trunc_exp>
trunc_exp(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_trunc_exp>(A.get_ref());
  }



//
// abs


template<typename T1>
arma_inline
const eOp<T1, eop_abs>
abs(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_abs>(X.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_abs>
abs(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_abs>(X.get_ref());
  }



template<typename T1>
inline
Mat<typename T1::pod_type>
abs(const Base< std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> A(X.get_ref());

  // if T1 is a complex matrix,
  // pod_type is the underlying type used by std::complex;
  // otherwise pod_type is the same as elem_type
  
  typedef typename T1::elem_type  in_eT;
  typedef typename T1::pod_type  out_eT;

  Mat<out_eT> out(A.n_rows, A.n_cols);
  
  out_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::abs(A[i]);
    }
  
  return out;
  }



template<typename T1>
inline
Mat<typename T1::pod_type>
abs(const BaseCube< std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<T1> A(X.get_ref());

  // if T1 is a complex matrix,
  // pod_type is the underlying type used by std::complex;
  // otherwise pod_type is the same as elem_type
  
  typedef typename T1::elem_type  in_eT;
  typedef typename T1::pod_type  out_eT;

  Cube<out_eT> out(A.n_rows, A.n_cols, A.n_slices);
  
  out_eT* out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::abs(A[i]);
    }
  
  return out;
  }



//
// fabs

template<typename T1>
arma_inline
const eOp<T1, eop_abs>
fabs(const Base<typename T1::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_abs>(X.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_abs>
fabs(const BaseCube<typename T1::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_abs>(X.get_ref());
  }



template<typename T1>
arma_inline
Mat<typename T1::pod_type>
fabs(const Base< std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return abs(X);
  }



template<typename T1>
arma_inline
Mat<typename T1::pod_type>
fabs(const BaseCube< std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return abs(X);
  }



//
// square

template<typename T1>
arma_inline
const eOp<T1, eop_square>
square(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_square>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_square>
square(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_square>(A.get_ref());
  }



//
// sqrt

template<typename T1>
arma_inline
const eOp<T1, eop_sqrt>
sqrt(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_sqrt>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_sqrt>
sqrt(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_sqrt>(A.get_ref());
  }



//
// conj

template<typename T1>
arma_inline
const T1&
conj(const Base<typename T1::pod_type,T1>& A)
  {
  arma_extra_debug_sigprint();

  return A.get_ref();
  }



template<typename T1>
arma_inline
const T1&
conj(const BaseCube<typename T1::pod_type,T1>& A)
  {
  arma_extra_debug_sigprint();

  return A.get_ref();
  }



template<typename T1>
arma_inline
const eOp<T1, eop_conj>
conj(const Base<std::complex<typename T1::pod_type>,T1>& A)
  {
  arma_extra_debug_sigprint();

  return eOp<T1, eop_conj>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_conj>
conj(const BaseCube<std::complex<typename T1::pod_type>,T1>& A)
  {
  arma_extra_debug_sigprint();

  return eOpCube<T1, eop_cube_conj>(A.get_ref());
  }



template<typename T1>
arma_inline
const T1&
conj(const eOp<T1, eop_conj>& A)
  {
  arma_extra_debug_sigprint();
  
  return A.m;
  }



template<typename T1>
arma_inline
const T1&
conj(const eOpCube<T1, eop_cube_conj>& A)
  {
  arma_extra_debug_sigprint();
  
  return A.m;
  }



// TODO: this needs a more elaborate template restriction mechanism to work properly,
//       i.e. an overloaded version of thus function should do nothing if the input type is non-complex
// 
// //! the conjugate of the transpose of a complex matrix is the same as the hermitian transpose
// template<typename T1>
// arma_inline
// const Op<T1, op_htrans>
// conj(const Op<T1, op_trans>& A)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<T1, op_htrans>(A.m);
//   }



// pow

template<typename T1>
arma_inline
const eOp<T1, eop_pow>
pow(const Base<typename T1::elem_type,T1>& A, const typename T1::elem_type exponent)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_pow>(A.get_ref(), exponent);
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_pow>
pow(const BaseCube<typename T1::elem_type,T1>& A, const typename T1::elem_type exponent)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_cube_pow>(A.get_ref(), exponent);
  }



// pow, specialised handling (non-complex exponent for complex matrices)

template<typename T1>
arma_inline
const eOp<T1, eop_pow>
pow(const Base<typename T1::elem_type,T1>& A, const typename T1::elem_type::value_type exponent)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  return eOp<T1, eop_pow>(A.get_ref(), eT(exponent));
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_pow>
pow(const BaseCube<typename T1::elem_type,T1>& A, const typename T1::elem_type::value_type exponent)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  return eOpCube<T1, eop_cube_pow>(A.get_ref(), eT(exponent));
  }



#if defined(ARMA_GOOD_COMPILER)


// pow_s32  (integer exponent)

template<typename T1>
arma_inline
const eOp<T1, eop_pow_int>
pow(const Base<typename T1::elem_type,T1>& A, const int exponent)
  {
  arma_extra_debug_sigprint();
  
  if(exponent >= 0)
    {
    return eOp<T1, eop_pow_int>(A.get_ref(), exponent, 0);
    }
  else
    {
    return eOp<T1, eop_pow_int>(A.get_ref(), -exponent, 1);
    }
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_cube_pow_int>
pow(const BaseCube<typename T1::elem_type,T1>& A, const int exponent)
  {
  arma_extra_debug_sigprint();
  
  if(exponent >= 0)
    {
    return eOpCube<T1, eop_cube_pow_int>(A.get_ref(),  exponent, 0);
    }
  else
    {
    return eOpCube<T1, eop_cube_pow_int>(A.get_ref(), -exponent, 1);
    }
  }



#endif



//! @}
