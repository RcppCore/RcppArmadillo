// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
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
  const typename arma_Mat_Col_Row_only<vec_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  arma_debug_check( (num < 2), "linspace(): num must be >= 2");
  
  typedef typename vec_type::elem_type eT;
  typedef typename vec_type::pod_type   T;
  
  const u32 n_rows = (is_Row<vec_type>::value == true) ? 1   : num;
  const u32 n_cols = (is_Row<vec_type>::value == true) ? num : 1;
  
  Mat<eT> x(n_rows, n_cols);
  eT* x_mem = x.memptr();
  
  const u32 num_m1 = num - 1;
  
  if(is_non_integral<T>::value == true)
    {
    const T delta = (end-start)/T(num_m1);
    
    for(u32 i=0; i<num_m1; ++i)
      {
      x_mem[i] = eT(start + i*delta);
      }
    
    x_mem[num_m1] = eT(end);
    }
  else
    {
    const double delta = (end >= start) ? double(end-start)/double(num_m1) : -double(start-end)/double(num_m1);
    
    for(u32 i=0; i<num_m1; ++i)
      {
      x_mem[i] = eT(double(start) + i*delta);
      }
    
    x_mem[num_m1] = eT(end);
    }
  
  return x;
  }



inline
mat
linspace(const double start, const double end, const u32 num)
  {
  arma_extra_debug_sigprint();
  return linspace<mat>(start, end, num);
  }



template<typename eT, typename T1>
inline
const mtOp<u32, T1, op_find>
find(const Base<eT,T1>& X, const u32 k = 0, const char* direction = "first")
  {
  arma_extra_debug_sigprint();
  
  const char sig = direction[0];
  
  arma_debug_check( (sig != 'f' && sig != 'F' && sig != 'l' && sig != 'L'), "find(): 3rd input argument must be \"first\" or \"last\"" );
  
  const u32 type = (sig == 'f' || sig == 'F') ? 0 : 1;
  
  return mtOp<u32, T1, op_find>(X.get_ref(), k, type);
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
const mtOp<typename T1::pod_type, T1, op_real>
real(const Base<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_real>( X.get_ref() );
  }



template<typename T1>
inline
const mtOpCube<typename T1::pod_type, T1, op_real>
real(const BaseCube<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<typename T1::pod_type, T1, op_real>( X.get_ref() );
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
  
  return eOp<Mat<typename T1::pod_type>, eop_zeros>(A.get_n_rows(), A.get_n_cols());
  }



template<typename T1>
inline
const eOpCube<Cube<typename T1::pod_type>, eop_zeros>
imag(const BaseCube<typename T1::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<T1> A(X.get_ref());
  
  return eOpCube<Cube<typename T1::pod_type>, eop_zeros>(A.get_n_rows(), A.get_n_cols(), A.get_n_slices());
  }



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_imag>
imag(const Base<std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename T1::pod_type, T1, op_imag>( X.get_ref() );
  }



template<typename T1>
inline
const mtOpCube<typename T1::pod_type, T1, op_imag>
imag(const BaseCube<std::complex<typename T1::pod_type>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return mtOpCube<typename T1::pod_type, T1, op_imag>( X.get_ref() );
  }



//
// log_add

template<typename eT>
inline
typename arma_float_only<eT>::result
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
const eOpCube<T1, eop_log>
log(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_log>(A.get_ref());
  }



//
// log2

template<typename T1>
arma_inline
const eOp<T1, eop_log2>
log2(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_log2>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_log2>
log2(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_log2>(A.get_ref());
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
const eOpCube<T1, eop_log10>
log10(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_log10>(A.get_ref());
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
const eOpCube<T1, eop_exp>
exp(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_exp>(A.get_ref());
  }



// exp2

template<typename T1>
arma_inline
const eOp<T1, eop_exp2>
exp2(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_exp2>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_exp2>
exp2(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_exp2>(A.get_ref());
  }



// exp10

template<typename T1>
arma_inline
const eOp<T1, eop_exp10>
exp10(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_exp10>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_exp10>
exp10(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_exp10>(A.get_ref());
  }



//
// abs


template<typename T1>
arma_inline
const eOp<T1, eop_abs>
abs(const Base<typename T1::elem_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOp<T1, eop_abs>(X.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_abs>
abs(const BaseCube<typename T1::elem_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOpCube<T1, eop_abs>(X.get_ref());
  }



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_abs>
abs(const Base<std::complex<typename T1::pod_type>, T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return mtOp<typename T1::pod_type, T1, op_abs>( X.get_ref() );
  }



template<typename T1>
inline
const mtOpCube<typename T1::pod_type, T1, op_abs>
abs(const BaseCube< std::complex<typename T1::pod_type>,T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return mtOpCube<typename T1::pod_type, T1, op_abs>( X.get_ref() );
  }



//
// fabs

template<typename T1>
arma_inline
const eOp<T1, eop_abs>
fabs(const Base<typename T1::pod_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOp<T1, eop_abs>(X.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_abs>
fabs(const BaseCube<typename T1::pod_type,T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOpCube<T1, eop_abs>(X.get_ref());
  }



template<typename T1>
inline
const mtOp<typename T1::pod_type, T1, op_abs>
fabs(const Base<std::complex<typename T1::pod_type>, T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return mtOp<typename T1::pod_type, T1, op_abs>( X.get_ref() );
  }



template<typename T1>
arma_inline
const mtOpCube<typename T1::pod_type, T1, op_abs>
fabs(const BaseCube< std::complex<typename T1::pod_type>,T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
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
const eOpCube<T1, eop_square>
square(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_square>(A.get_ref());
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
const eOpCube<T1, eop_sqrt>
sqrt(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_sqrt>(A.get_ref());
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
const eOpCube<T1, eop_conj>
conj(const BaseCube<std::complex<typename T1::pod_type>,T1>& A)
  {
  arma_extra_debug_sigprint();

  return eOpCube<T1, eop_conj>(A.get_ref());
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
conj(const eOpCube<T1, eop_conj>& A)
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
const eOpCube<T1, eop_pow>
pow(const BaseCube<typename T1::elem_type,T1>& A, const typename T1::elem_type exponent)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_pow>(A.get_ref(), exponent);
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
const eOpCube<T1, eop_pow>
pow(const BaseCube<typename T1::elem_type,T1>& A, const typename T1::elem_type::value_type exponent)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  return eOpCube<T1, eop_pow>(A.get_ref(), eT(exponent));
  }



//! @}
