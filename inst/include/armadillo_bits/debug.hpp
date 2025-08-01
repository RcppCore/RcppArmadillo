// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup debug
//! @{



inline
std::ostream&
get_cout_stream()
  {
  return (ARMA_COUT_STREAM);
  }



inline
std::ostream&
get_cerr_stream()
  {
  return (ARMA_CERR_STREAM);
  }



arma_deprecated
inline
std::ostream&
get_stream_err1()
  {
  return get_cerr_stream();
  }



arma_deprecated
inline
std::ostream&
get_stream_err2()
  {
  return get_cerr_stream();
  }



arma_frown("this function does nothing; instead use ARMA_COUT_STREAM or ARMA_WARN_LEVEL; see documentation")
inline
void
set_cout_stream(const std::ostream&)
  {
  }



arma_frown("this function does nothing; instead use ARMA_CERR_STREAM or ARMA_WARN_LEVEL; see documentation")
inline
void
set_cerr_stream(const std::ostream&)
  {
  }



arma_frown("this function does nothing; instead use ARMA_CERR_STREAM or ARMA_WARN_LEVEL; see documentation")
inline
void
set_stream_err1(const std::ostream&)
  {
  }



arma_frown("this function does nothing; instead use ARMA_CERR_STREAM or ARMA_WARN_LEVEL; see documentation")
inline
void
set_stream_err2(const std::ostream&)
  {
  }



template<typename T>
arma_frown("this function does nothing; instead use ARMA_COUT_STREAM or ARMA_WARN_LEVEL; see documentation")
inline
std::ostream&
arma_cout_stream(std::ostream*)
  {
  return (ARMA_COUT_STREAM);
  }



template<typename T>
arma_frown("this function does nothing; instead use ARMA_CERR_STREAM or ARMA_WARN_LEVEL; see documentation")
inline
std::ostream&
arma_cerr_stream(std::ostream*)
  {
  return (ARMA_CERR_STREAM);
  }



//! print a message to get_cerr_stream() and throw logic_error exception
template<typename T1>
arma_cold
arma_noinline
static
void
arma_stop_logic_error(const T1& x)
  {
  #if defined(ARMA_PRINT_EXCEPTIONS)
    {
    get_cerr_stream() << "\nerror: " << x << std::endl;
    }
  #endif
  
  throw std::logic_error( std::string(x) );
  }



arma_cold
arma_noinline
static
void
arma_stop_logic_error(const char* x, const char* y)
  {
  arma_stop_logic_error( std::string(x) + std::string(y) );
  }



//! print a message to get_cerr_stream() and throw out_of_range exception
template<typename T1>
arma_cold
arma_noinline
static
void
arma_stop_bounds_error(const T1& x)
  {
  #if defined(ARMA_PRINT_EXCEPTIONS)
    {
    get_cerr_stream() << "\nerror: " << x << std::endl;
    }
  #endif
  
  throw std::out_of_range( std::string(x) );
  }



//! print a message to get_cerr_stream() and throw bad_alloc exception
template<typename T1>
arma_cold
arma_noinline
static
void
arma_stop_bad_alloc(const T1& x)
  {
  #if defined(ARMA_PRINT_EXCEPTIONS)
    {
    get_cerr_stream() << "\nerror: " << x << std::endl;
    }
  #else
    {
    arma_ignore(x);
    }
  #endif
  
  throw std::bad_alloc();
  }



//! print a message to get_cerr_stream() and throw runtime_error exception
template<typename T1>
arma_cold
arma_noinline
static
void
arma_stop_runtime_error(const T1& x)
  {
  #if defined(ARMA_PRINT_EXCEPTIONS)
    {
    get_cerr_stream() << "\nerror: " << x << std::endl;
    }
  #endif
  
  throw std::runtime_error( std::string(x) );
  }



//
// arma_print


arma_cold
inline
void
arma_print()
  {
  get_cerr_stream() << std::endl;
  }


template<typename T1>
arma_cold
arma_noinline
static
void
arma_print(const T1& x)
  {
  get_cerr_stream() << x << std::endl;
  }



template<typename T1, typename T2>
arma_cold
arma_noinline
static
void
arma_print(const T1& x, const T2& y)
  {
  get_cerr_stream() << x << y << std::endl;
  }



template<typename T1, typename T2, typename T3>
arma_cold
arma_noinline
static
void
arma_print(const T1& x, const T2& y, const T3& z)
  {
  get_cerr_stream() << x << y << z << std::endl;
  }






//
// arma_sigprint

//! print a message to the cerr stream with a preceding @ character.
//! used for printing the signature of a function
//! (see the arma_debug_sigprint macro) 
inline
void
arma_sigprint(const char* x)
  {
  get_cerr_stream() << "@ " << x;
  }



//
// arma_bktprint


inline
void
arma_bktprint()
  {
  get_cerr_stream() << std::endl;
  }


template<typename T1>
inline
void
arma_bktprint(const T1& x)
  {
  get_cerr_stream() << " [" << x << ']' << std::endl;
  }



template<typename T1, typename T2>
inline
void
arma_bktprint(const T1& x, const T2& y)
  {
  get_cerr_stream() << " [" << x << y << ']' << std::endl;
  }






//
// arma_thisprint

inline
void
arma_thisprint(const void* this_ptr)
  {
  get_cerr_stream() << " [this: " << this_ptr << ']' << std::endl;
  }



//
// arma_plain_warn


//! print a message to the warn stream
template<typename T1>
arma_cold
arma_noinline
static
void
arma_plain_warn(const T1& arg1)
  {
  get_cerr_stream() << "\nwarning: " << arg1 << std::endl;
  }


template<typename T1, typename T2>
arma_cold
arma_noinline
static
void
arma_plain_warn(const T1& arg1, const T2& arg2)
  {
  get_cerr_stream() << "\nwarning: " << arg1 << arg2 << std::endl;
  }


template<typename T1, typename T2, typename T3>
arma_cold
arma_noinline
static
void
arma_plain_warn(const T1& arg1, const T2& arg2, const T3& arg3)
  {
  get_cerr_stream() << "\nwarning: " << arg1 << arg2 << arg3 << std::endl;
  }


template<typename T1, typename T2, typename T3, typename T4>
arma_cold
arma_noinline
static
void
arma_plain_warn(const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4)
  {
  get_cerr_stream() << "\nwarning: " << arg1 << arg2 << arg3 << arg4 << std::endl;
  }



//
// arma_warn


template<typename T1>
inline
void
arma_warn(const uword level, const T1& arg1)
  {
  constexpr uword config_level = (sword(ARMA_WARN_LEVEL) > 0) ? uword(ARMA_WARN_LEVEL) : uword(0);
  
  if((config_level > 0) && (level <= config_level))  { arma_plain_warn(arg1); }
  }


template<typename T1, typename T2>
inline
void
arma_warn(const uword level, const T1& arg1, const T2& arg2)
  {
  constexpr uword config_level = (sword(ARMA_WARN_LEVEL) > 0) ? uword(ARMA_WARN_LEVEL) : uword(0);
  
  if((config_level > 0) && (level <= config_level))  { arma_plain_warn(arg1,arg2); }
  }


template<typename T1, typename T2, typename T3>
inline
void
arma_warn(const uword level, const T1& arg1, const T2& arg2, const T3& arg3)
  {
  constexpr uword config_level = (sword(ARMA_WARN_LEVEL) > 0) ? uword(ARMA_WARN_LEVEL) : uword(0);
  
  if((config_level > 0) && (level <= config_level))  { arma_plain_warn(arg1,arg2,arg3); }
  }


template<typename T1, typename T2, typename T3, typename T4>
inline
void
arma_warn(const uword level, const T1& arg1, const T2& arg2, const T3& arg3, const T4& arg4)
  {
  constexpr uword config_level = (sword(ARMA_WARN_LEVEL) > 0) ? uword(ARMA_WARN_LEVEL) : uword(0);
  
  if((config_level > 0) && (level <= config_level))  { arma_plain_warn(arg1,arg2,arg3,arg4); }
  }



//
// arma_check

//! if state is true, abort program
template<typename T1>
arma_hot
inline
void
arma_check(const bool state, const T1& x)
  {
  if(state)  { arma_stop_logic_error(arma_str::str_wrapper(x)); }
  }


template<typename Functor>
arma_hot
inline
void
arma_check(const bool state, const char* x, const Functor& fn)
  {
  if(state)  { fn(); arma_stop_logic_error(x); }
  }


arma_hot
inline
void
arma_check(const bool state, const char* x, const char* y)
  {
  if(state)  { arma_stop_logic_error(x,y); }
  }


template<typename Functor>
arma_hot
inline
void
arma_check(const bool state, const char* x, const char* y, const Functor& fn)
  {
  if(state)  { fn(); arma_stop_logic_error(x,y); }
  }


template<typename T1>
arma_hot
inline
void
arma_check_bounds(const bool state, const T1& x)
  {
  if(state)  { arma_stop_bounds_error(arma_str::str_wrapper(x)); }
  }


template<typename T1>
arma_hot
inline
void
arma_check_bad_alloc(const bool state, const T1& x)
  {
  if(state)  { arma_stop_bad_alloc(x); }
  }



//
// arma_set_error


arma_hot
arma_inline
void
arma_set_error(bool& err_state, char*& err_msg, const bool expression, const char* message)
  {
  if(expression)
    {
    err_state = true;
    err_msg   = const_cast<char*>(message);
    }
  }




//
// functions for generating strings indicating size errors

arma_cold
arma_noinline
static
std::string
arma_incompat_size_string(const uword A_n_rows, const uword A_n_cols, const uword B_n_rows, const uword B_n_cols, const char* x)
  {
  std::ostringstream tmp;
  
  tmp << x << ": incompatible matrix dimensions: " << A_n_rows << 'x' << A_n_cols << " and " << B_n_rows << 'x' << B_n_cols;
  
  return tmp.str();
  }



arma_cold
arma_noinline
static
std::string
arma_incompat_size_string(const uword A_n_rows, const uword A_n_cols, const uword A_n_slices, const uword B_n_rows, const uword B_n_cols, const uword B_n_slices, const char* x)
  {
  std::ostringstream tmp;
  
  tmp << x << ": incompatible cube dimensions: " << A_n_rows << 'x' << A_n_cols << 'x' << A_n_slices << " and " << B_n_rows << 'x' << B_n_cols << 'x' << B_n_slices;
  
  return tmp.str();
  }



template<typename eT>
arma_cold
arma_noinline
static
std::string
arma_incompat_size_string(const subview_cube<eT>& Q, const Mat<eT>& A, const char* x)
  {
  std::ostringstream tmp;
  
  tmp << x
      << ": interpreting matrix as cube with dimensions: "
      << A.n_rows << 'x' << A.n_cols << 'x' << 1
      << " or "
      << A.n_rows << 'x' << 1        << 'x' << A.n_cols
      << " or "
      << 1        << 'x' << A.n_rows << 'x' << A.n_cols
      << " is incompatible with cube dimensions: "
      << Q.n_rows << 'x' << Q.n_cols << 'x' << Q.n_slices;
      
  return tmp.str();
  }



//
// functions for checking whether two dense matrices have the same dimensions



arma_hot
arma_inline
void
arma_assert_same_size(const uword A_n_rows, const uword A_n_cols, const uword B_n_rows, const uword B_n_cols, const char* x)
  {
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



//! stop if given matrices do not have the same size
template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Mat<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



//! stop if given proxies do not have the same size
template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Proxy<eT1>& A, const Proxy<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.get_n_rows();
  const uword A_n_cols = A.get_n_cols();
  
  const uword B_n_rows = B.get_n_rows();
  const uword B_n_cols = B.get_n_cols();
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const subview<eT1>& A, const subview<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Mat<eT1>& A, const subview<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const subview<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Mat<eT1>& A, const Proxy<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.get_n_rows();
  const uword B_n_cols = B.get_n_cols();
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Proxy<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.get_n_rows();
  const uword A_n_cols = A.get_n_cols();
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Proxy<eT1>& A, const subview<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.get_n_rows();
  const uword A_n_cols = A.get_n_cols();
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const subview<eT1>& A, const Proxy<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.get_n_rows();
  const uword B_n_cols = B.get_n_cols();
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



//
// functions for checking whether two sparse matrices have the same dimensions



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const SpMat<eT1>& A, const SpMat<eT2>& B, const char* x)
  {
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



//
// functions for checking whether two cubes have the same dimensions



arma_hot
inline
void
arma_assert_same_size(const uword A_n_rows, const uword A_n_cols, const uword A_n_slices, const uword B_n_rows, const uword B_n_cols, const uword B_n_slices, const char* x)
  {
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) || (A_n_slices != B_n_slices) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, A_n_slices, B_n_rows, B_n_cols, B_n_slices, x) );
    }
  }



//! stop if given cubes do not have the same size
template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Cube<eT1>& A, const Cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Cube<eT1>& A, const subview_cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const subview_cube<eT1>& A, const Cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const subview_cube<eT1>& A, const subview_cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices))
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT, typename T1>
arma_hot
inline
void
arma_assert_same_size(const subview_cube<eT>& A, const ProxyCube<T1>& B, const char* x)
  {
  const uword A_n_rows   = A.n_rows;
  const uword A_n_cols   = A.n_cols;
  const uword A_n_slices = A.n_slices;
  
  const uword B_n_rows   = B.get_n_rows();
  const uword B_n_cols   = B.get_n_cols();
  const uword B_n_slices = B.get_n_slices();
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) || (A_n_slices != B_n_slices) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, A_n_slices, B_n_rows, B_n_cols, B_n_slices, x) );
    }
  }



//! stop if given cube proxies do not have the same size
template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const ProxyCube<eT1>& A, const ProxyCube<eT2>& B, const char* x)
  {
  const uword A_n_rows   = A.get_n_rows();
  const uword A_n_cols   = A.get_n_cols();
  const uword A_n_slices = A.get_n_slices();
  
  const uword B_n_rows   = B.get_n_rows();
  const uword B_n_cols   = B.get_n_cols();
  const uword B_n_slices = B.get_n_slices();
  
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) || (A_n_slices != B_n_slices))
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, A_n_slices, B_n_rows, B_n_cols, B_n_slices, x) );
    }
  }



//
// functions for checking whether a cube or subcube can be interpreted as a matrix (ie. single slice)



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Cube<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != 1) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, 1, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Mat<eT1>& A, const Cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (1 != B.n_slices) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, 1, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const subview_cube<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != 1) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, 1, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_same_size(const Mat<eT1>& A, const subview_cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (1 != B.n_slices) )
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, 1, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT, typename T1>
inline
void
arma_assert_cube_as_mat(const Mat<eT>& M, const T1& Q, const char* x, const bool check_compat_size)
  {
  const uword Q_n_rows   = Q.n_rows;
  const uword Q_n_cols   = Q.n_cols;
  const uword Q_n_slices = Q.n_slices;
  
  const uword M_vec_state = M.vec_state;
  
  if(M_vec_state == 0)
    {
    if( ( (Q_n_rows == 1) || (Q_n_cols == 1) || (Q_n_slices == 1) ) == false )
      {
      std::ostringstream tmp;
        
      tmp << x
          << ": can't interpret cube with dimensions "
          << Q_n_rows << 'x' << Q_n_cols << 'x' << Q_n_slices 
          << " as a matrix; one of the dimensions must be 1";
      
      arma_stop_logic_error( tmp.str() );
      }
    }
  else
    {
    if(Q_n_slices == 1)
      {
      if( (M_vec_state == 1) && (Q_n_cols != 1) )
        {
        std::ostringstream tmp;
        
        tmp << x
            << ": can't interpret cube with dimensions "
            << Q_n_rows << 'x' << Q_n_cols << 'x' << Q_n_slices
            << " as a column vector";
        
        arma_stop_logic_error( tmp.str() );
        }
      
      if( (M_vec_state == 2) && (Q_n_rows != 1) )
        {
        std::ostringstream tmp;
        
        tmp << x
            << ": can't interpret cube with dimensions "
            << Q_n_rows << 'x' << Q_n_cols << 'x' << Q_n_slices
            << " as a row vector";
        
        arma_stop_logic_error( tmp.str() );
        }
      }
    else
      {
      if( (Q_n_cols != 1) && (Q_n_rows != 1) )
        {
        std::ostringstream tmp;
        
        tmp << x
            << ": can't interpret cube with dimensions "
            << Q_n_rows << 'x' << Q_n_cols << 'x' << Q_n_slices
            << " as a vector";
        
        arma_stop_logic_error( tmp.str() );
        }
      }
    }
  
  
  if(check_compat_size)
    {
    const uword M_n_rows = M.n_rows;
    const uword M_n_cols = M.n_cols;
    
    if(M_vec_state == 0)
      {
      if(
          (
          ( (Q_n_rows == M_n_rows) && (Q_n_cols   == M_n_cols) )
          ||
          ( (Q_n_rows == M_n_rows) && (Q_n_slices == M_n_cols) )
          ||
          ( (Q_n_cols == M_n_rows) && (Q_n_slices == M_n_cols) )
          )
          == false
        )
        {
        std::ostringstream tmp;
        
        tmp << x
            << ": can't interpret cube with dimensions "
            << Q_n_rows << 'x' << Q_n_cols << 'x' << Q_n_slices
            << " as a matrix with dimensions "
            << M_n_rows << 'x' << M_n_cols;
        
        arma_stop_logic_error( tmp.str() );
        }
      }
    else
      {
      if(Q_n_slices == 1)
        {
        if( (M_vec_state == 1) && (Q_n_rows != M_n_rows) )
          {
          std::ostringstream tmp;
          
          tmp << x
              << ": can't interpret cube with dimensions "
              << Q_n_rows << 'x' << Q_n_cols << 'x' << Q_n_slices
              << " as a column vector with dimensions "
              << M_n_rows << 'x' << M_n_cols;
          
          arma_stop_logic_error( tmp.str() );
          }
        
        if( (M_vec_state == 2) && (Q_n_cols != M_n_cols) )
          {
          std::ostringstream tmp;
          
          tmp << x
              << ": can't interpret cube with dimensions "
              << Q_n_rows << 'x' << Q_n_cols << 'x' << Q_n_slices
              << " as a row vector with dimensions "
              << M_n_rows << 'x' << M_n_cols;
          
          arma_stop_logic_error( tmp.str() );
          }
        }
      else
        {
        if( ( (M_n_cols == Q_n_slices) || (M_n_rows == Q_n_slices) ) == false )
          {
          std::ostringstream tmp;
          
          tmp << x
              << ": can't interpret cube with dimensions "
              << Q_n_rows << 'x' << Q_n_cols << 'x' << Q_n_slices
              << " as a vector with dimensions "
              << M_n_rows << 'x' << M_n_cols;
          
          arma_stop_logic_error( tmp.str() );
          }
        }
      }
    }
  }



//
// functions for checking whether two matrices have dimensions that are compatible with the matrix multiply operation



arma_hot
inline
void
arma_assert_mul_size(const uword A_n_rows, const uword A_n_cols, const uword B_n_rows, const uword B_n_cols, const char* x)
  {
  if(A_n_cols != B_n_rows)
    {
    arma_stop_logic_error( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



//! stop if given matrices are incompatible for multiplication
template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_mul_size(const Mat<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  const uword A_n_cols = A.n_cols;
  const uword B_n_rows = B.n_rows;
  
  if(A_n_cols != B_n_rows)
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A_n_cols, B_n_rows, B.n_cols, x) );
    }
  }



//! stop if given matrices are incompatible for multiplication
template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_mul_size(const Mat<eT1>& A, const Mat<eT2>& B, const bool do_trans_A, const bool do_trans_B, const char* x)
  {
  const uword final_A_n_cols = (do_trans_A == false) ? A.n_cols : A.n_rows;
  const uword final_B_n_rows = (do_trans_B == false) ? B.n_rows : B.n_cols;
    
  if(final_A_n_cols != final_B_n_rows)
    {
    const uword final_A_n_rows = (do_trans_A == false) ? A.n_rows : A.n_cols;
    const uword final_B_n_cols = (do_trans_B == false) ? B.n_cols : B.n_rows;
    
    arma_stop_logic_error( arma_incompat_size_string(final_A_n_rows, final_A_n_cols, final_B_n_rows, final_B_n_cols, x) );
    }
  }



template<const bool do_trans_A, const bool do_trans_B>
arma_hot
inline
void
arma_assert_trans_mul_size(const uword A_n_rows, const uword A_n_cols, const uword B_n_rows, const uword B_n_cols, const char* x)
  {
  const uword final_A_n_cols = (do_trans_A == false) ? A_n_cols : A_n_rows;
  const uword final_B_n_rows = (do_trans_B == false) ? B_n_rows : B_n_cols;
    
  if(final_A_n_cols != final_B_n_rows)
    {
    const uword final_A_n_rows = (do_trans_A == false) ? A_n_rows : A_n_cols;
    const uword final_B_n_cols = (do_trans_B == false) ? B_n_cols : B_n_rows;
    
    arma_stop_logic_error( arma_incompat_size_string(final_A_n_rows, final_A_n_cols, final_B_n_rows, final_B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_mul_size(const Mat<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_mul_size(const subview<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
arma_hot
inline
void
arma_assert_mul_size(const subview<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    arma_stop_logic_error( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename T1>
arma_hot
inline
void
arma_assert_blas_size(const T1& A)
  {
  if(sizeof(uword) >= sizeof(blas_int))
    {
    bool overflow;
    
    overflow = (A.n_rows > ARMA_MAX_BLAS_INT);
    overflow = (A.n_cols > ARMA_MAX_BLAS_INT) || overflow;
    
    if(overflow)
      {
      arma_stop_runtime_error("integer overflow: matrix dimensions are too large for integer type used by BLAS and LAPACK");
      }
    }
  }



template<typename T1, typename T2>
arma_hot
inline
void
arma_assert_blas_size(const T1& A, const T2& B)
  {
  if(sizeof(uword) >= sizeof(blas_int))
    {
    bool overflow;
    
    overflow = (A.n_rows > ARMA_MAX_BLAS_INT);
    overflow = (A.n_cols > ARMA_MAX_BLAS_INT) || overflow;
    overflow = (B.n_rows > ARMA_MAX_BLAS_INT) || overflow;
    overflow = (B.n_cols > ARMA_MAX_BLAS_INT) || overflow;
    
    if(overflow)
      {
      arma_stop_runtime_error("integer overflow: matrix dimensions are too large for integer type used by BLAS and LAPACK");
      }
    }
  }



// TODO: remove support for ATLAS in next major version
template<typename T1>
arma_hot
inline
void
arma_assert_atlas_size(const T1& A)
  {
  if(sizeof(uword) >= sizeof(int))
    {
    bool overflow;
    
    overflow = (A.n_rows > INT_MAX);
    overflow = (A.n_cols > INT_MAX) || overflow;
    
    if(overflow)
      {
      arma_stop_runtime_error("integer overflow: matrix dimensions are too large for integer type used by ATLAS");
      }
    }
  }



// TODO: remove support for ATLAS in next major version
template<typename T1, typename T2>
arma_hot
inline
void
arma_assert_atlas_size(const T1& A, const T2& B)
  {
  if(sizeof(uword) >= sizeof(int))
    {
    bool overflow;
    
    overflow = (A.n_rows > INT_MAX);
    overflow = (A.n_cols > INT_MAX) || overflow;
    overflow = (B.n_rows > INT_MAX) || overflow;
    overflow = (B.n_cols > INT_MAX) || overflow;
    
    if(overflow)
      {
      arma_stop_runtime_error("integer overflow: matrix dimensions are too large for integer type used by ATLAS");
      }
    }
  }



//
// macros


// #define ARMA_STRING1(x) #x
// #define ARMA_STRING2(x) ARMA_STRING1(x)
// #define ARMA_FILELINE  __FILE__ ": " ARMA_STRING2(__LINE__)


#if defined(ARMA_CHECK_CONFORMANCE)
  
  #define arma_conform_check                 arma_check
  #define arma_conform_check_bounds          arma_check_bounds
  #define arma_conform_set_error             arma_set_error
  #define arma_conform_assert_same_size      arma_assert_same_size
  #define arma_conform_assert_mul_size       arma_assert_mul_size
  #define arma_conform_assert_trans_mul_size arma_assert_trans_mul_size
  #define arma_conform_assert_cube_as_mat    arma_assert_cube_as_mat
  #define arma_conform_assert_blas_size      arma_assert_blas_size
  #define arma_conform_assert_atlas_size     arma_assert_atlas_size
  
#else
  
  #define arma_conform_check                   true ? (void)0 : arma_check
  #define arma_conform_check_bounds            true ? (void)0 : arma_check_bounds
  #define arma_conform_set_error               true ? (void)0 : arma_set_error
  #define arma_conform_assert_same_size        true ? (void)0 : arma_assert_same_size
  #define arma_conform_assert_mul_size         true ? (void)0 : arma_assert_mul_size
  #define arma_conform_assert_trans_mul_size   true ? (void)0 : arma_assert_trans_mul_size
  #define arma_conform_assert_cube_as_mat      true ? (void)0 : arma_assert_cube_as_mat
  #define arma_conform_assert_blas_size        true ? (void)0 : arma_assert_blas_size
  #define arma_conform_assert_atlas_size       true ? (void)0 : arma_assert_atlas_size
  
#endif



#if defined(ARMA_DEBUG)
  
  #define arma_debug_sigprint       arma_sigprint(ARMA_FNSIG); arma_bktprint
  #define arma_debug_sigprint_this  arma_sigprint(ARMA_FNSIG); arma_thisprint
  #define arma_debug_print          arma_print
  
  // for compatibility with earlier versions of Armadillo
  #define arma_extra_debug_sigprint       arma_sigprint(ARMA_FNSIG); arma_bktprint
  #define arma_extra_debug_sigprint_this  arma_sigprint(ARMA_FNSIG); arma_thisprint
  #define arma_extra_debug_print          arma_print
  
#else
  
  #define arma_debug_sigprint        true ? (void)0 : arma_bktprint
  #define arma_debug_sigprint_this   true ? (void)0 : arma_thisprint
  #define arma_debug_print           true ? (void)0 : arma_print
  
  // for compatibility with earlier versions of Armadillo
  #define arma_extra_debug_sigprint        true ? (void)0 : arma_bktprint
  #define arma_extra_debug_sigprint_this   true ? (void)0 : arma_thisprint
  #define arma_extra_debug_print           true ? (void)0 : arma_print
  
#endif


// for compatibility with earlier versions of Armadillo
arma_frown("use arma_conform_check() instead")
inline void arma_debug_check(bool state, const char* msg)  { arma_conform_check(state, msg); }


#if defined(ARMA_DEBUG)

  namespace junk
    {
    class arma_first_debug_message
      {
      public:
      
      inline
      arma_first_debug_message()
        {
        const char* nickname = ARMA_VERSION_NAME;
        
        std::ostream& out = get_cerr_stream();
        
        out << "@ ---" << '\n';
        out << "@ Armadillo "
            << arma_version::major << '.' << arma_version::minor << '.' << arma_version::patch
            << " (" << nickname << ')';
        
        out << "\n@ arma_config::wrapper          = " << arma_config::wrapper;
        out << "\n@ arma_config::cxx14            = " << arma_config::cxx14;
        out << "\n@ arma_config::cxx17            = " << arma_config::cxx17;
        out << "\n@ arma_config::cxx20            = " << arma_config::cxx20;
        out << "\n@ arma_config::cxx23            = " << arma_config::cxx23;
        out << "\n@ arma_config::std_mutex        = " << arma_config::std_mutex;
        out << "\n@ arma_config::posix            = " << arma_config::posix;
        out << "\n@ arma_config::openmp           = " << arma_config::openmp;
        out << "\n@ arma_config::lapack           = " << arma_config::lapack;
        out << "\n@ arma_config::blas             = " << arma_config::blas;
        out << "\n@ arma_config::newarp           = " << arma_config::newarp;
        out << "\n@ arma_config::arpack           = " << arma_config::arpack;
        out << "\n@ arma_config::superlu          = " << arma_config::superlu;
        out << "\n@ arma_config::atlas            = " << arma_config::atlas;
        out << "\n@ arma_config::hdf5             = " << arma_config::hdf5;
        out << "\n@ arma_config::good_comp        = " << arma_config::good_comp;
        out << "\n@ arma_config::extra_code       = " << arma_config::extra_code;
        out << "\n@ arma_config::hidden_args      = " << arma_config::hidden_args;
        out << "\n@ arma_config::mat_prealloc     = " << arma_config::mat_prealloc;
        out << "\n@ arma_config::mp_threshold     = " << arma_config::mp_threshold;
        out << "\n@ arma_config::mp_threads       = " << arma_config::mp_threads;
        out << "\n@ arma_config::optimise_band    = " << arma_config::optimise_band;
        out << "\n@ arma_config::optimise_sym     = " << arma_config::optimise_sym;
        out << "\n@ arma_config::optimise_invexpr = " << arma_config::optimise_invexpr;
        out << "\n@ arma_config::optimise_powexpr = " << arma_config::optimise_powexpr;
        out << "\n@ arma_config::check_conform    = " << arma_config::check_conform;
        out << "\n@ arma_config::check_nonfinite  = " << arma_config::check_nonfinite;
        out << "\n@ arma_config::fast_math        = " << arma_config::fast_math;
        out << "\n@ sizeof(void*)    = " << sizeof(void*);
        out << "\n@ sizeof(int)      = " << sizeof(int);
        out << "\n@ sizeof(long)     = " << sizeof(long);
        out << "\n@ sizeof(uword)    = " << sizeof(uword);
        out << "\n@ sizeof(blas_int) = " << sizeof(blas_int);
        out << "\n@ ---" << std::endl;
        }
      };
    
    static arma_first_debug_message arma_first_debug_message_run;
    }

#endif



//! @}
