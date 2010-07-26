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


//! \addtogroup debug
//! @{


//
// arma_stop

//! print a message to std::cout and/or throw a run-time error exception
template<typename T1>
inline
void
arma_cold
arma_stop(const T1& x)
  {
  std::cerr.flush();
  std::cout.flush();
  
  std::cout << '\n';
  std::cout << "run-time error: " << x << '\n';
  std::cout << '\n';
  std::cout.flush();
  
  throw std::runtime_error("");
  }



//
// arma_print


inline
void
arma_cold
arma_print()
  {
  std::cout << std::endl;
  }


template<typename T1>
inline
void
arma_cold
arma_print(const T1& x)
  {
  std::cout << x << std::endl;
  }



template<typename T1, typename T2>
inline
void
arma_cold
arma_print(const T1& x, const T2& y)
  {
  std::cout << x << y << std::endl;
  }



template<typename T1, typename T2, typename T3>
inline
void
arma_cold
arma_print(const T1& x, const T2& y, const T3& z)
  {
  std::cout << x << y << z << std::endl;
  }






//
// arma_sigprint

//! print a message on cout, with a preceding @ character.
//! used for printing the signature of a function
//! (see the arma_extra_debug_sigprint macro) 
inline
void
arma_sigprint(const char* x)
  {
  std::cout << "@ " << x;
  }



//
// arma_bktprint


inline
void
arma_bktprint()
  {
  std::cout << std::endl;
  }


template<typename T1>
inline
void
arma_bktprint(const T1& x)
  {
  std::cout << " [" << x << ']' << std::endl;
  }



template<typename T1, typename T2>
inline
void
arma_bktprint(const T1& x, const T2& y)
  {
  std::cout << " [" << x << y << ']' << std::endl;
  }






//
// arma_thisprint


inline
void
arma_thisprint(void* this_ptr)
  {
  std::cout << " [this = " << this_ptr << ']' << std::endl;
  }



//
// arma_warn

//! if state is true, print a message on cout
template<typename T1>
inline
void
arma_hot
arma_warn(const bool state, const T1& x)
  {
  if(state==true)
    {
    arma_print(x);
    }
  }


template<typename T1, typename T2>
inline
void
arma_hot
arma_warn(const bool state, const T1& x, const T2& y)
  {
  if(state==true)
    {
    arma_print(x,y);
    }
  }





//
// arma_check

//! if state is true, abort program
template<typename T1>
inline
void
arma_hot
arma_check(const bool state, const T1& x)
  {
  if(state==true)
    {
    arma_stop(arma_boost::str_wrapper(x));
    }
  }


template<typename T1, typename T2>
inline
void
arma_hot
arma_check(const bool state, const T1& x, const T2& y)
  {
  if(state==true)
    {
    arma_stop( std::string(x) + std::string(y) );
    }
  }






//
// functions for checking whether two matrices have the same dimensions



inline
std::string
arma_cold
arma_incompat_size_string(const u32 A_n_rows, const u32 A_n_cols, const u32 B_n_rows, const u32 B_n_cols, const char* x)
  {
  std::stringstream tmp;
  
  tmp << x << ": incompatible matrix dimensions: (" << A_n_rows << ',' << A_n_cols << ") and (" << B_n_rows << ',' << B_n_cols << ')';
  
  return tmp.str();
  }



inline
void
arma_hot
arma_assert_same_size(const u32 A_n_rows, const u32 A_n_cols, const u32 B_n_rows, const u32 B_n_cols, const char* x)
  {
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    arma_stop( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



//! stop if given matrices have different sizes
template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Mat<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



//! stop if given proxies have different sizes
template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Proxy<eT1>& A, const Proxy<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const subview<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Mat<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const subview<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Mat<eT1>& A, const Proxy<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop ( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Proxy<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Proxy<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const subview<eT1>& A, const Proxy<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



//
// functions for checking whether two cubes have the same dimensions



inline
arma_cold
std::string
arma_incompat_size_string(const u32 A_n_rows, const u32 A_n_cols, const u32 A_n_slices, const u32 B_n_rows, const u32 B_n_cols, const u32 B_n_slices, const char* x)
  {
  std::stringstream tmp;
  
  tmp << x << ": incompatible cube dimensions: (" << A_n_rows << ',' << A_n_cols << ',' << A_n_slices << ") and (" << B_n_rows << ',' << B_n_cols << ',' << B_n_slices << ')';
  
  return tmp.str();
  }



inline
void
arma_hot
arma_assert_same_size(const u32 A_n_rows, const u32 A_n_cols, const u32 A_n_slices, const u32 B_n_rows, const u32 B_n_cols, const u32 B_n_slices, const char* x)
  {
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) || (A_n_slices != B_n_slices) )
    {
    arma_stop( arma_incompat_size_string(A_n_rows, A_n_cols, A_n_slices, B_n_rows, B_n_cols, B_n_slices, x) );
    }
  }



//! stop if given cubes have different sizes
template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Cube<eT1>& A, const Cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Cube<eT1>& A, const subview_cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const subview_cube<eT1>& A, const Cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const subview_cube<eT1>& A, const subview_cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices))
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



//! stop if given cube proxies have different sizes
template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const ProxyCube<eT1>& A, const ProxyCube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != B.n_slices))
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



//
// functions for checking whether a cube or subcube can be interpreted as a matrix (i.e. single slice)



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Cube<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != 1) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, 1, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Mat<eT1>& A, const Cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (1 != B.n_slices) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, 1, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const subview_cube<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (A.n_slices != 1) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, A.n_slices, B.n_rows, B.n_cols, 1, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Mat<eT1>& A, const subview_cube<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) || (1 != B.n_slices) )
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, 1, B.n_rows, B.n_cols, B.n_slices, x) );
    }
  }



//
// functions for checking whether two matrices have dimensions that are compatible with the matrix multiply operation



inline
void
arma_hot
arma_assert_mul_size(const u32 A_n_rows, const u32 A_n_cols, const u32 B_n_rows, const u32 B_n_cols, const char* x)
  {
  if(A_n_cols != B_n_rows)
    {
    arma_stop( arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x) );
    }
  }



//! stop if given matrices are incompatible for multiplication
template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const Mat<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



//! stop if given matrices are incompatible for multiplication
template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const Mat<eT1>& A, const Mat<eT2>& B, const bool do_trans_A, const bool do_trans_B, const char* x)
  {
  const u32 final_A_n_cols = (do_trans_A == false) ? A.n_cols : A.n_rows;
  const u32 final_B_n_rows = (do_trans_B == false) ? B.n_rows : B.n_cols;
    
  if(final_A_n_cols != final_B_n_rows)
    {
    const u32 final_A_n_rows = (do_trans_A == false) ? A.n_rows : A.n_cols;
    const u32 final_B_n_cols = (do_trans_B == false) ? B.n_cols : B.n_rows;
    
    arma_stop( arma_incompat_size_string(final_A_n_rows, final_A_n_cols, final_B_n_rows, final_B_n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const Mat<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const subview<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const subview<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    arma_stop( arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x) );
    }
  }



//
// macros


#define ARMA_STRING1(x) #x
#define ARMA_STRING2(x) ARMA_STRING1(x)
#define ARMA_FILELINE  __FILE__ ": " ARMA_STRING2(__LINE__)


#if defined (__GNUG__)
  #define ARMA_FNSIG  __PRETTY_FUNCTION__
#elif defined (_MSC_VER)
  #define ARMA_FNSIG  __FUNCSIG__ 
#elif defined (ARMA_USE_BOOST)
  #define ARMA_FNSIG  BOOST_CURRENT_FUNCTION  
#else 
  #define ARMA_FNSIG  "(unknown)"
#endif



#if !defined(ARMA_NO_DEBUG) && !defined(NDEBUG)
  
  #define arma_debug_print            arma_print
  #define arma_debug_warn             arma_warn
  #define arma_debug_check            arma_check
  #define arma_debug_assert_same_size arma_assert_same_size
  #define arma_debug_assert_mul_size  arma_assert_mul_size
  
#else
  
  #undef ARMA_EXTRA_DEBUG
  
  #define arma_debug_print            true ? (void)0 : arma_print
  #define arma_debug_warn             true ? (void)0 : arma_warn
  #define arma_debug_check            true ? (void)0 : arma_check
  #define arma_debug_assert_same_size true ? (void)0 : arma_assert_same_size
  #define arma_debug_assert_mul_size  true ? (void)0 : arma_assert_mul_size

#endif


#if defined(ARMA_EXTRA_DEBUG)
  
  #define arma_extra_debug_sigprint       arma_sigprint(ARMA_FNSIG); arma_bktprint
  #define arma_extra_debug_sigprint_this  arma_sigprint(ARMA_FNSIG); arma_thisprint
  #define arma_extra_debug_print          arma_print
  #define arma_extra_debug_warn           arma_warn
  #define arma_extra_debug_check          arma_check

#else
  
  #define arma_extra_debug_sigprint        true ? (void)0 : arma_bktprint
  #define arma_extra_debug_sigprint_this   true ? (void)0 : arma_thisprint
  #define arma_extra_debug_print           true ? (void)0 : arma_print
  #define arma_extra_debug_warn            true ? (void)0 : arma_warn
  #define arma_extra_debug_check           true ? (void)0 : arma_check
 
#endif




#if defined(ARMA_EXTRA_DEBUG)

  namespace junk
    {
    class arma_first_extra_debug_message
      {
      public:
      
      inline
      arma_cold
      arma_first_extra_debug_message()
        {
        std::cout << "@ ---" << '\n';
        std::cout << "@ Armadillo " << arma_version::major << '.' << arma_version::minor << '.' << arma_version::patch << '\n';
        std::cout << "@ arma_config::atlas      = " << arma_config::atlas      << '\n';
        std::cout << "@ arma_config::lapack     = " << arma_config::lapack     << '\n';
        std::cout << "@ arma_config::blas       = " << arma_config::blas       << '\n';
        std::cout << "@ arma_config::boost      = " << arma_config::boost      << '\n';
        std::cout << "@ arma_config::boost_date = " << arma_config::boost_date << '\n';
        std::cout << "@ arma_config::good_comp  = " << arma_config::good_comp  << '\n';
        std::cout << "@ sizeof(int)  = " << sizeof(int)  << '\n';
        std::cout << "@ sizeof(int*) = " << sizeof(int*) << '\n';
        std::cout << "@ sizeof(long) = " << sizeof(long) << '\n';
        std::cout << "@ ---" << std::endl;
        }
      
      };
    
    static arma_first_extra_debug_message arma_first_extra_debug_message_run;
    }

#endif


//! @}
