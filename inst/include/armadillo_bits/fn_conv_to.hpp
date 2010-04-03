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


//! \addtogroup fn_conv_to
//! @{


//
//
// conversions between various mat types

template<typename out_eT, typename in_eT>
arma_inline
void
copy_complex_elem(out_eT& out, const in_eT& in)
  {
  //arma_extra_debug_sigprint();
  
  out = out_eT(in);
  }



template<typename out_eT, typename in_T>
arma_inline
void
copy_complex_elem(out_eT& out, const std::complex<in_T>& in)
  {
  //arma_extra_debug_sigprint();
  
  out = out_eT( in.real() );
  }



template<typename out_T, typename in_T>
arma_inline
void
copy_complex_elem(std::complex<out_T>& out, const std::complex<in_T>& in)
  {
  //arma_extra_debug_sigprint();
  
  typedef std::complex<out_T> out_eT;
  out = out_eT(in);
  }



//
// scalar family


template<typename out_eT>
class conv_to
  {
  public:
  
  inline static out_eT from(const Mat< out_eT >& in);

  template<typename in_eT>
  inline static out_eT from(const Mat< in_eT >& in);

  template<typename in_eT, typename T1>
  inline static out_eT from(const Base<in_eT,T1>& in);

  //

  inline static out_eT from(const Cube< out_eT >& in);

  template<typename in_eT>
  inline static out_eT from(const Cube< in_eT >& in);

  template<typename in_eT, typename T1>
  inline static out_eT from(const BaseCube<in_eT,T1>& in);
  };



template<typename out_eT>
inline
out_eT
conv_to<out_eT>::from(const Mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_elem != 1), "conv_to<>: matrix doesn't have exactly one element" );
  
  return in.mem[0];
  }



template<typename out_eT>
template<typename in_eT>
inline
out_eT
conv_to<out_eT>::from(const Mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< is_supported_elem_type<out_eT>::value == false >::apply();
  
  arma_debug_check( (in.n_elem != 1), "conv_to<>: matrix doesn't have exactly one element" );
  
  out_eT out;
  copy_complex_elem(out, in.mem[0]);
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
out_eT
conv_to<out_eT>::from(const Base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< is_supported_elem_type<out_eT>::value == false >::apply();
  
  const unwrap<T1> tmp(in.get_ref());

  return conv_to<out_eT>::from( tmp.M );
  }



template<typename out_eT>
inline
out_eT
conv_to<out_eT>::from(const Cube<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_elem != 1), "conv_to<>: cube doesn't have exactly one element" );
  
  return in.mem[0];
  }



template<typename out_eT>
template<typename in_eT>
inline
out_eT
conv_to<out_eT>::from(const Cube<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< is_supported_elem_type<out_eT>::value == false >::apply();
  
  arma_debug_check( (in.n_elem != 1), "conv_to<>: cube doesn't have exactly one element" );
  
  out_eT out;
  copy_complex_elem(out, in.mem[0]);
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
out_eT
conv_to<out_eT>::from(const BaseCube<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< is_supported_elem_type<out_eT>::value == false >::apply();
  
  const unwrap_cube<T1> tmp(in.get_ref());

  return conv_to<out_eT>::from( tmp.M );
  }



//
//
// Mat family

template<typename out_eT>
class conv_to< Mat<out_eT> >
  {
  public:
  
  template<typename in_eT>
  inline static Mat<out_eT> from(const Mat< in_eT >& in);
  
  template<typename in_T>
  inline static Mat<out_eT> from(const Mat< std::complex<in_T> >& in);
  
  template<typename in_eT, typename T1>
  inline static Mat<out_eT> from(const Base<in_eT,T1>& in);
  
  
  
//   template<typename in_eT>
//   inline static Mat<out_eT> from(const std::vector< in_eT >& in);
//   
//   template<typename in_eT>
//   inline static Mat<out_eT> from(const std::vector< std::complex<in_eT> >& in);
  
  
  template<typename in_eT>
  inline static Mat<out_eT> from(const itpp::Mat< in_eT >& in);
  
  template<typename in_T>
  inline static Mat<out_eT> from(const itpp::Mat< std::complex<in_T> >& in);
  };



template<typename out_eT>
template<typename in_eT>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const Mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  Mat<out_eT> out(in.n_rows, in.n_cols);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const Mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  Mat<out_eT> out(in.n_rows, in.n_cols);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const Base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());

  return conv_to< Mat<out_eT> >::from( tmp.M );
  }



template<typename out_eT>
template<typename in_eT>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const itpp::Mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  Mat<out_eT> out(in.rows(), in.cols());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const itpp::Mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  Mat<out_eT> out(in.rows(), in.cols());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



//
//
// Row family

template<typename out_eT>
class conv_to< Row<out_eT> >
  {
  public:
  
  inline static Row<out_eT> from(const Mat< out_eT >& in);
  
  template<typename in_eT>
  inline static Row<out_eT> from(const Mat< in_eT >& in);
  
  template<typename in_T>  
  inline static Row<out_eT> from(const Mat< std::complex<in_T> >& in);
  
  template<typename in_eT, typename T1>
  inline static Row<out_eT> from(const Base<in_eT,T1>& in);
  
  
  
  template<typename in_eT>
  inline static Row<out_eT> from(const itpp::Vec< in_eT >& in);
  
  template<typename in_T>  
  inline static Row<out_eT> from(const itpp::Vec< std::complex<in_T> >& in);
  
  //inline static Row<out_eT> from(const Col< out_eT >& in);
  //template<typename in_eT> inline static Row<out_eT> from(const Col< in_eT >& in);
  //template<typename in_T>  inline static Row<out_eT> from(const Col< std::complex<in_T> >& in);
  
  };



template<typename out_eT>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const Mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_rows > 1), "conv_to<>: given matrix has more than one row");
  Row<out_eT> out(in.n_cols);
  
  syslib::copy_elem(out.memptr(), in.mem, out.n_elem);
  
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const Mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_rows > 1), "conv_to<>: given matrix has more than one row");
  Row<out_eT> out(in.n_cols);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const Mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  arma_debug_check( (in.n_rows > 1), "conv_to<>: given matrix has more than one row");
  Row<out_eT> out(in.n_cols);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const Base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();

  const unwrap<T1> tmp(in.get_ref());
  
  return conv_to< Row<out_eT> >::from( tmp.M );
  }



template<typename out_eT>
template<typename in_eT>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const itpp::Vec<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  Row<out_eT> out(in.length());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const itpp::Vec< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  Row<out_eT> out(in.length());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



// template<typename out_eT>
// inline
// Row<out_eT>
// conv_to< Row<out_eT> >::from(const Col<out_eT>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   return trans(in);
//   }
// 
// 
// 
// template<typename out_eT>
// template<typename in_eT>
// inline
// Row<out_eT>
// conv_to< Row<out_eT> >::from(const Col<in_eT>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   Row<out_eT> out(in.n_rows);
//   
//   const in_eT*  in_mem = in.mem;
//   out_eT*      out_mem = out.memptr();
//   
//   for(u32 i=0; i<out.n_elem; ++i)
//     {
//     out_mem[i] = out_eT( in_mem[i] );
//     }
//   
//   return out;
//   }
// 
// 
// 
// template<typename out_eT>
// template<typename in_T>
// inline
// Row<out_eT>
// conv_to< Row<out_eT> >::from(const Col< std::complex<in_T> >& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   typedef std::complex<in_T> in_eT; 
//   
//   Row<out_eT> out(in.n_rows);
//   
//   const in_eT*  in_mem = in.mem;
//   out_eT*      out_mem = out.memptr();
//   
//   for(u32 i=0; i<out.n_elem; ++i)
//     {
//     copy_complex_elem(out_mem[i], in_mem[i]);
//     }
//   
//   return out;
//   }



//
//
// Col family

template<typename out_eT>
class conv_to< Col<out_eT> >
  {
  public:
  
  inline static Col<out_eT> from(const Mat< out_eT >& in);
  
  template<typename in_eT>
  inline static Col<out_eT> from(const Mat< in_eT >& in);
  
  template<typename in_T>
  inline static Col<out_eT> from(const Mat< std::complex<in_T> >& in);
  
  template<typename in_eT, typename T1>
  inline static Col<out_eT> from(const Base<in_eT,T1>& in);
  
  
  
  template<typename in_eT>
  inline static Col<out_eT> from(const itpp::Vec< in_eT >& in);
  
  template<typename in_T>
  inline static Col<out_eT> from(const itpp::Vec< std::complex<in_T> >& in);

//   inline static Col<out_eT> from(const Row< out_eT >& in);
//   template<typename in_eT> inline static Col<out_eT> from(const Row< in_eT >& in);
//   template<typename in_T>  inline static Col<out_eT> from(const Row< std::complex<in_T> >& in);
  
  };



template<typename out_eT>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const Mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_cols > 1), "conv_to<>: given matrix has more than one column");
  Col<out_eT> out(in.n_rows);
  
  syslib::copy_elem(out.memptr(), in.mem, out.n_elem);
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const Mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_cols > 1), "conv_to<>: given matrix has more than one column");
  Col<out_eT> out(in.n_rows);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const Mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  arma_debug_check( (in.n_cols > 1), "conv_to<>: given matrix has more than one column");
  Col<out_eT> out(in.n_rows);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const Base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());
  
  return conv_to< Col<out_eT> >::from( tmp.M );
  }



template<typename out_eT>
template<typename in_eT>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const itpp::Vec<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  Col<out_eT> out(in.length());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const itpp::Vec< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  Col<out_eT> out(in.length());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



// template<typename out_eT>
// inline
// Col<out_eT>
// conv_to< Col<out_eT> >::from(const Row<out_eT>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   return trans(in);
//   }
// 
// 
// 
// template<typename out_eT>
// template<typename in_eT>
// inline
// Col<out_eT>
// conv_to< Col<out_eT> >::from(const Row<in_eT>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   Col<out_eT> out(in.n_cols);
//   
//   const in_eT*  in_mem = in.mem;
//   out_eT*      out_mem = out.memptr();
//   
//   for(u32 i=0; i<out.n_elem; ++i)
//     {
//     out_mem[i] = out_eT( in_mem[i] );
//     }
//   
//   return out;
//   }
// 
// 
// 
// template<typename out_eT>
// template<typename in_T>
// inline
// Col<out_eT>
// conv_to< Col<out_eT> >::from(const Row< std::complex<in_T> >& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   typedef std::complex<in_T> in_eT; 
//   
//   Col<out_eT> out(in.n_cols);
//   
//   const in_eT*  in_mem = in.mem;
//   out_eT*      out_mem = out.memptr();
//   
//   for(u32 i=0; i<out.n_elem; ++i)
//     {
//     copy_complex_elem(out_mem[i], in_mem[i]);
//     }
//   
//   return out;
//   }



//
//
// Cube family

template<typename out_eT>
class conv_to< Cube<out_eT> >
  {
  public:
  
  template<typename in_eT>
  inline static Cube<out_eT> from(const Cube< in_eT >& in);
  
  template<typename in_T>
  inline static Cube<out_eT> from(const Cube< std::complex<in_T> >& in);
  
  template<typename in_eT, typename T1>
  inline static Cube<out_eT> from(const BaseCube<in_eT,T1>& in);
  };



template<typename out_eT>
template<typename in_eT>
inline
Cube<out_eT>
conv_to< Cube<out_eT> >::from(const Cube<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  Cube<out_eT> out(in.n_rows, in.n_cols, in.n_slices);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
Cube<out_eT>
conv_to< Cube<out_eT> >::from(const Cube< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  Cube<out_eT> out(in.n_rows, in.n_cols, in.n_slices);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
Cube<out_eT>
conv_to< Cube<out_eT> >::from(const BaseCube<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> tmp(in.get_ref());

  return conv_to< Mat<out_eT> >::from( tmp.M );
  }



//
//
// itpp::Mat family



template<typename out_eT>
class conv_to< itpp::Mat<out_eT> >
  {
  public:
  
  inline static itpp::Mat<out_eT> from(const Mat< out_eT >& in);
  
  inline static itpp::Mat<out_eT> from(const Col< out_eT >& in);
  
  inline static itpp::Mat<out_eT> from(const Row< out_eT >& in);

  
  
  template<typename in_eT>
  inline static itpp::Mat<out_eT> from(const Mat< in_eT >& in);
  
  template<typename in_eT>
  inline static itpp::Mat<out_eT> from(const Col< in_eT >& in);
  
  template<typename in_eT>
  inline static itpp::Mat<out_eT> from(const Row< in_eT >& in);
  
  
  template<typename in_T>
  inline static itpp::Mat<out_eT> from(const Mat< std::complex<in_T> >& in);
  
  template<typename in_T>
  inline static itpp::Mat<out_eT> from(const Col< std::complex<in_T> >& in);
  
  template<typename in_T>
  inline static itpp::Mat<out_eT> from(const Row< std::complex<in_T> >& in);
  
  
  
  template<typename in_eT, typename T1>
  inline static itpp::Mat<out_eT> from(const Base<in_eT,T1>& in);
  
  };



template<typename out_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  itpp::Mat<out_eT> out(in.n_rows, in.n_cols);
  
  syslib::copy_elem(out._data(), in.mem, in.n_elem);
  
  return out;
  }



template<typename out_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Col<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const Mat<out_eT>& >(in) );
  }



template<typename out_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Row<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const Mat<out_eT>& >(in) );
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  itpp::Mat<out_eT> out(in.n_rows, in.n_cols);
  
  const in_eT* in_mem = in.memptr();
  out_eT*     out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Col<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const Mat<in_eT>& >(in) );
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Row<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const Mat<in_eT>& >(in) );
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  itpp::Mat<out_eT> out(in.n_rows, in.n_cols);
  
  const in_eT* in_mem = in.memptr();
  out_eT*     out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Col< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const Mat< std::complex<in_T> >& >(in) );
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Row< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const Mat< std::complex<in_T> >& >(in) );
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const Base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());

  return conv_to< itpp::Mat<out_eT> >::from( tmp.M );
  }



//
//
// itpp::Vec family

template<typename out_eT>
class conv_to< itpp::Vec<out_eT> >
  {
  public:
  
  inline static itpp::Vec<out_eT> from(const Mat< out_eT >& in);
  
  template<typename in_eT>
  inline static itpp::Vec<out_eT> from(const Mat< in_eT >& in);
  
  template<typename in_T>  
  inline static itpp::Vec<out_eT> from(const Mat< std::complex<in_T> >& in);
  
  
  
  template<typename in_eT>
  inline static itpp::Vec<out_eT> from(const Col< in_eT >& in);
  
  template<typename in_eT>
  inline static itpp::Vec<out_eT> from(const Row< in_eT >& in);
  
  
  
  template<typename in_T>  
  inline static itpp::Vec<out_eT> from(const Col< std::complex<in_T> >& in);
  
  template<typename in_T>  
  inline static itpp::Vec<out_eT> from(const Row< std::complex<in_T> >& in);
  
  
  
  template<typename in_eT, typename T1>
  inline static itpp::Vec<out_eT> from(const Base<in_eT,T1>& in);
  };



template<typename out_eT>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const Mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in.n_cols != 1) && (in.n_rows != 1) ), "conv_to<>: given matrix can't be interpreted as a vector");
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  syslib::copy_elem(out._data(), in.mem, in.n_elem);
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const Mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in.n_cols != 1) && (in.n_rows != 1) ), "conv_to<>: given matrix can't be interpreted as a vector");
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const Mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  arma_debug_check( ( (in.n_cols != 1) && (in.n_rows != 1) ), "conv_to<>: given matrix can't be interpreted as a vector");
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const Col<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const Col< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const Row<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const Row< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const Base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());

  return conv_to< itpp::Vec<out_eT> >::from( tmp.M );
  }


//! @}
