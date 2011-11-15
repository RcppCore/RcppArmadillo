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


//! \addtogroup arma_ostream
//! @{



class arma_ostream_state
  {
  private:

  const ios::fmtflags   orig_flags;
  const std::streamsize orig_precision;
  const std::streamsize orig_width;
  const char            orig_fill;


  public:

  inline arma_ostream_state(const std::ostream& o);
  
  inline void restore(std::ostream& o) const;
  };



class arma_ostream
  {
  public:
  
  template<typename eT> inline static std::streamsize modify_stream(std::ostream& o, const eT*              data, const uword n_elem);
  template<typename  T> inline static std::streamsize modify_stream(std::ostream& o, const std::complex<T>* data, const uword n_elem);
  
  template<typename eT> inline static void print_elem_zero(std::ostream& o);
  
  template<typename eT> arma_inline static void print_elem(std::ostream& o, const eT&              x);
  template<typename  T>      inline static void print_elem(std::ostream& o, const std::complex<T>& x);

  template<typename eT> inline static void print(std::ostream& o, const  Mat<eT>& m, const bool modify);
  template<typename eT> inline static void print(std::ostream& o, const Cube<eT>& m, const bool modify);
  
  template<typename oT> inline static void print(std::ostream& o, const field<oT>&         m);
  template<typename oT> inline static void print(std::ostream& o, const subview_field<oT>& m);
  };



//! @}
