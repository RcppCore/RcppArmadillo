// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup arma_ostream
//! @{



struct arma_ostream_state
  {
  const ios::fmtflags   orig_flags;
  const std::streamsize orig_precision;
  const std::streamsize orig_width;
  const char            orig_fill;
  
  inline arma_ostream_state(const std::ostream& o);
  
  inline void restore(std::ostream& o) const;
  };



struct arma_ostream
  {
  template<typename eT> inline static std::streamsize modify_stream(std::ostream& o, const eT*              data, const uword n_elem);
  template<typename  T> inline static std::streamsize modify_stream(std::ostream& o, const std::complex<T>* data, const uword n_elem);
  template<typename eT> inline static std::streamsize modify_stream(std::ostream& o, typename SpMat<eT>::const_iterator begin, const uword n_elem, const typename  arma_not_cx<eT>::result* junk = nullptr);
  template<typename eT> inline static std::streamsize modify_stream(std::ostream& o, typename SpMat<eT>::const_iterator begin, const uword n_elem, const typename arma_cx_only<eT>::result* junk = nullptr);
  
  template<typename eT> inline static void print_elem_zero(std::ostream& o, const bool modify);
  
  template<typename eT> inline static void     print_elem(std::ostream& o, const eT& x, const bool modify);
  template<typename eT> inline static void raw_print_elem(std::ostream& o, const eT& x);
  
  template<typename  T> inline static void     print_elem(std::ostream& o, const std::complex<T>& x, const bool modify);
  template<typename  T> inline static void raw_print_elem(std::ostream& o, const std::complex<T>& x);
  
  template<typename eT> arma_cold inline static void print(std::ostream& o, const  Mat<eT>& m, const bool modify);
  template<typename eT> arma_cold inline static void print(std::ostream& o, const Cube<eT>& m, const bool modify);
  
  template<typename oT> arma_cold inline static void print(std::ostream& o, const field<oT>&         m);
  template<typename oT> arma_cold inline static void print(std::ostream& o, const subview_field<oT>& m);
  
  template<typename eT> arma_cold inline static void print_dense(std::ostream& o, const SpMat<eT>& m, const bool modify);
  template<typename eT> arma_cold inline static void       print(std::ostream& o, const SpMat<eT>& m, const bool modify);
  
  arma_cold inline static void print(std::ostream& o, const SizeMat&  S);
  arma_cold inline static void print(std::ostream& o, const SizeCube& S);
  
  template<typename eT> arma_cold inline static void brief_print(std::ostream& o, const   Mat<eT>& m, const bool print_size = true);
  template<typename eT> arma_cold inline static void brief_print(std::ostream& o, const  Cube<eT>& m);
  template<typename eT> arma_cold inline static void brief_print(std::ostream& o, const SpMat<eT>& m);
  };



//! @}
