// Copyright (C) 2012 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup SpBase
//! @{



template<typename elem_type, typename derived>
struct SpBase
  {
  arma_inline const derived& get_ref() const;
  
  inline const SpOp<derived,spop_htrans>  t() const;  //!< Hermitian transpose
  inline const SpOp<derived,spop_htrans> ht() const;  //!< Hermitian transpose
  inline const SpOp<derived,spop_strans> st() const;  //!< simple transpose
  
  inline void print(const std::string extra_text = "") const;
  inline void print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void raw_print(const std::string extra_text = "") const;
  inline void raw_print(std::ostream& user_stream, const std::string extra_text = "") const;

  inline void print_dense(const std::string extra_text = "") const;
  inline void print_dense(std::ostream& user_stream, const std::string extra_text = "") const;

  inline void raw_print_dense(const std::string extra_text = "") const;
  inline void raw_print_dense(std::ostream& user_stream, const std::string extra_text = "") const;
  };



//! @}
