// Copyright (C) 2012 Conrad Sanderson
// Copyright (C) 2012 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


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
