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


//! \addtogroup Base
//! @{



//! Class for static polymorphism, modelled after the "Curiously Recurring Template Pattern" (CRTP).
//! Used for type-safe downcasting in functions that restrict their input(s) to be classes that are
//! derived from Base (e.g. Mat, Op, Glue, diagview, subview).
//! A Base object can be converted to a Mat object by the unwrap class.

template<typename elem_type, typename derived>
struct Base
  {
  arma_inline const derived& get_ref() const;
  
  arma_inline const Op<derived,op_htrans>  t() const;
  arma_inline const Op<derived,op_strans> st() const;
  
  inline void print(const std::string extra_text = "") const;
  inline void print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void raw_print(const std::string extra_text = "") const;
  inline void raw_print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_deprecated inline void print_trans(const std::string extra_text = "") const;
  arma_deprecated inline void print_trans(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_deprecated inline void raw_print_trans(const std::string extra_text = "") const;
  arma_deprecated inline void raw_print_trans(std::ostream& user_stream, const std::string extra_text = "") const;
  };



//! @}
