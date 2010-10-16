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


//! \addtogroup Base
//! @{



//! Class for static polymorphism, modelled after the "Curiously Recurring Template Pattern" (CRTP).
//! Used for type-safe downcasting in functions that restrict their input(s) to be classes that are
//! derived from Base (e.g. Mat, Op, Glue, diagview, subview).
//! A Base object can be converted to a Mat object by the unwrap class.

template<typename elem_type, typename derived>
struct Base
  {
  
  arma_inline
  const derived&
  get_ref() const
    {
    return static_cast<const derived&>(*this);
    }

  };



template<typename elem_type, typename derived>
struct BaseVec
  {
  
  arma_inline
  const derived&
  get_ref() const
    {
    return static_cast<const derived&>(*this);
    }

  };



//! @}
