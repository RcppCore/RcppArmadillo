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


//! \addtogroup BaseCube
//! @{



//! Analog of the Base class, intended for cubes
template<typename elem_type, typename derived>
struct BaseCube
  {
  
  arma_inline
  const derived&
  get_ref() const
    {
    return static_cast<const derived&>(*this);
    }

  };



//! @}
