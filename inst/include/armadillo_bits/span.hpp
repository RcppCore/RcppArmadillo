// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup span
//! @{



// an "anonymous" namespace is used to avoid possible multiple definitions of span::all
// clashing with each other during linking
namespace
  {
  
  class span
    {
    public:
    
    const u32  a;
    const u32  b;
    const bool whole;
    
    static const span all;
    
    inline
    span()
      : a(a)
      , b(b)
      , whole(true)
      {
      }
    
    
    // TODO:
    // if the "explicit" keyword is removed or commented out,
    // the compiler will be able to automatically convert integers to an instance of the span class.
    // this is useful for Cube::operator()(span&, span&, span&),
    // but it might have unintended consequences or interactions elsewhere.
    // as such, removal of "explicit" needs thorough testing.
    inline
    explicit
    span(const u32 in_a)
      : a(in_a)
      , b(in_a)
      , whole(false)
      {
      }
    
    inline
    span(const u32 in_a, const u32 in_b)
      : a(in_a)
      , b(in_b)
      , whole(false)
      {
      }

    };

  const span span::all = span();
  }



//! @}
