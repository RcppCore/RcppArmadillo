// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// Copyright (C) 2011 Stanislav Funiak
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


struct span_alt {};


template<typename Dummy = int>
class span_base
  {
  public:
  static const span_alt all;
  };


template<typename Dummy>
const span_alt span_base<Dummy>::all = span_alt();


class span : public span_base<>
  {
  public:

  uword a;
  uword b;
  bool  whole;
  
  inline
  span()
    : whole(true)
    {
    }
  
  
  inline
  span(const span_alt&)
    : whole(true)
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
  span(const uword in_a)
    : a(in_a)
    , b(in_a)
    , whole(false)
    {
    }
  
  inline
  span(const uword in_a, const uword in_b)
    : a(in_a)
    , b(in_b)
    , whole(false)
    {
    }

  };



//! @}
