// Copyright (C) 2013 Conrad Sanderson
// Copyright (C) 2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup distr_param
//! @{



class distr_param
  {
  public:
  
  int a_int;
  int b_int;
  
  double a_double;
  double b_double;
  
  const bool init_int;
  const bool init_double;
  
  
  inline distr_param()
    : init_int   (false)
    , init_double(false)
    {
    }
  
  
  inline explicit distr_param(const int a, const int b)
    : a_int(a)
    , b_int(b)
    , init_int   (true )
    , init_double(false)
    {
    }
  
  
  inline explicit distr_param(const double a, const double b)
    : a_double(a)
    , b_double(b)
    , init_int   (false)
    , init_double(true )
    {
    }
  };



//! @}
