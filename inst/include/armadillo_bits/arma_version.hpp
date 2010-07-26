// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup arma_version
//! @{



struct arma_version
  {
  static const unsigned int major = 0;
  static const unsigned int minor = 9;
  static const unsigned int patch = 52;
  
  static
  inline
  std::string
  as_string()
    {
    const char* nickname = "Monkey Wrench";
    
    std::stringstream ss;
    ss << arma_version::major
       << '.'
       << arma_version::minor
       << '.'
       << arma_version::patch
       << " ("
       << nickname
       << ')';
    
    return ss.str();
    }
  };



//! @}
