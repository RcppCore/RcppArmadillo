// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
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



#define ARMA_VERSION_MAJOR 2
#define ARMA_VERSION_MINOR 4
#define ARMA_VERSION_PATCH 2
#define ARMA_VERSION_NAME  "Loco Lounge Lizard"



struct arma_version
  {
  static const unsigned int major = ARMA_VERSION_MAJOR;
  static const unsigned int minor = ARMA_VERSION_MINOR;
  static const unsigned int patch = ARMA_VERSION_PATCH;
  
  static
  inline
  std::string
  as_string()
    {
    const char* nickname = ARMA_VERSION_NAME;
    
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
