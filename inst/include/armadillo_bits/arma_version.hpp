// Copyright (C) 2009-2013 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2013 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup arma_version
//! @{



#define ARMA_VERSION_MAJOR 3
#define ARMA_VERSION_MINOR 900
#define ARMA_VERSION_PATCH 7
#define ARMA_VERSION_NAME  "Bavarian Sunflower"



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
