// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup arma_version
//! @{



#define ARMA_VERSION_MAJOR 11
#define ARMA_VERSION_MINOR 2
#define ARMA_VERSION_PATCH 4
#define ARMA_VERSION_NAME  "Classic Roast"



struct arma_version
  {
  static constexpr unsigned int major = ARMA_VERSION_MAJOR;
  static constexpr unsigned int minor = ARMA_VERSION_MINOR;
  static constexpr unsigned int patch = ARMA_VERSION_PATCH;
  
  static
  inline
  std::string
  as_string()
    {
    const char* nickname = ARMA_VERSION_NAME;
    
    std::ostringstream ss;
    
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
