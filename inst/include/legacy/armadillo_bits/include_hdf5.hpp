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


#if defined(ARMA_USE_HDF5)
  
  #undef  H5_USE_110_API
  #define H5_USE_110_API
  
  #if defined(__has_include)
    #if __has_include(<hdf5.h>)
      #include <hdf5.h>
    #else
      #undef ARMA_USE_HDF5
      #pragma message ("WARNING: use of HDF5 disabled; hdf5.h header not found")
    #endif
  #else
    #include <hdf5.h>
  #endif
  
  #if defined(H5_USE_16_API) || defined(H5_USE_16_API_DEFAULT)
    #pragma message ("WARNING: use of HDF5 disabled; incompatible configuration: H5_USE_16_API or H5_USE_16_API_DEFAULT")
    #undef ARMA_USE_HDF5
  #endif
  
  // // TODO
  // #if defined(H5_USE_18_API) || defined(H5_USE_18_API_DEFAULT)
  //   #pragma message ("WARNING: detected possibly incompatible configuration of HDF5: H5_USE_18_API or H5_USE_18_API_DEFAULT")
  // #endif
  
#endif
