// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup wall_clock
//! @{


//! Class for measuring time intervals
class wall_clock
  {
  public:
  
  inline  wall_clock();
  inline ~wall_clock();
  
                   inline void   tic();  //!< reset and start the timer
  arma_warn_unused inline double toc();  //!< return the number of seconds since the last call to tic()
  
  inline void   freeze();  //!<   freeze the timer
  inline void unfreeze();  //!< unfreeze the timer
  
  private:
  
  std::chrono::steady_clock::time_point    tic_point;
  std::chrono::steady_clock::time_point freeze_point;
  
  std::chrono::duration<double> frozen_span = std::chrono::duration<double>::zero();
  
  bool is_started = false;
  bool is_frozen  = false;
  };


//! @}
