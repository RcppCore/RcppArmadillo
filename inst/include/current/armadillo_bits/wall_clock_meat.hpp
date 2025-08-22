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


inline
wall_clock::wall_clock()
  {
  arma_debug_sigprint();
  
  tic_point = std::chrono::steady_clock::now();  // warmup
  }



inline
wall_clock::~wall_clock()
  {
  arma_debug_sigprint();
  }



inline
void
wall_clock::tic()
  {
  arma_debug_sigprint();
  
  if(is_frozen)
    {
    is_frozen = false;
    
    frozen_span = std::chrono::duration<double>::zero();
    }
  
  is_started = true;
  
  tic_point = std::chrono::steady_clock::now();
  }



inline
double
wall_clock::toc()
  {
  arma_debug_sigprint();
  
  typedef std::chrono::duration<double> duration_type;
  
  const std::chrono::steady_clock::time_point toc_point = std::chrono::steady_clock::now();
  
  duration_type local_frozen_span = frozen_span;
  
  if(is_frozen)
    {
    // treat toc_point as equivalent to thaw_point
    
    const duration_type chrono_span = std::chrono::duration_cast< duration_type >(toc_point - freeze_point);
    
    local_frozen_span += chrono_span;
    }
  
  duration_type chrono_span = std::chrono::duration_cast< duration_type >(toc_point - tic_point);
  
  chrono_span -= local_frozen_span;
  
  return (is_started) ? double(chrono_span.count()) : double(0);
  }


inline
void
wall_clock::freeze()
  {
  arma_debug_sigprint();
  
  if( (is_started) && (is_frozen == false) )
    {
    is_frozen = true;
    
    freeze_point = std::chrono::steady_clock::now();
    }
  }



inline
void
wall_clock::unfreeze()
  {
  arma_debug_sigprint();
  
  typedef std::chrono::duration<double> duration_type;
  
  const std::chrono::steady_clock::time_point thaw_point = std::chrono::steady_clock::now();
  
  if(is_frozen)
    {
    is_frozen = false;
    
    const duration_type chrono_span = std::chrono::duration_cast< duration_type >(thaw_point - freeze_point);
    
    frozen_span += chrono_span;
    }
  }



//! @}
