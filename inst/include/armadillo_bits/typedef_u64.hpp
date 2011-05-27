// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup typedef_u64
//! @{


// "better than nothing" definition for u64, until the next C++ standard

template<const bool size_t_is_greater_or_equal_to_8_bytes>
struct deduce_u64
  {
  };
  

template<>
struct deduce_u64<true>
  {
  typedef std::size_t u64;
  
  static const u64  max   = (sizeof(u64) >= 8) ? 0xFFFFFFFFFFFFFFFF : 0xFFFFFFFF;  // check required for silly compilers
  static const bool trunc = false;
  };


template<>
struct deduce_u64<false>
  {
  #if (ULONG_MAX >= 0xFFFFFFFFFFFFFFFF)
    
    typedef unsigned long u64;
    
    static const u64  max   = 0xFFFFFFFFFFFFFFFF;
    static const bool trunc = false;
    
  #elif defined(ULLONG_MAX)
    
    typedef unsigned long long u64;
    
    static const u64  max   = 0xFFFFFFFFFFFFFFFF;
    static const bool trunc = false;
  
  #elif (_MSC_VER >= 1200)
  //#elif (_MSC_VER >= 1310) && defined(_MSC_EXTENSIONS)
    
    typedef unsigned __int64 u64;
    
    static const u64  max   = 0xFFFFFFFFFFFFFFFF;
    static const bool trunc = false;
  
  #else
    
    #error "don't know how to typedef 'u64' on this system"
    
  #endif
  };


typedef deduce_u64< (sizeof(std::size_t) >= 8) >::u64 u64;



//! @}
