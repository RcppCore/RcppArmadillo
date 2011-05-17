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



template<const bool size_t_is_greater_or_equal_to_8_bytes>
struct deduce_u64_helper_b
  {
  };
  

template<> struct deduce_u64_helper_b<true>
  {
  typedef std::size_t u64;
  
  static const u64  max   = (sizeof(u64) >= 8) ? 0xFFFFFFFFFFFFFFFF : 0xFFFFFFFF;  // check required for silly compilers
  static const bool trunc = false;
  };


template<> struct deduce_u64_helper_b<false>
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
    
    // #error "don't know how to typedef 'u64' on this system"
    
    // use u32 as a last resort
    
    typedef u32 u64;
    
    static const u64  max   = 0xFFFFFFFF;
    static const bool trunc = true;
    
  #endif
  };


template<const bool size_of_voidptr_is_less_than_8_bytes>
struct deduce_u64_helper_a
  {
  };


template<>
struct deduce_u64_helper_a<true>
  {
  // on 32 bit systems it's not possible to allocate more than 4 Gb of memory,
  // so we don't need to worry about indexing more than 4 Gb.
  typedef u32 u64;
  
  static const u64  max   = 0xFFFFFFFF;
  static const bool trunc = true;
  };


template<>
struct deduce_u64_helper_a<false> : public deduce_u64_helper_b< (sizeof(std::size_t) >= 8) >
  {
  };



struct deduce_u64 : public deduce_u64_helper_a< (sizeof(void*) < 8) >
  {
  };



typedef deduce_u64::u64 u64;



//! @}
