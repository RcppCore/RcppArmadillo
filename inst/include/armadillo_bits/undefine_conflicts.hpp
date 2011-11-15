// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



#if defined(log2)
  #undef log2
  
  #if defined(__GNUG__)
    #warning         "detected 'log2' macro and undefined it"
  #elif defined(_MSC_VER)
    #pragma message ("detected 'log2' macro and undefined it")
  #endif
#endif



// 
// whoever defined macros with the names "min" and "max" should be permanently removed from the gene pool

#if defined(min)
  #undef min
  
  #if defined(__GNUG__)
    #warning         "detected 'min' macro and undefined it; you may wish to define NOMINMAX before including any windows header"
  #elif defined(_MSC_VER)
    #pragma message ("detected 'min' macro and undefined it; you may wish to define NOMINMAX before including any windows header")
  #endif
#endif

#if defined(max)
  #undef max
  
  #if defined(__GNUG__)
    #warning         "detected 'max' macro and undefined it; you may wish to define NOMINMAX before including any windows header"
  #elif defined(_MSC_VER)
    #pragma message ("detected 'max' macro and undefined it; you may wish to define NOMINMAX before including any windows header")
  #endif
#endif

