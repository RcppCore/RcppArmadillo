// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



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

