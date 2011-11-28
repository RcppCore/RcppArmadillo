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


#if defined(ARMA_USE_ATLAS)
  #if !defined(ARMA_ATLAS_INCLUDE_DIR)
    extern "C"
      {
      #include <cblas.h>
      #include <clapack.h>
      }
  #else
    #define ARMA_STR1(x) x
    #define ARMA_STR2(x) ARMA_STR1(x)
    
    #define ARMA_CBLAS   ARMA_STR2(ARMA_ATLAS_INCLUDE_DIR)ARMA_STR2(cblas.h)
    #define ARMA_CLAPACK ARMA_STR2(ARMA_ATLAS_INCLUDE_DIR)ARMA_STR2(clapack.h)
    
    extern "C"
      {
      #include ARMA_INCFILE_WRAP(ARMA_CBLAS)
      #include ARMA_INCFILE_WRAP(ARMA_CLAPACK)
      }
    
    #undef ARMA_STR1
    #undef ARMA_STR2
    #undef ARMA_CBLAS
    #undef ARMA_CLAPACK
  #endif
#endif
