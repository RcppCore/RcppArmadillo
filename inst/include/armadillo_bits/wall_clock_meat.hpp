// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup wall_clock
//! @{


inline
wall_clock::wall_clock()
  : valid(false)
  {
  arma_extra_debug_sigprint();
  }



inline
wall_clock::~wall_clock()
  {
  arma_extra_debug_sigprint();
  }



inline
void
wall_clock::tic()
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_BOOST_DATE)
    {
    boost_time1 = boost::posix_time::microsec_clock::local_time();
    valid = true;
    }
  #else
    #if defined(ARMA_HAVE_GETTIMEOFDAY)
      {
      gettimeofday(&posix_time1, 0);
      valid = true;
      }
    #else
      {
      arma_stop("wall_clock::tic(): need Boost libraries or POSIX gettimeofday()");
      }
    #endif
  #endif
  }



inline
double
wall_clock::toc()
  {
  arma_extra_debug_sigprint();
  
  if(valid)
    {
    #if defined(ARMA_USE_BOOST_DATE)
      {
      boost_duration = boost::posix_time::microsec_clock::local_time() - boost_time1;
      return boost_duration.total_microseconds() * 1e-6;
      }
    #else
      #if defined(ARMA_HAVE_GETTIMEOFDAY)
        {
        gettimeofday(&posix_time2, 0);
        
        const double tmp_time1 = posix_time1.tv_sec + posix_time1.tv_usec * 1.0e-6;
        const double tmp_time2 = posix_time2.tv_sec + posix_time2.tv_usec * 1.0e-6;
        
        return tmp_time2 - tmp_time1;
        }
      #else
        {
        arma_stop("wall_clock::toc(): need Boost libraries or POSIX gettimeofday()");
        return 0.0;
        }
      #endif
    #endif
    }
  else
    {  
    return 0.0;
    }
  }

//! @}

