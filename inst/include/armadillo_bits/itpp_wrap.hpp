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


//! \addtogroup itpp_wrap
//! @{


#ifdef ARMA_USE_ITPP

  #include <itpp/base/mat.h>
  #include <itpp/base/vec.h>

#else

  namespace itpp
    {
    
    //! dummy itpp::Mat class, to prevent compilation errors
    template<typename eT>
    class Mat
      {
      public:
      
      ~Mat()                      {}
      Mat()                       {}
      Mat(int n_rows, int n_cols) {}
      Mat(const Mat& m)           {}
      
      const Mat& operator=(const Mat& m) { return *this; }

            int rows()  const { return 0; }
            int cols()  const { return 0; }
            int size()  const { return 0; }      
      const eT* _data() const { return 0; }
            eT* _data()       { return 0; }
     };
  
  
    //! dummy itpp::Vec class, to prevent compilation errors
    template<typename eT>
    class Vec
      {
      public:
      
      ~Vec()            {}
      Vec()             {}
      Vec(int length)   {}
      Vec(const Vec& m) {}
      
      const Vec& operator=(const Vec& m) { return *this; }
      
            int size()   const { return 0; }      
            int length() const { return 0; }      
      const eT* _data()  const { return 0; }
            eT* _data()        { return 0; }
      };
    
    typedef Mat<short int> smat;
    typedef Vec<short int> svec;
    
    typedef Mat<int> imat;
    typedef Vec<int> ivec;
  
    typedef Mat<double> mat;
    typedef Vec<double> vec;
    
    typedef Mat< std::complex<double> > cmat;
    typedef Vec< std::complex<double> > cvec;
    }

#endif


//! @}
