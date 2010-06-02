// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup eglue_cube_misc
//! @{



class eglue_cube_plus : public eglue_cube_core<eglue_cube_plus>
  {
  public:
  
  inline static const char* text() { return "cube addition"; }
  };



class eglue_cube_minus : public eglue_cube_core<eglue_cube_minus>
  {
  public:
  
  inline static const char* text() { return "cube subtraction"; }
  };



class eglue_cube_div : public eglue_cube_core<eglue_cube_div>
  {
  public:
  
  inline static const char* text() { return "element-wise cube division"; }
  };



class eglue_cube_schur : public eglue_cube_core<eglue_cube_schur>
  {
  public:
  
  inline static const char* text() { return "element-wise cube multiplication"; }
  };



//! @}
