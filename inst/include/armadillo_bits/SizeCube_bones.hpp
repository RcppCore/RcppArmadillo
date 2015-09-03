// Copyright (C) 2013-2015 Conrad Sanderson
// Copyright (C) 2013-2015 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup SizeCube
//! @{



class SizeCube
  {
  public:
  
  const uword n_rows;
  const uword n_cols;
  const uword n_slices;
  
  inline explicit SizeCube(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices);
  
  inline uword operator[](const uword dim) const;
  inline uword operator()(const uword dim) const;
  
  inline bool operator==(const SizeCube& s) const;
  inline bool operator!=(const SizeCube& s) const;
  
  inline SizeCube operator+(const SizeCube& s) const;
  inline SizeCube operator-(const SizeCube& s) const;
  
  inline SizeCube operator+(const uword val) const;
  inline SizeCube operator-(const uword val) const;
  
  inline SizeCube operator*(const uword val) const;
  inline SizeCube operator/(const uword val) const;
  };



//! @}
