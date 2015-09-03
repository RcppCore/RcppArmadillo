// Copyright (C) 2013-2015 Conrad Sanderson
// Copyright (C) 2013-2015 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup SizeMat
//! @{



class SizeMat
  {
  public:
  
  const uword n_rows;
  const uword n_cols;
  
  inline explicit SizeMat(const uword in_n_rows, const uword in_n_cols);
  
  inline uword operator[](const uword dim) const;
  inline uword operator()(const uword dim) const;
  
  inline bool operator==(const SizeMat& s) const;
  inline bool operator!=(const SizeMat& s) const;
  
  inline SizeMat operator+(const SizeMat& s) const;
  inline SizeMat operator-(const SizeMat& s) const;
  
  inline SizeMat operator+(const uword val) const;
  inline SizeMat operator-(const uword val) const;
  
  inline SizeMat operator*(const uword val) const;
  inline SizeMat operator/(const uword val) const;
  };



//! @}
