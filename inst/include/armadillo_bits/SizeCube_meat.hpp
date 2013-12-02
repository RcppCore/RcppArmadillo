// Copyright (C) 2013 Conrad Sanderson
// Copyright (C) 2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup SizeCube
//! @{



inline
SizeCube::SizeCube(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices)
  : n_rows  (in_n_rows  )
  , n_cols  (in_n_cols  )
  , n_slices(in_n_slices)
  {
  arma_extra_debug_sigprint();
  }



inline
SizeCube::operator SizeMat () const 
  {
  arma_debug_check( (n_slices != 1), "SizeCube: n_slices != 1, hence cube size cannot be interpreted as matrix size" );
  
  return SizeMat(n_rows, n_cols);
  }



inline
bool
SizeCube::operator==(const SizeCube& s) const
  {
  if(n_rows   != s.n_rows  )  { return false; }
  
  if(n_cols   != s.n_cols  )  { return false; }
  
  if(n_slices != s.n_slices)  { return false; }
  
  return true;
  }



inline
bool
SizeCube::operator!=(const SizeCube& s) const
  {
  if(n_rows   != s.n_rows  )  { return true; }
  
  if(n_cols   != s.n_cols  )  { return true; }
  
  if(n_slices != s.n_slices)  { return true; }
  
  return false;
  }



inline
bool
SizeCube::operator==(const SizeMat& s) const
  {
  if(n_rows   != s.n_rows)  { return false; }
  
  if(n_cols   != s.n_cols)  { return false; }
  
  if(n_slices != uword(1))  { return false; }
  
  return true;
  }



inline
bool
SizeCube::operator!=(const SizeMat& s) const
  {
  if(n_rows   != s.n_rows)  { return true; }
  
  if(n_cols   != s.n_cols)  { return true; }
  
  if(n_slices != uword(1))  { return true; }
  
  return false;
  }



//! @}
