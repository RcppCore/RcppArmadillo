// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// Copyright (C) 2009-2010 Ian Cullinan
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup diskio
//! @{


//! Generate the first line of the header used for saving matrices in text format.
//! Format: "ARMA_MAT_TXT_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not appliable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, e.g. "008" indicates eight bytes.
template<typename eT>
inline
std::string
diskio::gen_txt_header(const Mat<eT>& x)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));

  arma_ignore(x);
  
  if(is_u8<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IU001");
    }
  else
  if(is_s8<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IS001");
    }
  else
  if(is_u16<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IU002");
    }
  else
  if(is_s16<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IS002");
    }
  else
  if(is_u32<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IU004");
    }
  else
  if(is_s32<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IS004");
    }
#if defined(ARMA_64BIT_WORD)
  else
  if(is_u64<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IU008");
    }
  else
  if(is_s64<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IS008");
    }
#endif
  else
  if(is_float<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_FN004");
    }
  else
  if(is_double<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_FN008");
    }
  else
  if(is_complex_float<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_FC008");
    }
  else
  if(is_complex_double<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_FC016");
    }
  else
    {
    return std::string();
    }
  
  }



//! Generate the first line of the header used for saving matrices in binary format.
//! Format: "ARMA_MAT_BIN_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not appliable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, e.g. "008" indicates eight bytes.
template<typename eT>
inline
std::string
diskio::gen_bin_header(const Mat<eT>& x)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  
  arma_ignore(x);
  
  if(is_u8<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IU001");
    }
  else
  if(is_s8<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IS001");
    }
  else
  if(is_u16<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IU002");
    }
  else
  if(is_s16<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IS002");
    }
  else
  if(is_u32<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IU004");
    }
  else
  if(is_s32<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IS004");
    }
#if defined(ARMA_64BIT_WORD)
  else
  if(is_u64<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IU008");
    }
  else
  if(is_s64<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IS008");
    }
#endif
  else
  if(is_float<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_FN004");
    }
  else
  if(is_double<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_FN008");
    }
  else
  if(is_complex_float<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_FC008");
    }
  else
  if(is_complex_double<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_FC016");
    }
  else
    {
    return std::string();
    }
  
  }



//! Generate the first line of the header used for saving cubes in text format.
//! Format: "ARMA_CUB_TXT_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not appliable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, e.g. "008" indicates eight bytes.
template<typename eT>
inline
std::string
diskio::gen_txt_header(const Cube<eT>& x)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  
  arma_ignore(x);

  if(is_u8<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_IU001");
    }
  else
  if(is_s8<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_IS001");
    }
  else
  if(is_u16<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_IU002");
    }
  else
  if(is_s16<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_IS002");
    }
  else
  if(is_u32<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_IU004");
    }
  else
  if(is_s32<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_IS004");
    }
#if defined(ARMA_64BIT_WORD)
  else
  if(is_u64<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_IU008");
    }
  else
  if(is_s64<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_IS008");
    }
#endif
  else
  if(is_float<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_FN004");
    }
  else
  if(is_double<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_FN008");
    }
  else
  if(is_complex_float<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_FC008");
    }
  else
  if(is_complex_double<eT>::value == true)
    {
    return std::string("ARMA_CUB_TXT_FC016");
    }
  else
    {
    return std::string();
    }
  
  }



//! Generate the first line of the header used for saving cubes in binary format.
//! Format: "ARMA_CUB_BIN_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not appliable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, e.g. "008" indicates eight bytes.
template<typename eT>
inline
std::string
diskio::gen_bin_header(const Cube<eT>& x)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  
  arma_ignore(x);
  
  if(is_u8<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_IU001");
    }
  else
  if(is_s8<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_IS001");
    }
  else
  if(is_u16<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_IU002");
    }
  else
  if(is_s16<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_IS002");
    }
  else
  if(is_u32<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_IU004");
    }
  else
  if(is_s32<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_IS004");
    }
#if defined(ARMA_64BIT_WORD)
  else
  if(is_u64<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_IU008");
    }
  else
  if(is_s64<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_IS008");
    }
#endif
  else
  if(is_float<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_FN004");
    }
  else
  if(is_double<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_FN008");
    }
  else
  if(is_complex_float<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_FC008");
    }
  else
  if(is_complex_double<eT>::value == true)
    {
    return std::string("ARMA_CUB_BIN_FC016");
    }
  else
    {
    return std::string();
    }
  
  }



inline
file_type
diskio::guess_file_type(std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  f.clear();
  f.seekg(0, ios::end);
  
  f.clear();
  const std::fstream::pos_type pos2 = f.tellg();
  
  const uword N = ( (pos1 >= 0) && (pos2 >= 0) ) ? uword(pos2 - pos1) : 0;
  
  f.clear();
  f.seekg(pos1);
  
  podarray<unsigned char> data(N);
  
  unsigned char* ptr = data.memptr();
  
  f.clear();
  f.read( reinterpret_cast<char*>(ptr), std::streamsize(N) );
  
  const bool load_okay = f.good();
  
  f.clear();
  f.seekg(pos1);
  
  bool has_binary = false;
  bool has_comma  = false;
  
  if(load_okay == true)
    {
    uword i = 0;
    uword j = (N >= 2) ? 1 : 0;
    
    for(; j<N; i+=2, j+=2)
      {
      const unsigned char val_i = ptr[i];
      const unsigned char val_j = ptr[j];
      
      // the range checking can be made more elaborate
      if( ((val_i <= 8) || (val_i >= 123)) || ((val_j <= 8) || (val_j >= 123)) )
        {
        has_binary = true;
        break;
        }
      
      if( (val_i == ',') || (val_j == ',') )
        {
        has_comma = true;
        break;
        }
      }
    }
  else
    {
    return file_type_unknown;
    }
  
  if(has_binary)
    {
    return raw_binary;
    }
  
  if(has_comma)
    {
    return csv_ascii;
    }
  
  return raw_ascii;
  }



inline
char
diskio::conv_to_hex_char(const u8 x)
  {
  char out;

  switch(x)
    {
    case  0: out = '0'; break;
    case  1: out = '1'; break;
    case  2: out = '2'; break;
    case  3: out = '3'; break;
    case  4: out = '4'; break;
    case  5: out = '5'; break;
    case  6: out = '6'; break;
    case  7: out = '7'; break;
    case  8: out = '8'; break;
    case  9: out = '9'; break;
    case 10: out = 'a'; break;
    case 11: out = 'b'; break;
    case 12: out = 'c'; break;
    case 13: out = 'd'; break;
    case 14: out = 'e'; break;
    case 15: out = 'f'; break;
    default: out = '-'; break;
    }

  return out;  
  }



inline
void
diskio::conv_to_hex(char* out, const u8 x)
  {
  const u8 a = x / 16;
  const u8 b = x - 16*a;

  out[0] = conv_to_hex_char(a);
  out[1] = conv_to_hex_char(b);
  }



//! Append a quasi-random string to the given filename.
//! The rand() function is deliberately not used,
//! as rand() has an internal state that changes
//! from call to call. Such states should not be
//! modified in scientific applications, where the
//! results should be reproducable and not affected 
//! by saving data.
inline
std::string
diskio::gen_tmp_name(const std::string& x)
  {
  const std::string* ptr_x     = &x;
  const u8*          ptr_ptr_x = reinterpret_cast<const u8*>(&ptr_x);
  
  const char* extra      = ".tmp_";
  const uword extra_size = 5;
  
  const uword tmp_size   = 2*sizeof(u8*) + 2*2;
        char  tmp[tmp_size];
  
  uword char_count = 0;
  
  for(uword i=0; i<sizeof(u8*); ++i)
    {
    conv_to_hex(&tmp[char_count], ptr_ptr_x[i]);
    char_count += 2;
    }
  
  const uword x_size = static_cast<uword>(x.size());
  u8 sum = 0;
  
  for(uword i=0; i<x_size; ++i)
    {
    sum += u8(x[i]);
    }
  
  conv_to_hex(&tmp[char_count], sum);
  char_count += 2;
  
  conv_to_hex(&tmp[char_count], u8(x_size));
  
  
  std::string out;
  out.resize(x_size + extra_size + tmp_size);
  
  
  for(uword i=0; i<x_size; ++i)
    {
    out[i] = x[i];
    }
  
  for(uword i=0; i<extra_size; ++i)
    {
    out[x_size + i] = extra[i];
    }
  
  for(uword i=0; i<tmp_size; ++i)
    {
    out[x_size + extra_size + i] = tmp[i];
    }
  
  return out;
  }



//! Safely rename a file.
//! Before renaming, test if we can write to the final file.
//! This should prevent:
//! (i)  overwriting files that are write protected,
//! (ii) overwriting directories.
inline
bool
diskio::safe_rename(const std::string& old_name, const std::string& new_name)
  {
  std::fstream f(new_name.c_str(), std::fstream::out | std::fstream::app);
  f.put(' ');
  
  bool save_okay = f.good();
  f.close();
  
  if(save_okay == true)
    {
    std::remove(new_name.c_str());
    
    const int mv_result = std::rename(old_name.c_str(), new_name.c_str());
    
    save_okay = (mv_result == 0);
    }
  
  return save_okay;
  }



//! Save a matrix as raw text (no header, human readable).
//! Matrices can be loaded in Matlab and Octave, as long as they don't have complex elements.
template<typename eT>
inline
bool
diskio::save_raw_ascii(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out);
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_raw_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



//! Save a matrix as raw text (no header, human readable).
//! Matrices can be loaded in Matlab and Octave, as long as they don't have complex elements.
template<typename eT>
inline
bool
diskio::save_raw_ascii(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  uword cell_width;
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(10);
    cell_width = 18;
    }
  
  for(uword row=0; row < x.n_rows; ++row)
    {
    for(uword col=0; col < x.n_cols; ++col)
      {
      f.put(' ');
      
      if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
        {
        f.width(cell_width);
        }
      
      f << x.at(row,col);
      }
      
    f.put('\n');
    }
  
  return f.good();
  }



//! Save a matrix as raw binary (no header)
template<typename eT>
inline
bool
diskio::save_raw_binary(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_raw_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



template<typename eT>
inline
bool
diskio::save_raw_binary(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f.write( reinterpret_cast<const char*>(x.mem), std::streamsize(x.n_elem*sizeof(eT)) );
  
  return f.good();
  }



//! Save a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_ascii(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();

  if(save_okay == true)  
    {
    save_okay = diskio::save_arma_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



//! Save a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_ascii(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const ios::fmtflags orig_flags = f.flags();
  
  f << diskio::gen_txt_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << '\n';
  
  uword cell_width;
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(10);
    cell_width = 18;
    }
    
  for(uword row=0; row < x.n_rows; ++row)
    {
    for(uword col=0; col < x.n_cols; ++col)
      {
      f.put(' ');
      
      if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )        
        {
        f.width(cell_width);
        }
      
      f << x.at(row,col);
      }
    
    f.put('\n');
    }
  
  const bool save_okay = f.good();
  
  f.flags(orig_flags);
  
  return save_okay;
  }



//! Save a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::save_csv_ascii(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)  
    {
    save_okay = diskio::save_csv_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



//! Save a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::save_csv_ascii(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const ios::fmtflags orig_flags = f.flags();
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(10);
    }
  
  uword x_n_rows = x.n_rows;
  uword x_n_cols = x.n_cols;
  
  for(uword row=0; row < x_n_rows; ++row)
    {
    for(uword col=0; col < x_n_cols; ++col)
      {
      f << x.at(row,col);
      
      if( col < (x_n_cols-1) )
        {
        f.put(',');
        }
      }
    
    f.put('\n');
    }
  
  const bool save_okay = f.good();
  
  f.flags(orig_flags);
  
  return save_okay;
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_arma_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();

  f << diskio::gen_bin_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << '\n';
  
  f.write( reinterpret_cast<const char*>(x.mem), std::streamsize(x.n_elem*sizeof(eT)) );
  
  return f.good();
  }



//! Save a matrix as a PGM greyscale image
template<typename eT>
inline
bool
diskio::save_pgm_binary(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out | std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_pgm_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



//
// TODO:
// add functionality to save the image in a normalised format,
// i.e. scaled so that every value falls in the [0,255] range.

//! Save a matrix as a PGM greyscale image
template<typename eT>
inline
bool
diskio::save_pgm_binary(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f << "P5" << '\n';
  f << x.n_cols << ' ' << x.n_rows << '\n';
  f << 255 << '\n';
  
  const uword n_elem = x.n_rows * x.n_cols;
  podarray<u8> tmp(n_elem);
  
  uword i = 0;
  
  for(uword row=0; row < x.n_rows; ++row)
    {
    for(uword col=0; col < x.n_cols; ++col)
      {
      tmp[i] = u8( x.at(row,col) );  // TODO: add round() ?
      ++i;
      }
    }
  
  f.write(reinterpret_cast<const char*>(tmp.mem), std::streamsize(n_elem) );
  
  return f.good();
  }



//! Save a matrix as a PGM greyscale image
template<typename T>
inline
bool
diskio::save_pgm_binary(const Mat< std::complex<T> >& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const uchar_mat tmp = conv_to<uchar_mat>::from(x);
  
  return diskio::save_pgm_binary(tmp, final_name);
  }



//! Save a matrix as a PGM greyscale image
template<typename T>
inline
bool
diskio::save_pgm_binary(const Mat< std::complex<T> >& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const uchar_mat tmp = conv_to<uchar_mat>::from(x);
  
  return diskio::save_pgm_binary(tmp, f);
  }



//! Load a matrix as raw text (no header, human readable).
//! Can read matrices saved as text in Matlab and Octave.
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
bool
diskio::load_raw_ascii(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();

  std::fstream f;
  f.open(name.c_str(), std::fstream::in);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_raw_ascii(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Load a matrix as raw text (no header, human readable).
//! Can read matrices saved as text in Matlab and Octave.
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
bool
diskio::load_raw_ascii(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = f.good();
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  //
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  bool f_n_cols_found = false;
  
  std::string line_string;
  std::string token;
  
  while( (f.good() == true) && (load_okay == true) )
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)
      {
      break;
      }
    
    std::stringstream line_stream(line_string);
    
    uword line_n_cols = 0;
    
    while (line_stream >> token)
      {
      ++line_n_cols;
      }
    
    if(f_n_cols_found == false)
      {
      f_n_cols = line_n_cols;
      f_n_cols_found = true;
      }
    else
      {
      if(line_n_cols != f_n_cols)
        {
        err_msg = "inconsistent number of columns in ";
        load_okay = false;
        }
      }
    
    ++f_n_rows;
    }
    
  if(load_okay == true)
    {
    f.clear();
    f.seekg(pos1);
    
    x.set_size(f_n_rows, f_n_cols);
    
    eT val;
    
    for(uword row=0; (row < x.n_rows) && (load_okay == true); ++row)
      {
      for(uword col=0; (col < x.n_cols) && (load_okay == true); ++col)
        {
        f >> val;
        
        if(f.fail() == false)
          {
          x.at(row,col) = val;
          }
        else
          {
          load_okay = false;
          err_msg = "couldn't interpret data in ";
          //break;
          }
        }
      }
    }
  
  
  // an empty file indicates an empty matrix
  if( (f_n_cols_found == false) && (load_okay == true) )
    {
    x.reset();
    }
  
  
  return load_okay;
  }



//! Load a matrix in binary format (no header);
//! the matrix is assumed to have one column
template<typename eT>
inline
bool
diskio::load_raw_binary(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_raw_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_raw_binary(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  arma_ignore(err_msg);
  
  f.clear();
  const std::streampos pos1 = f.tellg();
  
  f.clear();
  f.seekg(0, ios::end);

  f.clear();
  const std::streampos pos2 = f.tellg();
  
  const uword N = ( (pos1 >= 0) && (pos2 >= 0) ) ? uword(pos2 - pos1) : 0;
  
  f.clear();
  //f.seekg(0, ios::beg);
  f.seekg(pos1);
  
  x.set_size(N / sizeof(eT), 1);
  
  f.clear();
  f.read( reinterpret_cast<char *>(x.memptr()), std::streamsize(N) );
  
  return f.good();
  }



//! Load a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_ascii(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f(name.c_str());
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_arma_ascii(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Load a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_ascii(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  uword f_n_rows;
  uword f_n_cols;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  
  if(f_header == diskio::gen_txt_header(x))
    {
    x.set_size(f_n_rows, f_n_cols);
    
    for(uword row=0; row < x.n_rows; ++row)
      {
      for(uword col=0; col < x.n_cols; ++col)
        {
        f >> x.at(row,col);
        }
      }
    
    load_okay = f.good();
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header in ";
    }
  
  return load_okay;
  }



//! Load a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::load_csv_ascii(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_csv_ascii(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Load a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::load_csv_ascii(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = f.good();
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  //
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  std::string line_string;
  std::string token;
  
  while( (f.good() == true) && (load_okay == true) )
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)
      {
      break;
      }
    
    std::stringstream line_stream(line_string);
    
    uword line_n_cols = 0;
    
    while(line_stream.good() == true)
      {
      getline(line_stream, token, ',');
      ++line_n_cols;
      }
    
    if(f_n_cols < line_n_cols)
      {
      f_n_cols = line_n_cols;
      }
    
    ++f_n_rows;
    }
  
  f.clear();
  f.seekg(pos1);
  
  x.zeros(f_n_rows, f_n_cols);
  
  uword row = 0;
  
  while(f.good() == true)
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)
      {
      break;
      }
    
    std::stringstream line_stream(line_string);
    
    uword col = 0;
    
    while(line_stream.good() == true)
      {
      getline(line_stream, token, ',');
      
      eT val;
      
      std::stringstream ss(token);
      
      ss >> val;
      
      if(ss.fail() == false)
        {
        x.at(row,col) = val;
        }
      
      ++col;
      }
    
    ++row;
    }
  
  return load_okay;
  }



//! Load a matrix in binary format,
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_binary(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_arma_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_arma_binary(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  uword f_n_rows;
  uword f_n_cols;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  
  if(f_header == diskio::gen_bin_header(x))
    {
    //f.seekg(1, ios::cur);  // NOTE: this may not be portable, as on a Windows machine a newline could be two characters
    f.get();
    
    x.set_size(f_n_rows,f_n_cols);
    f.read( reinterpret_cast<char *>(x.memptr()), std::streamsize(x.n_elem*sizeof(eT)) );
    
    load_okay = f.good();
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header in ";
    }
  
  return load_okay;
  }



inline
void
diskio::pnm_skip_comments(std::istream& f)
  {
  while( isspace(f.peek()) )
    {
    while( isspace(f.peek()) )
      {
      f.get();
      }
  
    if(f.peek() == '#')
      {
      while( (f.peek() != '\r') && (f.peek()!='\n') )
        {
        f.get();
        }
      }
    }
  }



//! Load a PGM greyscale image as a matrix
template<typename eT>
inline
bool
diskio::load_pgm_binary(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_pgm_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Load a PGM greyscale image as a matrix
template<typename eT>
inline
bool
diskio::load_pgm_binary(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  bool load_okay = true;
  
  std::string f_header;
  f >> f_header;
  
  if(f_header == "P5")
    {
    uword f_n_rows = 0;
    uword f_n_cols = 0;
    int f_maxval = 0;
  
    diskio::pnm_skip_comments(f);
  
    f >> f_n_cols;
    diskio::pnm_skip_comments(f);
  
    f >> f_n_rows;
    diskio::pnm_skip_comments(f);
  
    f >> f_maxval;
    f.get();
    
    if( (f_maxval > 0) || (f_maxval <= 65535) )
      {
      x.set_size(f_n_rows,f_n_cols);
      
      if(f_maxval <= 255)
        {
        const uword n_elem = f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), std::streamsize(n_elem) );
        
        uword i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        
        for(uword row=0; row < f_n_rows; ++row)
          {
          for(uword col=0; col < f_n_cols; ++col)
            {
            x.at(row,col) = eT(tmp[i]);
            ++i;
            }
          }
          
        }
      else
        {
        const uword n_elem = f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), std::streamsize(n_elem*2) );
        
        uword i = 0;
        
        for(uword row=0; row < f_n_rows; ++row)
          {
          for(uword col=0; col < f_n_cols; ++col)
            {
            x.at(row,col) = eT(tmp[i]);
            ++i;
            }
          }
        
        }
      
      }
    else
      {
      load_okay = false;
      err_msg = "currently no code available to handle loading ";
      }
    
    if(f.good() == false)
      {
      load_okay = false;
      }
    }
  else
    {
    load_okay = false;
    err_msg = "unsupported header in ";
    }
  
  return load_okay;
  }



//! Load a PGM greyscale image as a matrix
template<typename T>
inline
bool
diskio::load_pgm_binary(Mat< std::complex<T> >& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  uchar_mat tmp;
  const bool load_okay = diskio::load_pgm_binary(tmp, name, err_msg);
  
  x = conv_to< Mat< std::complex<T> > >::from(tmp);
  
  return load_okay;
  }



//! Load a PGM greyscale image as a matrix
template<typename T>
inline
bool
diskio::load_pgm_binary(Mat< std::complex<T> >& x, std::istream& is, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  uchar_mat tmp;
  const bool load_okay = diskio::load_pgm_binary(tmp, is, err_msg);
  
  x = conv_to< Mat< std::complex<T> > >::from(tmp);
  
  return load_okay;
  }



//! Try to load a matrix by automatically determining its type
template<typename eT>
inline
bool
diskio::load_auto_detect(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_auto_detect(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Try to load a matrix by automatically determining its type
template<typename eT>
inline
bool
diskio::load_auto_detect(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  static const std::string ARMA_MAT_TXT = "ARMA_MAT_TXT";
  static const std::string ARMA_MAT_BIN = "ARMA_MAT_BIN";
  static const std::string           P5 = "P5";
  
  podarray<char> raw_header(ARMA_MAT_TXT.length() + 1);
  
  std::streampos pos = f.tellg();
    
  f.read( raw_header.memptr(), std::streamsize(ARMA_MAT_TXT.length()) );
  raw_header[ARMA_MAT_TXT.length()] = '\0';
  
  f.clear();
  f.seekg(pos);
  
  const std::string header = raw_header.mem;
  
  if(ARMA_MAT_TXT == header.substr(0,ARMA_MAT_TXT.length()))
    {
    return load_arma_ascii(x, f, err_msg);
    }
  else
  if(ARMA_MAT_BIN == header.substr(0,ARMA_MAT_BIN.length()))
    {
    return load_arma_binary(x, f, err_msg);
    }
  else
  if(P5 == header.substr(0,P5.length()))
    {
    return load_pgm_binary(x, f, err_msg);
    }
  else
    {
    const file_type ft = guess_file_type(f);
    
    switch(ft)
      {
      case csv_ascii:
        return load_csv_ascii(x, f, err_msg);
        break;
      
      case raw_binary:
        return load_raw_binary(x, f, err_msg);
        break;
        
      case raw_ascii:
        return load_raw_ascii(x, f, err_msg);
        break;
      
      default:
        err_msg = "unknown data in ";
        return false;
      }
    }
  
  return false;
  }



// cubes



//! Save a cube as raw text (no header, human readable).
template<typename eT>
inline
bool
diskio::save_raw_ascii(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out);
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = save_raw_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



//! Save a cube as raw text (no header, human readable).
template<typename eT>
inline
bool
diskio::save_raw_ascii(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  uword cell_width;
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(10);
    cell_width = 18;
    }
  
  for(uword slice=0; slice < x.n_slices; ++slice)
    {
    for(uword row=0; row < x.n_rows; ++row)
      {
      for(uword col=0; col < x.n_cols; ++col)
        {
        f.put(' ');
        
        if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
          {
          f.width(cell_width);
          }
        
        f << x.at(row,col,slice);
        }
        
      f.put('\n');
      }
    }
  
  return f.good();
  }



//! Save a cube as raw binary (no header)
template<typename eT>
inline
bool
diskio::save_raw_binary(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_raw_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



template<typename eT>
inline
bool
diskio::save_raw_binary(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f.write( reinterpret_cast<const char*>(x.mem), std::streamsize(x.n_elem*sizeof(eT)) );
  
  return f.good();
  }



//! Save a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_ascii(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_arma_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



//! Save a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_ascii(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const ios::fmtflags orig_flags = f.flags();
  
  f << diskio::gen_txt_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << ' ' << x.n_slices << '\n';
  
  uword cell_width;
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(10);
    cell_width = 18;
    }
    
  for(uword slice=0; slice < x.n_slices; ++slice)
    {
    for(uword row=0; row < x.n_rows; ++row)
      {
      for(uword col=0; col < x.n_cols; ++col)
        {
        f.put(' ');
        
        if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )        
          {
          f.width(cell_width);
          }
        
        f << x.at(row,col,slice);
        }
      
      f.put('\n');
      }
    }
  
  const bool save_okay = f.good();
  
  f.flags(orig_flags);
  
  return save_okay;
  }



//! Save a cube in binary format,
//! with a header that stores the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_arma_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



//! Save a cube in binary format,
//! with a header that stores the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f << diskio::gen_bin_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << ' ' << x.n_slices << '\n';
  
  f.write( reinterpret_cast<const char*>(x.mem), std::streamsize(x.n_elem*sizeof(eT)) );
  
  return f.good();
  }



//! Load a cube as raw text (no header, human readable).
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
bool
diskio::load_raw_ascii(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  const bool load_okay = diskio::load_raw_ascii(tmp, name, err_msg);
  
  if(load_okay == true)
    {
    if(tmp.is_empty() == false)
      {
      x.set_size(tmp.n_rows, tmp.n_cols, 1);
      
      x.slice(0) = tmp;
      }
    else
      {
      x.reset();
      }
    }
  
  return load_okay;
  }



//! Load a cube as raw text (no header, human readable).
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
bool
diskio::load_raw_ascii(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  const bool load_okay = diskio::load_raw_ascii(tmp, f, err_msg);
  
  if(load_okay == true)
    {
    if(tmp.is_empty() == false)
      {
      x.set_size(tmp.n_rows, tmp.n_cols, 1);
      
      x.slice(0) = tmp;
      }
    else
      {
      x.reset();
      }
    }
  
  return load_okay;
  }



//! Load a cube in binary format (no header);
//! the cube is assumed to have one slice with one column
template<typename eT>
inline
bool
diskio::load_raw_binary(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_raw_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_raw_binary(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  arma_ignore(err_msg);
  
  f.clear();
  const std::streampos pos1 = f.tellg();
  
  f.clear();
  f.seekg(0, ios::end);
  
  f.clear();
  const std::streampos pos2 = f.tellg();
  
  const uword N = ( (pos1 >= 0) && (pos2 >= 0) ) ? uword(pos2 - pos1) : 0;
  
  f.clear();
  //f.seekg(0, ios::beg);
  f.seekg(pos1);
  
  x.set_size(N / sizeof(eT), 1, 1);
  
  f.clear();
  f.read( reinterpret_cast<char *>(x.memptr()), std::streamsize(N) );
  
  return f.good();
  }



//! Load a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_ascii(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f(name.c_str());
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_arma_ascii(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }
  


//! Load a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_ascii(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  uword f_n_rows;
  uword f_n_cols;
  uword f_n_slices;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  f >> f_n_slices;
  
  if(f_header == diskio::gen_txt_header(x))
    {
    x.set_size(f_n_rows, f_n_cols, f_n_slices);

    for(uword slice=0; slice < x.n_slices; ++slice)
      {
      for(uword row=0; row < x.n_rows; ++row)
        {
        for(uword col=0; col < x.n_cols; ++col)
          {
          f >> x.at(row,col,slice);
          }
        }
      }
    
    load_okay = f.good();
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header in ";
    }
  
  return load_okay;
  }



//! Load a cube in binary format,
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_binary(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_arma_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_arma_binary(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  uword f_n_rows;
  uword f_n_cols;
  uword f_n_slices;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  f >> f_n_slices;
  
  if(f_header == diskio::gen_bin_header(x))
    {
    //f.seekg(1, ios::cur);  // NOTE: this may not be portable, as on a Windows machine a newline could be two characters
    f.get();
    
    x.set_size(f_n_rows, f_n_cols, f_n_slices);
    f.read( reinterpret_cast<char *>(x.memptr()), std::streamsize(x.n_elem*sizeof(eT)) );
    
    load_okay = f.good();
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header in ";
    }
  
  return load_okay;
  }



//! Try to load a cube by automatically determining its type
template<typename eT>
inline
bool
diskio::load_auto_detect(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_auto_detect(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Try to load a cube by automatically determining its type
template<typename eT>
inline
bool
diskio::load_auto_detect(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  static const std::string ARMA_CUB_TXT = "ARMA_CUB_TXT";
  static const std::string ARMA_CUB_BIN = "ARMA_CUB_BIN";
  static const std::string           P6 = "P6";
  
  podarray<char> raw_header(ARMA_CUB_TXT.length() + 1);
  
  std::streampos pos = f.tellg();
  
  f.read( raw_header.memptr(), std::streamsize(ARMA_CUB_TXT.length()) );
  raw_header[ARMA_CUB_TXT.length()] = '\0';
  
  f.clear();
  f.seekg(pos);
  
  const std::string header = raw_header.mem;
  
  if(ARMA_CUB_TXT == header.substr(0, ARMA_CUB_TXT.length()))
    {
    return load_arma_ascii(x, f, err_msg);
    }
  else
  if(ARMA_CUB_BIN == header.substr(0, ARMA_CUB_BIN.length()))
    {
    return load_arma_binary(x, f, err_msg);
    }
  else
  if(P6 == header.substr(0, P6.length()))
    {
    return load_ppm_binary(x, f, err_msg);
    }
  else
    {
    const file_type ft = guess_file_type(f);
    
    switch(ft)
      {
      // case csv_ascii:
      //   return load_csv_ascii(x, f, err_msg);
      //   break;
      
      case raw_binary:
        return load_raw_binary(x, f, err_msg);
        break;
        
      case raw_ascii:
        return load_raw_ascii(x, f, err_msg);
        break;
        
      default:
        err_msg = "unknown data in ";
        return false;
      }
    }
  
  return false;
  }





// fields



template<typename T1>
inline
bool
diskio::save_arma_binary(const field<T1>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_arma_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



template<typename T1>
inline
bool
diskio::save_arma_binary(const field<T1>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( (is_Mat<T1>::value == false) && (is_Cube<T1>::value == false) ));
  
  f << "ARMA_FLD_BIN" << '\n';
  f << x.n_rows << '\n';
  f << x.n_cols << '\n';
  
  bool save_okay = true;
  
  for(uword i=0; i<x.n_elem; ++i)
    {
    save_okay = diskio::save_arma_binary(x[i], f);
    
    if(save_okay == false)
      {
      break;
      }
    }
  
  return save_okay;
  }



template<typename T1>
inline
bool
diskio::load_arma_binary(field<T1>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f( name.c_str(), std::fstream::binary );
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_arma_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename T1>
inline
bool
diskio::load_arma_binary(field<T1>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( (is_Mat<T1>::value == false) && (is_Cube<T1>::value == false) ));
  
  bool load_okay = true;
  
  std::string f_type;
  f >> f_type;
  
  if(f_type != "ARMA_FLD_BIN")
    {
    load_okay = false;
    err_msg = "unsupported field type in ";
    }
  else
    {
    uword f_n_rows;
    uword f_n_cols;
  
    f >> f_n_rows;
    f >> f_n_cols;
    
    x.set_size(f_n_rows, f_n_cols);
    
    f.get();      
    
    for(uword i=0; i<x.n_elem; ++i)
      {
      load_okay = diskio::load_arma_binary(x[i], f, err_msg);
      
      if(load_okay == false)
        {
        break;
        }
      }
    }
  
  return load_okay;
  }



inline
bool
diskio::save_std_string(const field<std::string>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_std_string(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



inline
bool
diskio::save_std_string(const field<std::string>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  for(uword row=0; row<x.n_rows; ++row)
  for(uword col=0; col<x.n_cols; ++col)
    {
    f << x.at(row,col);
    
    if(col < x.n_cols-1)
      {
      f << ' ';
      }
    else
      {
      f << '\n';
      }
    }
  
  return f.good();
  }



inline
bool
diskio::load_std_string(field<std::string>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f( name.c_str() );
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_std_string(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



inline
bool
diskio::load_std_string(field<std::string>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  //
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  bool f_n_cols_found = false;
  
  std::string line_string;
  std::string token;
  
  while( (f.good() == true) && (load_okay == true) )
    {
    std::getline(f, line_string);
    if(line_string.size() == 0)
      break;
    
    std::stringstream line_stream(line_string);
    
    uword line_n_cols = 0;
    while (line_stream >> token)
      line_n_cols++;
    
    if(f_n_cols_found == false)
      {
      f_n_cols = line_n_cols;
      f_n_cols_found = true;
      }
    else
      {
      if(line_n_cols != f_n_cols)
        {
        load_okay = false;
        err_msg = "inconsistent number of columns in ";
        }
      }
    
    ++f_n_rows;
    }
    
  if(load_okay == true)
    {
    f.clear();
    f.seekg(0, ios::beg);
    //f.seekg(start);
    
    x.set_size(f_n_rows, f_n_cols);
  
    for(uword row=0; row < x.n_rows; ++row)
      {
      for(uword col=0; col < x.n_cols; ++col)
        {
        f >> x.at(row,col);
        }
      }
    }
  
  if(f.good() == false)
    {
    load_okay = false; 
    }
  
  return load_okay;
  }



//! Try to load a field by automatically determining its type
template<typename T1>
inline
bool
diskio::load_auto_detect(field<T1>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_auto_detect(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Try to load a field by automatically determining its type
template<typename T1>
inline
bool
diskio::load_auto_detect(field<T1>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_Mat<T1>::value == false ));
  
  static const std::string ARMA_FLD_BIN = "ARMA_FLD_BIN";
  static const std::string           P6 = "P6";
  
  podarray<char> raw_header(ARMA_FLD_BIN.length() + 1);
  
  std::streampos pos = f.tellg();
  
  f.read( raw_header.memptr(), std::streamsize(ARMA_FLD_BIN.length()) );
  
  f.clear();
  f.seekg(pos);
  
  raw_header[ARMA_FLD_BIN.length()] = '\0';
  
  const std::string header = raw_header.mem;
  
  if(ARMA_FLD_BIN == header.substr(0, ARMA_FLD_BIN.length()))
    {
    return load_arma_binary(x, f, err_msg);
    }
  else
  if(P6 == header.substr(0, P6.length()))
    {
    return load_ppm_binary(x, f, err_msg);
    }
  else
    {
    err_msg = "unsupported header in ";
    return false;
    }
  }



//
// handling of PPM images by cubes


template<typename eT>
inline
bool
diskio::load_ppm_binary(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_ppm_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_ppm_binary(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  f >> f_header;
  
  if(f_header == "P6")
    {
    uword f_n_rows = 0;
    uword f_n_cols = 0;
    int f_maxval = 0;
  
    diskio::pnm_skip_comments(f);
  
    f >> f_n_cols;
    diskio::pnm_skip_comments(f);
  
    f >> f_n_rows;
    diskio::pnm_skip_comments(f);
  
    f >> f_maxval;
    f.get();
    
    if( (f_maxval > 0) || (f_maxval <= 65535) )
      {
      x.set_size(f_n_rows, f_n_cols, 3);
      
      if(f_maxval <= 255)
        {
        const uword n_elem = 3*f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), std::streamsize(n_elem) );
        
        uword i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        
        for(uword row=0; row < f_n_rows; ++row)
          {
          for(uword col=0; col < f_n_cols; ++col)
            {
            x.at(row,col,0) = eT(tmp[i+0]);
            x.at(row,col,1) = eT(tmp[i+1]);
            x.at(row,col,2) = eT(tmp[i+2]);
            i+=3;
            }
          
          }
        }
      else
        {
        const uword n_elem = 3*f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), std::streamsize(2*n_elem) );
        
        uword i = 0;
        
        for(uword row=0; row < f_n_rows; ++row)
          {
          for(uword col=0; col < f_n_cols; ++col)
            {
            x.at(row,col,0) = eT(tmp[i+0]);
            x.at(row,col,1) = eT(tmp[i+1]);
            x.at(row,col,2) = eT(tmp[i+2]);
            i+=3;
            }
          
          }
        
        }
      
      }
    else
      {
      load_okay = false;
      err_msg = "currently no code available to handle loading ";
      }
      
    if(f.good() == false)
      {
      load_okay = false;
      }
    
    }
  else
    {
    load_okay = false;
    err_msg = "unsupported header in ";
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::save_ppm_binary(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_ppm_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



template<typename eT>
inline
bool
diskio::save_ppm_binary(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (x.n_slices != 3), "diskio::save_ppm_binary(): given cube must have exactly 3 slices" );
  
  const uword n_elem = 3 * x.n_rows * x.n_cols;
  podarray<u8> tmp(n_elem);
  
  uword i = 0;
  for(uword row=0; row < x.n_rows; ++row)
    {
    for(uword col=0; col < x.n_cols; ++col)
      {
      tmp[i+0] = u8( access::tmp_real( x.at(row,col,0) ) );
      tmp[i+1] = u8( access::tmp_real( x.at(row,col,1) ) );
      tmp[i+2] = u8( access::tmp_real( x.at(row,col,2) ) );
      
      i+=3;
      }
    }
  
  f << "P6" << '\n';
  f << x.n_cols << '\n';
  f << x.n_rows << '\n';
  f << 255 << '\n';
  
  f.write( reinterpret_cast<const char*>(tmp.mem), std::streamsize(n_elem) );
  
  return f.good();
  }



//
// handling of PPM images by fields



template<typename T1>
inline
bool
diskio::load_ppm_binary(field<T1>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay == true)
    {
    load_okay = diskio::load_ppm_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename T1>
inline
bool
diskio::load_ppm_binary(field<T1>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_Mat<T1>::value == false ));
  typedef typename T1::elem_type eT;
  
  bool load_okay = true;
  
  std::string f_header;
  f >> f_header;
  
  if(f_header == "P6")
    {
    uword f_n_rows = 0;
    uword f_n_cols = 0;
    int f_maxval = 0;
  
    diskio::pnm_skip_comments(f);
  
    f >> f_n_cols;
    diskio::pnm_skip_comments(f);
  
    f >> f_n_rows;
    diskio::pnm_skip_comments(f);
  
    f >> f_maxval;
    f.get();
    
    if( (f_maxval > 0) || (f_maxval <= 65535) )
      {
      x.set_size(3);
      Mat<eT>& R = x(0);
      Mat<eT>& G = x(1);
      Mat<eT>& B = x(2);
      
      R.set_size(f_n_rows,f_n_cols);
      G.set_size(f_n_rows,f_n_cols);
      B.set_size(f_n_rows,f_n_cols);
      
      if(f_maxval <= 255)
        {
        const uword n_elem = 3*f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), std::streamsize(n_elem) );
        
        uword i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        
        for(uword row=0; row < f_n_rows; ++row)
          {
          for(uword col=0; col < f_n_cols; ++col)
            {
            R.at(row,col) = eT(tmp[i+0]);
            G.at(row,col) = eT(tmp[i+1]);
            B.at(row,col) = eT(tmp[i+2]);
            i+=3;
            }
          
          }
        }
      else
        {
        const uword n_elem = 3*f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), std::streamsize(2*n_elem) );
        
        uword i = 0;
        
        for(uword row=0; row < f_n_rows; ++row)
          {
          for(uword col=0; col < f_n_cols; ++col)
            {
            R.at(row,col) = eT(tmp[i+0]);
            G.at(row,col) = eT(tmp[i+1]);
            B.at(row,col) = eT(tmp[i+2]);
            i+=3;
            }
          
          }
        
        }
      
      }
    else
      {
      load_okay = false;
      err_msg = "currently no code available to handle loading ";
      }
    
    if(f.good() == false)
      {
      load_okay = false;
      }
    
    }
  else
    {
    load_okay = false;
    err_msg = "unsupported header in ";
    }
  
  return load_okay;
  }



template<typename T1>
inline
bool
diskio::save_ppm_binary(const field<T1>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  bool save_okay = f.is_open();
  
  if(save_okay == true)
    {
    save_okay = diskio::save_ppm_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay == true)
      {
      save_okay = diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  return save_okay;
  }



template<typename T1>
inline
bool
diskio::save_ppm_binary(const field<T1>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_Mat<T1>::value == false ));
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (x.n_elem != 3), "diskio::save_ppm_binary(): given field must have exactly 3 matrices of equal size" );
  
  bool same_size = true;
  for(uword i=1; i<3; ++i)
    {
    if( (x(0).n_rows != x(i).n_rows) || (x(0).n_cols != x(i).n_cols) )
      {
      same_size = false;
      break;
      }
    }
  
  arma_debug_check( (same_size != true), "diskio::save_ppm_binary(): given field must have exactly 3 matrices of equal size" );
  
  const Mat<eT>& R = x(0);
  const Mat<eT>& G = x(1);
  const Mat<eT>& B = x(2);
  
  f << "P6" << '\n';
  f << R.n_cols << '\n';
  f << R.n_rows << '\n';
  f << 255 << '\n';

  const uword n_elem = 3 * R.n_rows * R.n_cols;
  podarray<u8> tmp(n_elem);

  uword i = 0;
  for(uword row=0; row < R.n_rows; ++row)
    {
    for(uword col=0; col < R.n_cols; ++col)
      {
      tmp[i+0] = u8( access::tmp_real( R.at(row,col) ) );
      tmp[i+1] = u8( access::tmp_real( G.at(row,col) ) );
      tmp[i+2] = u8( access::tmp_real( B.at(row,col) ) );
      
      i+=3;
      }
    }
  
  f.write( reinterpret_cast<const char*>(tmp.mem), std::streamsize(n_elem) );
  
  return f.good();
  }



//! @}

