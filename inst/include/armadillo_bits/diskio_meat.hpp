// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Ian Cullinan (ian dot cullinan at nicta dot com dot au)
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
  arma_type_check<diskio::is_supported_type<eT>::value == false>::apply();

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
  arma_type_check<diskio::is_supported_type<eT>::value == false>::apply();
  
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
  arma_type_check<diskio::is_supported_type<eT>::value == false>::apply();

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
  arma_type_check<diskio::is_supported_type<eT>::value == false>::apply();
  
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
  const u32   extra_size = 5;
  
  const u32   tmp_size   = 2*sizeof(u8*) + 2*2;
        char  tmp[tmp_size];
  
  u32 char_count = 0;
  
  for(u32 i=0; i<sizeof(u8*); ++i)
    {
    conv_to_hex(&tmp[char_count], ptr_ptr_x[i]);
    char_count += 2;
    }
  
  const u32 x_size = x.size();
  u8 sum = 0;
  
  for(u32 i=0; i<x_size; ++i)
    {
    sum += u8(x[i]);
    }
  
  conv_to_hex(&tmp[char_count], sum);
  char_count += 2;
  
  conv_to_hex(&tmp[char_count], u8(x_size));
  
  
  std::string out;
  out.resize(x_size + extra_size + tmp_size);
  
  
  for(u32 i=0; i<x_size; ++i)
    {
    out[i] = x[i];
    }
  
  for(u32 i=0; i<extra_size; ++i)
    {
    out[x_size + i] = extra[i];
    }
  
  for(u32 i=0; i<tmp_size; ++i)
    {
    out[x_size + extra_size + i] = tmp[i];
    }
  
  return out;
  }



//! Safely rename a file.
//! Before renaming, test if we can write to the final file.
//! This should prevent:
//! (i)  overwriting files that have been write protected,
//! (ii) overwriting directories.
inline
void
diskio::safe_rename(const std::string& old_name, const std::string& new_name)
  {
  std::fstream f(new_name.c_str(), std::fstream::out | std::fstream::app);
  f.put(' ');
  
  const bool writing_problem = (f.good() == false);
  f.close();
  
  arma_warn( writing_problem, "trouble writing ", new_name );
  
  if(writing_problem == false)
    {
    std::remove(new_name.c_str());
    
    const int mv_result = std::rename(old_name.c_str(), new_name.c_str());
    arma_warn( (mv_result != 0), "trouble writing ", new_name );
    }
  
  }



//! Save a matrix as raw text (no header, human readable).
//! Matrices can be loaded in Matlab and Octave, as long as they don't have complex elements.
template<typename eT>
inline
void
diskio::save_raw_ascii(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out);
  
  if(f.is_open() == false)
    {
    arma_print("unable to write ", tmp_name);
    }
  else
    {
    diskio::save_raw_ascii(x, tmp_name, f);
    const bool writing_problem = (f.good() == false);
    
    f.flush();
    f.close();
    
    arma_warn(writing_problem, "trouble writing ", tmp_name);
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



//! Save a matrix as raw text (no header, human readable).
//! Matrices can be loaded in Matlab and Octave, as long as they don't have complex elements.
template<typename eT>
inline
void
diskio::save_raw_ascii(const Mat<eT>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  u32 cell_width;
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(8);
    cell_width = 16;
    }
  
  for(u32 row=0; row < x.n_rows; ++row)
    {
    for(u32 col=0; col < x.n_cols; ++col)
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
  }



//! Save a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_ascii(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  if(f.is_open() == false)
    {
    arma_debug_print("unable to write ", tmp_name);
    }
  else
    {
    diskio::save_arma_ascii(x, tmp_name, f);
    
    const bool writing_problem = (f.good() == false);
    
    f.flush();
    f.close();
    
    arma_warn( writing_problem, "trouble writing ", tmp_name );
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



//! Save a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void 
diskio::save_arma_ascii(const Mat<eT>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const ios::fmtflags orig_flags = f.flags();
  
  f << diskio::gen_txt_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << '\n';
  
  u32 cell_width;
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(8);
    cell_width = 16;
    }
    
  for(u32 row=0; row < x.n_rows; ++row)
    {
    for(u32 col=0; col < x.n_cols; ++col)
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
  
  f.flags(orig_flags);
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_binary(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_print("unable to write ", tmp_name);
    }
  else
    {  
    diskio::save_arma_binary(x, tmp_name, f);
    
    const bool writing_problem = (f.good() == false);
    
    f.flush();
    f.close();
    
    arma_warn( writing_problem, "trouble writing ", tmp_name );
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_binary(const Mat<eT>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();

  f << diskio::gen_bin_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << '\n';
  
  f.write(reinterpret_cast<const char*>(x.mem), x.n_elem*sizeof(eT));
  }



//! Save a matrix as a PGM greyscale image
template<typename eT>
inline
void
diskio::save_pgm_binary(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_print("unable to write ", tmp_name);
    }
  else
    {
    diskio::save_pgm_binary(x, tmp_name, f);
    
    const bool writing_problem = (f.good() == false);
    
    arma_warn(writing_problem, "trouble writing ", tmp_name );
    
    f.flush();
    f.close();
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



//
// TODO:
// add functionality to save the image in a normalised format,
// i.e. scaled so that every value falls in the [0,255] range.

//! Save a matrix as a PGM greyscale image
template<typename eT>
inline
void
diskio::save_pgm_binary(const Mat<eT>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f << "P5" << '\n';
  f << x.n_cols << ' ' << x.n_rows << '\n';
  f << 255 << '\n';
  
  const u32 n_elem = x.n_rows * x.n_cols;
  podarray<u8> tmp(n_elem);
  
  u32 i = 0;
  
  for(u32 row=0; row < x.n_rows; ++row)
    {
    for(u32 col=0; col < x.n_cols; ++col)
      {
      tmp[i] = u8( x(row,col) );  // TODO: add round() ?
      ++i;
      }
    }
  
  f.write(reinterpret_cast<const char*>(tmp.mem), n_elem);
  }



//! Save a matrix as a PGM greyscale image
template<typename T>
inline
void
diskio::save_pgm_binary(const Mat< std::complex<T> >& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const uchar_mat tmp = conv_to<uchar_mat>::from(x);
  diskio::save_pgm_binary(tmp, final_name);
  }



//! Save a matrix as a PGM greyscale image
template<typename T>
inline
void
diskio::save_pgm_binary(const Mat< std::complex<T> >& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const uchar_mat tmp = conv_to<uchar_mat>::from(x);
  diskio::save_pgm_binary(tmp, name, f);
  }



//! Load a matrix as raw text (no header, human readable).
//! Can read matrices saved as text in Matlab and Octave.
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
void
diskio::load_raw_ascii(Mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();

  std::fstream f;
  f.open(name.c_str(), std::fstream::in);
  
  if(f.is_open() == false)
    {
    x.reset();
    arma_extra_debug_print("unable to read ", name);
    }
  else
    {
    diskio::load_raw_ascii(x, name, f);
    f.close();
    }
  }

//! Load a matrix as raw text (no header, human readable).
//! Can read matrices saved as text in Matlab and Octave.
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
void
diskio::load_raw_ascii(Mat<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();

  bool load_okay = true;
  
  //std::fstream::pos_type start = f.tellg();
  
  //
  // work out the size
  
  u32 f_n_rows = 0;
  u32 f_n_cols = 0;
  
  bool f_n_cols_found = false;
  
  std::string line_string;
  std::string token;
  
  while( (f.good() == true) && (load_okay == true) )
    {
    std::getline(f, line_string);
    if(line_string.size() == 0)
      break;
    
    std::stringstream line_stream(line_string);
    
    u32 line_n_cols = 0;
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
        arma_print("inconsistent number of columns in ", name );
        load_okay = false;
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
  
    eT val;
    
    for(u32 row=0; row < x.n_rows; ++row)
      {
      for(u32 col=0; col < x.n_cols; ++col)
        {
        // f >> token;
        // x.at(row,col) = eT( strtod(token.c_str(), 0) );
        
        f >> val;
        x.at(row,col) = val;
        }
      }
    }
  
  if(f.good() == false)
    {
    arma_print("trouble reading ", name );
    load_okay = false; 
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! Load a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_ascii(Mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f(name.c_str());
  if(f.is_open() == false)
    {
    x.reset();
    arma_extra_debug_print("unable to read ", name);
    }
  else
    {
    diskio::load_arma_ascii(x, name, f);
    f.close();
    }
  }
  


//! Load a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_ascii(Mat<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  u32 f_n_rows;
  u32 f_n_cols;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  
  if(f_header == diskio::gen_txt_header(x))
    {
    x.set_size(f_n_rows, f_n_cols);
    
    for(u32 row=0; row < x.n_rows; ++row)
      {
      for(u32 col=0; col < x.n_cols; ++col)
        {
        f >> x.at(row,col);
        }
      }
    
    if(f.good() == false)
      {
      arma_print("trouble reading ", name);
      load_okay = false;
      }
    }
  else
    {
    arma_print("incorrect header in ", name );
    load_okay = false;
    } 
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! Load a matrix in binary format,
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_binary(Mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  if(f.is_open() == false)
    {
    x.reset();
    arma_extra_debug_print("unable to read ", name);
    }
  else
    {
    diskio::load_arma_binary(x, name, f);
    f.close();
    }
  }



template<typename eT>
inline
void
diskio::load_arma_binary(Mat<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  u32 f_n_rows;
  u32 f_n_cols;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  
  if(f_header == diskio::gen_bin_header(x))
    {
    //f.seekg(1, ios::cur);  // NOTE: this may not be portable, as on a Windows machine a newline could be two characters
    f.get();
    
    x.set_size(f_n_rows,f_n_cols);
    f.read( reinterpret_cast<char *>(x.memptr()), x.n_elem*sizeof(eT));
    
    if(f.good() == false)
      {
      arma_print("trouble reading ", name);
      load_okay = false;
      }
    }
  else
    {
    arma_print("incorrect header in ", name);
    load_okay = false;
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
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
void
diskio::load_pgm_binary(Mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    x.reset();
    }
  else
    {
    diskio::load_pgm_binary(x, name, f); // Do the actual load
    f.close();
    }
  }



//! Load a PGM greyscale image as a matrix
template<typename eT>
inline
void
diskio::load_pgm_binary(Mat<eT>& x, const std::string& name, std::istream& f)
  {
  bool load_okay = true;
  
  std::string f_header;
  f >> f_header;
  
  if(f_header == "P5")
    {
    u32 f_n_rows = 0;
    u32 f_n_cols = 0;
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
        const u32 n_elem = f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), n_elem);
        
        u32 i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        
        for(u32 row=0; row < f_n_rows; ++row)
          {
          for(u32 col=0; col < f_n_cols; ++col)
            {
            x.at(row,col) = eT(tmp[i]);
            ++i;
            }
          }
          
        }
      else
        {
        const u32 n_elem = f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), n_elem*2);
        
        u32 i = 0;
        
        for(u32 row=0; row < f_n_rows; ++row)
          {
          for(u32 col=0; col < f_n_cols; ++col)
            {
            x.at(row,col) = eT(tmp[i]);
            ++i;
            }
          }
        
        }
      
      }
    
    if(f.good() == false)
      {
      arma_print("trouble reading ", name);
      load_okay = false;
      }
    }
  else
    {
    arma_print("unsupported header in ", name);
    load_okay = false;
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! Load a PGM greyscale image as a matrix
template<typename T>
inline
void
diskio::load_pgm_binary(Mat< std::complex<T> >& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  uchar_mat tmp;
  diskio::load_pgm_binary(tmp, name);
  x = conv_to< Mat< std::complex<T> > >::from(tmp);
  }



//! Load a PGM greyscale image as a matrix
template<typename T>
inline
void
diskio::load_pgm_binary(Mat< std::complex<T> >& x, const std::string& name, std::istream& is)
  {
  arma_extra_debug_sigprint();
  
  uchar_mat tmp;
  diskio::load_pgm_binary(tmp, name, is);
  x = conv_to< Mat< std::complex<T> > >::from(tmp);
  }



//! Try to load a matrix by automatically determining its type
template<typename eT>
inline
void
diskio::load_auto_detect(Mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    x.reset();
    }
  else
    {
    diskio::load_auto_detect(x, name, f); // Do the actual load
    f.close();
    }
  }


//! Try to load a matrix by automatically determining its type
template<typename eT>
inline
void
diskio::load_auto_detect(Mat<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  static const std::string ARMA_MAT_TXT = "ARMA_MAT_TXT";
  static const std::string ARMA_MAT_BIN = "ARMA_MAT_BIN";
  static const std::string           P5 = "P5";
  
  podarray<char> raw_header(ARMA_MAT_TXT.length() + 1);
  
  std::streampos pos = f.tellg();
    
  f.read(raw_header.memptr(), ARMA_MAT_TXT.length());
  raw_header[ARMA_MAT_TXT.length()] = '\0';
  
  f.clear();
  f.seekg(pos);
  
  const std::string header = raw_header.mem;
  
  if(ARMA_MAT_TXT == header.substr(0,ARMA_MAT_TXT.length()))
    {
    load_arma_ascii(x, name, f);
    }
  else
  if(ARMA_MAT_BIN == header.substr(0,ARMA_MAT_BIN.length()))
    {
    load_arma_binary(x, name, f);
    }
  else
  if(P5 == header.substr(0,P5.length()))
    {
    load_pgm_binary(x, name, f);
    }
  else
    {
    load_raw_ascii(x, name, f);
    }
  }



// cubes



//! Save a cube as raw text (no header, human readable).
template<typename eT>
inline
void
diskio::save_raw_ascii(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out);
  
  if(f.is_open() == false)
    {
    arma_print("unable to write ", tmp_name);
    }
  else
    {
    save_raw_ascii(x, tmp_name, f);
    
    const bool writing_problem = (f.good() == false);
    
    arma_warn(writing_problem, "trouble writing ", tmp_name );
    
    f.flush();
    f.close();
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



//! Save a cube as raw text (no header, human readable).
template<typename eT>
inline
void
diskio::save_raw_ascii(const Cube<eT>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  u32 cell_width;
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(8);
    cell_width = 16;
    }
  
  for(u32 slice=0; slice < x.n_slices; ++slice)
    {
    for(u32 row=0; row < x.n_rows; ++row)
      {
      for(u32 col=0; col < x.n_cols; ++col)
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
  }



//! Save a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_ascii(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  if(f.is_open() == false)
    {
    arma_debug_print("unable to write ", tmp_name);
    }
  else
    {
    diskio::save_arma_ascii(x, tmp_name, f);
    
    const bool writing_problem = (f.good() == false);
    
    f.flush();
    f.close();
    
    arma_warn( writing_problem, "trouble writing ", tmp_name );
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



//! Save a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
void 
diskio::save_arma_ascii(const Cube<eT>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const ios::fmtflags orig_flags = f.flags();
  
  f << diskio::gen_txt_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << ' ' << x.n_slices << '\n';
  
  u32 cell_width;
  
  // TODO: need sane values for complex numbers
  
  if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
    {
    f.setf(ios::scientific);
    f.precision(8);
    cell_width = 16;
    }
    
  for(u32 slice=0; slice < x.n_slices; ++slice)
    {
    for(u32 row=0; row < x.n_rows; ++row)
      {
      for(u32 col=0; col < x.n_cols; ++col)
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
  
  f.flags(orig_flags);
  }



//! Save a cube in binary format,
//! with a header that stores the cube type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_binary(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_print("unable to write ", tmp_name);
    }
  else
    {
    diskio::save_arma_binary(x, tmp_name, f);
    
    const bool writing_problem = (f.good() == false);
    
    f.flush();
    f.close();
    
    arma_warn( writing_problem, "trouble writing ", tmp_name );
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



//! Save a cube in binary format,
//! with a header that stores the cube type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_binary(const Cube<eT>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  f << diskio::gen_bin_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << ' ' << x.n_slices << '\n';
  
  f.write(reinterpret_cast<const char*>(x.mem), x.n_elem*sizeof(eT));
  }



//! Load a cube as raw text (no header, human readable).
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
void
diskio::load_raw_ascii(Cube<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();

  Mat<eT> tmp;
  diskio::load_raw_ascii(tmp, name);
  
  x.set_size(tmp.n_rows, tmp.n_cols, 1);

  if(x.n_slices > 0)
    {
    x.slice(0) = tmp;
    }
  }



//! Load a cube as raw text (no header, human readable).
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
void
diskio::load_raw_ascii(Cube<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();

  Mat<eT> tmp;
  diskio::load_raw_ascii(tmp, name, f);
  
  x.set_size(tmp.n_rows, tmp.n_cols, 1);

  if(x.n_slices > 0)
    {
    x.slice(0) = tmp;
    }
  }



//! Load a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_ascii(Cube<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f(name.c_str());
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    }
  else
    {
    diskio::load_arma_ascii(x, name, f);
    f.close();
    }
  }
  


//! Load a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_ascii(Cube<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  u32 f_n_rows;
  u32 f_n_cols;
  u32 f_n_slices;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  f >> f_n_slices;
  
  if(f_header == diskio::gen_txt_header(x))
    {
    x.set_size(f_n_rows, f_n_cols, f_n_slices);

    for(u32 slice=0; slice < x.n_slices; ++slice)
      {
      for(u32 row=0; row < x.n_rows; ++row)
        {
        for(u32 col=0; col < x.n_cols; ++col)
          {
          f >> x.at(row,col,slice);
          }
        }
      }
    
    if(f.good() == false)
      {
      arma_print("trouble reading ", name);
      load_okay = false;
      }
    }
  else
    {
    arma_print("incorrect header in ", name );
    load_okay = false;
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! Load a cube in binary format,
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_binary(Cube<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    }
  else
    {
    diskio::load_arma_binary(x, name, f);
    f.close();
    }
  }



template<typename eT>
inline
void
diskio::load_arma_binary(Cube<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  u32 f_n_rows;
  u32 f_n_cols;
  u32 f_n_slices;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  f >> f_n_slices;
  
  if(f_header == diskio::gen_bin_header(x))
    {
    //f.seekg(1, ios::cur);  // NOTE: this may not be portable, as on a Windows machine a newline could be two characters
    f.get();
    
    x.set_size(f_n_rows, f_n_cols, f_n_slices);
    f.read( reinterpret_cast<char *>(x.memptr()), x.n_elem*sizeof(eT));
    
    if(f.good() == false)
      {
      arma_print("trouble reading ", name);
      load_okay = false;
      }
    }
  else
    {
    arma_print("incorrect header in ", name);
    load_okay = false;
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! Try to load a cube by automatically determining its type
template<typename eT>
inline
void
diskio::load_auto_detect(Cube<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    x.reset();
    }
  else
    {
    diskio::load_auto_detect(x, name, f); // Do the actual load
    f.close();
    }
  }



//! Try to load a cube by automatically determining its type
template<typename eT>
inline
void
diskio::load_auto_detect(Cube<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  static const std::string ARMA_CUB_TXT = "ARMA_CUB_TXT";
  static const std::string ARMA_CUB_BIN = "ARMA_CUB_BIN";
  static const std::string           P6 = "P6";
  
  podarray<char> raw_header(ARMA_CUB_TXT.length() + 1);
  
  std::streampos pos = f.tellg();
  
  f.read(raw_header.memptr(), ARMA_CUB_TXT.length());
  raw_header[ARMA_CUB_TXT.length()] = '\0';
  
  f.clear();
  f.seekg(pos);
  
  const std::string header = raw_header.mem;
  
  if(ARMA_CUB_TXT == header.substr(0, ARMA_CUB_TXT.length()))
    {
    load_arma_ascii(x, name, f);
    }
  else
  if(ARMA_CUB_BIN == header.substr(0, ARMA_CUB_BIN.length()))
    {
    load_arma_binary(x, name, f);
    }
  else
  if(P6 == header.substr(0,P6.length()))
    {
    load_ppm_binary(x, name, f);
    }
  else
    {
    load_raw_ascii(x, name, f);
    }
  }





// fields



template<typename T1>
inline
void
diskio::save_arma_binary(const field<T1>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  if(f.is_open() == false)
    {
    arma_print("couldn't write ", tmp_name);
    }
  else
    {
    diskio::save_arma_binary(x, tmp_name, f);

    const bool writing_problem = (f.good() == false);

    arma_warn(writing_problem, "trouble writing ", tmp_name );

    f.flush();
    f.close();

    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }

template<typename T1>
inline
void
diskio::save_arma_binary(const field<T1>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< (is_Mat<T1>::value == false) && (is_Cube<T1>::value == false) >::apply();
  
  f << "ARMA_FLD_BIN" << '\n';
  f << x.n_rows << '\n';
  f << x.n_cols << '\n';
  
  for(u32 i=0; i<x.n_elem; ++i)
    {
    diskio::save_arma_binary(x[i], name, f);
    }
  }



template<typename T1>
inline
void
diskio::load_arma_binary(field<T1>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f( name.c_str(), std::fstream::binary );
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    }
  else
    {
    diskio::load_arma_binary(x, name, f);
    f.close();
    }
  }



template<typename T1>
inline
void
diskio::load_arma_binary(field<T1>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< (is_Mat<T1>::value == false) && (is_Cube<T1>::value == false) >::apply();
  
  bool load_okay = true;
  
  std::string f_type;
  f >> f_type;
  
  if(f_type != "ARMA_FLD_BIN")
    {
    arma_print("unsupported field type in ", name);
    load_okay = false;
    }
  else
    {
    u32 f_n_rows;
    u32 f_n_cols;
  
    f >> f_n_rows;
    f >> f_n_cols;
    
    x.set_size(f_n_rows, f_n_cols);
    
    f.get();      
    
    for(u32 i=0; i<x.n_elem; ++i)
      {
      diskio::load_arma_binary(x[i], name, f);
      
      if(f.good() == false)
        {
        arma_print("trouble reading ", name);
        load_okay = false;
        break;
        }
      }
    }

  if(load_okay == false)
    {
    x.reset();
    }
  }



inline
void
diskio::save_std_string(const field<std::string>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  if(f.is_open() == false)
    {
    arma_print("couldn't write ", tmp_name);
    }
  else
    {
    diskio::save_std_string(x, tmp_name, f);
    
    const bool writing_problem = (f.good() == false);
    
    f.flush();
    f.close();
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



inline
void
diskio::save_std_string(const field<std::string>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  for(u32 row=0; row<x.n_rows; ++row)
  for(u32 col=0; col<x.n_cols; ++col)
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
  
  const bool writing_problem = (f.good() == false);
  
  arma_warn(writing_problem, "trouble writing ", name );
  }



inline
void
diskio::load_std_string(field<std::string>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f( name.c_str() );
  
  if(f.is_open() == false)
    {
    arma_print("unable to read ", name);
    }
  else
    {
    diskio::load_std_string(x, name, f);
    
    f.close();
    }
  }



inline
void
diskio::load_std_string(field<std::string>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  //
  // work out the size
  
  u32 f_n_rows = 0;
  u32 f_n_cols = 0;
  
  bool f_n_cols_found = false;
  
  std::string line_string;
  std::string token;
  
  while( (f.good() == true) && (load_okay == true) )
    {
    std::getline(f, line_string);
    if(line_string.size() == 0)
      break;
    
    std::stringstream line_stream(line_string);
    
    u32 line_n_cols = 0;
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
        arma_print("inconsistent number of columns in ", name );
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
  
    for(u32 row=0; row < x.n_rows; ++row)
      {
      for(u32 col=0; col < x.n_cols; ++col)
        {
        f >> x.at(row,col);
        }
      }
    }
  
  if(f.good() == false)
    {
    load_okay = false; 
    arma_print("trouble reading ", name );
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! Try to load a field by automatically determining its type
template<typename T1>
inline
void
diskio::load_auto_detect(field<T1>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    x.reset();
    }
  else
    {
    diskio::load_auto_detect(x, name, f); // Do the actual load
    f.close();
    }
  }



//! Try to load a field by automatically determining its type
template<typename T1>
inline
void
diskio::load_auto_detect(field<T1>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_Mat<T1>::value == false>::apply();
  
  static const std::string ARMA_FLD_BIN = "ARMA_FLD_BIN";
  static const std::string           P6 = "P6";
  
  podarray<char> raw_header(ARMA_FLD_BIN.length() + 1);
  
  std::streampos pos = f.tellg();
  
  f.read(raw_header.memptr(), ARMA_FLD_BIN.length());
  
  f.clear();
  f.seekg(pos);
  
  raw_header[ARMA_FLD_BIN.length()] = '\0';
  
  const std::string header = raw_header.mem;
  
  if(ARMA_FLD_BIN == header.substr(0,ARMA_FLD_BIN.length()))
    {
    load_arma_binary(x, name, f);
    }
  else
  if(P6 == header.substr(0,P6.length()))
    {
    load_ppm_binary(x, name, f);
    }
  else
    {
    arma_print("unsupported header in ", name);
    x.reset();
    }
  }



//
// handling of PPM images


template<typename eT>
inline
void
diskio::load_ppm_binary(Cube<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    }
  else
    {
    diskio::load_ppm_binary(x, name, f);
    f.close();
    }
  }



template<typename eT>
inline
void
diskio::load_ppm_binary(Cube<eT>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  f >> f_header;
  
  if(f_header == "P6")
    {
    u32 f_n_rows = 0;
    u32 f_n_cols = 0;
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
        const u32 n_elem = 3*f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), n_elem);
        
        u32 i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        
        for(u32 row=0; row < f_n_rows; ++row)
          {
          for(u32 col=0; col < f_n_cols; ++col)
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
        const u32 n_elem = 3*f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), 2*n_elem);
        
        u32 i = 0;
        
        for(u32 row=0; row < f_n_rows; ++row)
          {
          for(u32 col=0; col < f_n_cols; ++col)
            {
            x.at(row,col,0) = eT(tmp[i+0]);
            x.at(row,col,1) = eT(tmp[i+1]);
            x.at(row,col,2) = eT(tmp[i+2]);
            i+=3;
            }
          
          }
        
        }
      
      }
    
    if(f.good() == false)
      {
      arma_print("trouble reading ", name);
      load_okay = false;
      }
    
    }
  else
    {
    arma_print("unsupported header in ", name);
    load_okay = false;
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



template<typename eT>
inline
void
diskio::save_ppm_binary(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();

  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  if(f.is_open() == false)
    {
    arma_print("couldn't write ", tmp_name);
    }
  else
    {
    diskio::save_ppm_binary(x, tmp_name, f);
    
    const bool writing_problem = (f.good() == false);
    f.flush();
    f.close();
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }



template<typename eT>
inline
void
diskio::save_ppm_binary(const Cube<eT>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (x.n_slices != 3), "diskio::save_ppm_binary(): given cube must have exactly 3 slices" );
  
  const u32 n_elem = 3 * x.n_rows * x.n_cols;
  podarray<u8> tmp(n_elem);

  u32 i = 0;
  for(u32 row=0; row < x.n_rows; ++row)
    {
    for(u32 col=0; col < x.n_cols; ++col)
      {
      tmp[i+0] = u8( x.at(row,col,0) );
      tmp[i+1] = u8( x.at(row,col,1) );
      tmp[i+2] = u8( x.at(row,col,2) );
      
      i+=3;
      }
    }
  
  f << "P6" << '\n';
  f << x.n_cols << '\n';
  f << x.n_rows << '\n';
  f << 255 << '\n';

  f.write(reinterpret_cast<const char*>(tmp.mem), n_elem);
  
  const bool writing_problem = (f.good() == false);
  
  arma_warn(writing_problem, "trouble writing ", name );
  }



template<typename T1>
inline
void
diskio::load_ppm_binary(field<T1>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_extra_debug_print("unable to read ", name);
    }
  else
    {
    diskio::load_ppm_binary(x, name, f);
    f.close();
    }
  }  



template<typename T1>
inline
void
diskio::load_ppm_binary(field<T1>& x, const std::string& name, std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_Mat<T1>::value == false>::apply();
  typedef typename T1::elem_type eT;
  
  bool load_okay = true;
  
  std::string f_header;
  f >> f_header;
  
  if(f_header == "P6")
    {
    u32 f_n_rows = 0;
    u32 f_n_cols = 0;
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
        const u32 n_elem = 3*f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), n_elem);
        
        u32 i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        
        for(u32 row=0; row < f_n_rows; ++row)
          {
          for(u32 col=0; col < f_n_cols; ++col)
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
        const u32 n_elem = 3*f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), 2*n_elem);
        
        u32 i = 0;
        
        for(u32 row=0; row < f_n_rows; ++row)
          {
          for(u32 col=0; col < f_n_cols; ++col)
            {
            R.at(row,col) = eT(tmp[i+0]);
            G.at(row,col) = eT(tmp[i+1]);
            B.at(row,col) = eT(tmp[i+2]);
            i+=3;
            }
          
          }
        
        }
      
      }
    
    if(f.good() == false)
      {
      arma_print("trouble reading ", name);
      load_okay = false;
      }
    
    }
  else
    {
    arma_print("unsupported header in ", name);
    load_okay = false;
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
  
  }



template<typename T1>
inline
void
diskio::save_ppm_binary(const field<T1>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  if(f.is_open() == false)
    {
    arma_print("couldn't write ", tmp_name);
    }
  else
    {
    diskio::save_ppm_binary(x, tmp_name, f);
    const bool writing_problem = (f.good() == false);
    
    f.flush();
    f.close();
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  }
  


template<typename T1>
inline
void
diskio::save_ppm_binary(const field<T1>& x, const std::string& name, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_Mat<T1>::value == false>::apply();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (x.n_elem != 3), "diskio::save_ppm_binary(): given field must have exactly 3 matrices of equal size" );
  
  bool same_size = true;
  for(u32 i=1; i<3; ++i)
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

  const u32 n_elem = 3 * R.n_rows * R.n_cols;
  podarray<u8> tmp(n_elem);

  u32 i = 0;
  for(u32 row=0; row < R.n_rows; ++row)
    {
    for(u32 col=0; col < R.n_cols; ++col)
      {
      tmp[i+0] = u8( access::tmp_real( R.at(row,col) ) );
      tmp[i+1] = u8( access::tmp_real( G.at(row,col) ) );
      tmp[i+2] = u8( access::tmp_real( B.at(row,col) ) );
      
      i+=3;
      }
    }
  
  f.write(reinterpret_cast<const char*>(tmp.mem), n_elem);
  
  const bool writing_problem = (f.good() == false);
  
  arma_warn(writing_problem, "trouble writing ", name );
  }



//! @}

