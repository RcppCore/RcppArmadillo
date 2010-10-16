// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
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


//! class for saving and loading matrices and fields
class diskio
  {
  public:
  
  template<typename eT>
  struct is_supported_type
    {
    static const bool value = 
      (
      false
      || is_u8<eT>::value
      || is_s8<eT>::value
      || is_u16<eT>::value
      || is_s16<eT>::value
      || is_u32<eT>::value
      || is_s32<eT>::value
      || is_float<eT>::value
      || is_double<eT>::value
      || is_complex_float<eT>::value
      || is_complex_double<eT>::value
      );
    };
  
  
  template<typename eT> inline static std::string gen_txt_header(const Mat<eT>& x);
  template<typename eT> inline static std::string gen_bin_header(const Mat<eT>& x);

  template<typename eT> inline static std::string gen_txt_header(const Cube<eT>& x);
  template<typename eT> inline static std::string gen_bin_header(const Cube<eT>& x);

  
  inline static char conv_to_hex_char(const u8 x);
  inline static void conv_to_hex(char* out, const u8 x);

  inline static std::string gen_tmp_name(const std::string& x);
  
  inline static bool safe_rename(const std::string& old_name, const std::string& new_name);
  
  
  //
  // matrix saving
  
  template<typename eT> inline static bool save_raw_ascii  (const Mat<eT>&                x, const std::string& final_name);
  template<typename eT> inline static bool save_arma_ascii (const Mat<eT>&                x, const std::string& final_name);
  template<typename eT> inline static bool save_arma_binary(const Mat<eT>&                x, const std::string& final_name);
  template<typename eT> inline static bool save_pgm_binary (const Mat<eT>&                x, const std::string& final_name);
  template<typename  T> inline static bool save_pgm_binary (const Mat< std::complex<T> >& x, const std::string& final_name);
  
  template<typename eT> inline static bool save_raw_ascii  (const Mat<eT>&                x, std::ostream& f);
  template<typename eT> inline static bool save_arma_ascii (const Mat<eT>&                x, std::ostream& f);
  template<typename eT> inline static bool save_arma_binary(const Mat<eT>&                x, std::ostream& f);
  template<typename eT> inline static bool save_pgm_binary (const Mat<eT>&                x, std::ostream& f);
  template<typename  T> inline static bool save_pgm_binary (const Mat< std::complex<T> >& x, std::ostream& f);
  
  
  //
  // matrix loading
  
  template<typename eT> inline static bool load_raw_ascii  (Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_ascii (Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_pgm_binary (Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename  T> inline static bool load_pgm_binary (Mat< std::complex<T> >& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_auto_detect(Mat<eT>&                x, const std::string& name, std::string& err_msg);
  
  template<typename eT> inline static bool load_raw_ascii  (Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  template<typename eT> inline static bool load_arma_ascii (Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  template<typename eT> inline static bool load_pgm_binary (Mat<eT>&                x, std::istream& is, std::string& err_msg);
  template<typename  T> inline static bool load_pgm_binary (Mat< std::complex<T> >& x, std::istream& is, std::string& err_msg);
  template<typename eT> inline static bool load_auto_detect(Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  
  inline static void pnm_skip_comments(std::istream& f);
  
  
  //
  // cube saving
  
  template<typename eT> inline static bool save_raw_ascii  (const Cube<eT>& x, const std::string& name);
  template<typename eT> inline static bool save_arma_ascii (const Cube<eT>& x, const std::string& name);
  template<typename eT> inline static bool save_arma_binary(const Cube<eT>& x, const std::string& name);
  
  template<typename eT> inline static bool save_raw_ascii  (const Cube<eT>& x, std::ostream& f);
  template<typename eT> inline static bool save_arma_ascii (const Cube<eT>& x, std::ostream& f);
  template<typename eT> inline static bool save_arma_binary(const Cube<eT>& x, std::ostream& f);
  
  
  //
  // cube loading
  
  template<typename eT> inline static bool load_raw_ascii  (Cube<eT>& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_ascii (Cube<eT>& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(Cube<eT>& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_auto_detect(Cube<eT>& x, const std::string& name, std::string& err_msg);
  
  template<typename eT> inline static bool load_raw_ascii  (Cube<eT>& x, std::istream& f, std::string& err_msg);
  template<typename eT> inline static bool load_arma_ascii (Cube<eT>& x, std::istream& f, std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(Cube<eT>& x, std::istream& f, std::string& err_msg);
  template<typename eT> inline static bool load_auto_detect(Cube<eT>& x, std::istream& f, std::string& err_msg);
  
  
  //
  // field saving and loading
  
  template<typename T1> inline static bool save_arma_binary(const field<T1>& x, const std::string&  name);
  template<typename T1> inline static bool save_arma_binary(const field<T1>& x,       std::ostream& f);
  
  template<typename T1> inline static bool load_arma_binary(      field<T1>& x, const std::string&  name, std::string& err_msg);
  template<typename T1> inline static bool load_arma_binary(      field<T1>& x,       std::istream& f,    std::string& err_msg);
  
  template<typename T1> inline static bool load_auto_detect(      field<T1>& x, const std::string&  name, std::string& err_msg);
  template<typename T1> inline static bool load_auto_detect(      field<T1>& x,       std::istream& f,    std::string& err_msg);
  
  inline static bool save_std_string(const field<std::string>& x, const std::string&  name);
  inline static bool save_std_string(const field<std::string>& x,       std::ostream& f);
  
  inline static bool load_std_string(      field<std::string>& x, const std::string&  name, std::string& err_msg);
  inline static bool load_std_string(      field<std::string>& x,       std::istream& f,    std::string& err_msg);
  


  //
  // handling of PPM images by cubes

  template<typename T1> inline static bool save_ppm_binary(const Cube<T1>& x, const std::string&  final_name);
  template<typename T1> inline static bool save_ppm_binary(const Cube<T1>& x,       std::ostream& f);
  
  template<typename T1> inline static bool load_ppm_binary(      Cube<T1>& x, const std::string&  final_name, std::string& err_msg);
  template<typename T1> inline static bool load_ppm_binary(      Cube<T1>& x,       std::istream& f,          std::string& err_msg);


  //
  // handling of PPM images by fields

  template<typename T1> inline static bool save_ppm_binary(const field<T1>& x, const std::string&  final_name);
  template<typename T1> inline static bool save_ppm_binary(const field<T1>& x,       std::ostream& f);
  
  template<typename T1> inline static bool load_ppm_binary(      field<T1>& x, const std::string&  final_name, std::string& err_msg);
  template<typename T1> inline static bool load_ppm_binary(      field<T1>& x,       std::istream& f,          std::string& err_msg);
  


  };



//! @}
