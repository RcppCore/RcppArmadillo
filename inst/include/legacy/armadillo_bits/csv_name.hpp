// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup diskio
//! @{


namespace csv_opts
  {
  typedef unsigned int flag_type;
  
  struct opts
    {
    const flag_type flags;
    
    inline constexpr explicit opts(const flag_type in_flags);
    
    inline const opts operator+(const opts& rhs) const;
    };
  
  inline
  constexpr
  opts::opts(const flag_type in_flags)
    : flags(in_flags)
    {}
  
  inline
  const opts
  opts::operator+(const opts& rhs) const
    {
    const opts result( flags | rhs.flags );
    
    return result;
    }
  
  // The values below (eg. 1u << 0) are for internal Armadillo use only.
  // The values can change without notice.
  
  static constexpr flag_type flag_none        = flag_type(0      );
  static constexpr flag_type flag_trans       = flag_type(1u << 0);
  static constexpr flag_type flag_no_header   = flag_type(1u << 1);
  static constexpr flag_type flag_with_header = flag_type(1u << 2);
  static constexpr flag_type flag_semicolon   = flag_type(1u << 3);
  static constexpr flag_type flag_strict      = flag_type(1u << 4);
  
  struct opts_none        : public opts { inline constexpr opts_none()        : opts(flag_none       ) {} };
  struct opts_trans       : public opts { inline constexpr opts_trans()       : opts(flag_trans      ) {} };
  struct opts_no_header   : public opts { inline constexpr opts_no_header()   : opts(flag_no_header  ) {} };
  struct opts_with_header : public opts { inline constexpr opts_with_header() : opts(flag_with_header) {} };
  struct opts_semicolon   : public opts { inline constexpr opts_semicolon()   : opts(flag_semicolon  ) {} };
  struct opts_strict      : public opts { inline constexpr opts_strict()      : opts(flag_strict     ) {} };
  
  static constexpr opts_none        none;
  static constexpr opts_trans       trans;
  static constexpr opts_no_header   no_header;
  static constexpr opts_with_header with_header;
  static constexpr opts_semicolon   semicolon;
  static constexpr opts_strict      strict;
  }


struct csv_name
  {
  typedef field<std::string> header_type;
  
  const std::string    filename;
  const csv_opts::opts opts;
  
        header_type  header_junk;
  const header_type& header_ro;
        header_type& header_rw;
  
  inline
  csv_name(const std::string& in_filename)
    : filename (in_filename        )
    , opts     (csv_opts::no_header)
    , header_ro(header_junk        )
    , header_rw(header_junk        )
    {}
  
  inline
  csv_name(const std::string& in_filename, const csv_opts::opts& in_opts)
    : filename (in_filename                  )
    , opts     (csv_opts::no_header + in_opts)
    , header_ro(header_junk                  )
    , header_rw(header_junk                  )
    {}
  
  inline
  csv_name(const std::string& in_filename,       field<std::string>& in_header)
    : filename (in_filename          )
    , opts     (csv_opts::with_header)
    , header_ro(in_header            )
    , header_rw(in_header            )
    {}
  
  inline
  csv_name(const std::string& in_filename, const field<std::string>& in_header)
    : filename  (in_filename          )
    , opts      (csv_opts::with_header)
    , header_ro(in_header             )
    , header_rw(header_junk           )
    {}
  
  inline
  csv_name(const std::string& in_filename,       field<std::string>& in_header, const csv_opts::opts& in_opts)
    : filename  (in_filename                    )
    , opts      (csv_opts::with_header + in_opts)
    , header_ro(in_header                       )
    , header_rw(in_header                       )
    {}
  
  inline
  csv_name(const std::string& in_filename, const field<std::string>& in_header, const csv_opts::opts& in_opts)
    : filename  (in_filename                    )
    , opts      (csv_opts::with_header + in_opts)
    , header_ro(in_header                       )
    , header_rw(header_junk                     )
    {}
  };


//! @}
