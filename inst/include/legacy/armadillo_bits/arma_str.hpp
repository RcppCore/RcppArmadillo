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


//! \addtogroup arma_str
//! @{


namespace arma_str
  {
  class char_buffer
    {
    public:
    
    static constexpr uword n_chars_prealloc = 1024;
    
    char* mem     = nullptr;
    uword n_chars = 0;
    
    char local_mem[n_chars_prealloc];
    
    inline
    ~char_buffer()
      {
      if(n_chars > n_chars_prealloc)  { std::free(mem); }
      
      mem     = nullptr; 
      n_chars = 0;
      }
    
    inline
    char_buffer()
      {
      mem     = &(local_mem[0]);
      n_chars = n_chars_prealloc;
      
      if(n_chars > 0)  { mem[0] = char(0); }
      }
    
    inline
    void
    set_size(const uword new_n_chars)
      {
      if(n_chars > n_chars_prealloc)  { std::free(mem); }
      
      mem     = (new_n_chars <= n_chars_prealloc) ? &(local_mem[0])  : (char*)std::malloc(new_n_chars);
      n_chars = (new_n_chars <= n_chars_prealloc) ? n_chars_prealloc : new_n_chars;
      
      if(n_chars > 0)  { mem[0] = char(0); }
      }
    };
  
  
  class format
    {
    public:
    
    const std::string fmt;
    
    inline format(const char*        in_fmt) : fmt(in_fmt) { }
    inline format(const std::string& in_fmt) : fmt(in_fmt) { }
    
    private:
    format();
    };
  
  
  
  template<typename T1, typename T2>
  class basic_format
    {
    public:
    
    const T1& A;
    const T2& B;
    
    inline basic_format(const T1& in_A, const T2& in_B) : A(in_A) , B(in_B) { }
    
    private:
    basic_format();
    };
  
  
  
  template<typename T2>
  inline
  basic_format< format, T2 >
  operator% (const format& X, const T2& arg)
    {
    return basic_format< format, T2 >(X, arg);
    }
  
  
  
  template<typename T1, typename T2, typename T3>
  inline
  basic_format< basic_format<T1,T2>, T3 >
  operator% (const basic_format<T1,T2>& X, const T3& arg)
    {
    return basic_format< basic_format<T1,T2>, T3 >(X, arg);
    }
  
  
  
  template<typename T2>
  inline
  std::string
  str(const basic_format< format, T2>& X)
    {
    std::string out;
    char_buffer buf;
    
    bool status = false;
    
    while(status == false)
      {
      int required_size = (std::snprintf)(buf.mem, size_t(buf.n_chars), X.A.fmt.c_str(), X.B);
      
      if(required_size < 0)  { break; }
      
      if(uword(required_size) >= buf.n_chars)
        {
        if(buf.n_chars > char_buffer::n_chars_prealloc)  { break; }
        
        buf.set_size(1 + uword(required_size));
        }
      else
        {
        status = true;
        }
      
      if(status)  { out = buf.mem; }
      }
    
    return out;
    }
  
  
  
  template<typename T2, typename T3>
  inline
  std::string
  str(const basic_format< basic_format< format, T2>, T3>& X)
    {
    char_buffer buf;
    std::string out;
    
    bool status = false;
    
    while(status == false)
      {
      int required_size = (std::snprintf)(buf.mem, size_t(buf.n_chars), X.A.A.fmt.c_str(), X.A.B, X.B);
      
      if(required_size < 0)  { break; }
      
      if(uword(required_size) >= buf.n_chars)
        {
        if(buf.n_chars > char_buffer::n_chars_prealloc)  { break; }
        
        buf.set_size(1 + uword(required_size));
        }
      else
        {
        status = true;
        }
      
      if(status)  { out = buf.mem; }
      }
    
    return out;
    }
  
  
  
  template<typename T2, typename T3, typename T4>
  inline
  std::string
  str(const basic_format< basic_format< basic_format< format, T2>, T3>, T4>& X)
    {
    char_buffer buf;
    std::string out;
    
    bool status = false;
    
    while(status == false)
      {
      int required_size = (std::snprintf)(buf.mem, size_t(buf.n_chars), X.A.A.A.fmt.c_str(), X.A.A.B, X.A.B, X.B);
      
      if(required_size < 0)  { break; }
      
      if(uword(required_size) >= buf.n_chars)
        {
        if(buf.n_chars > char_buffer::n_chars_prealloc)  { break; }
        
        buf.set_size(1 + uword(required_size));
        }
      else
        {
        status = true;
        }
      
      if(status)  { out = buf.mem; }
      }
    
    return out;
    }
  
  
  
  template<typename T2, typename T3, typename T4, typename T5>
  inline
  std::string
  str(const basic_format< basic_format< basic_format< basic_format< format, T2>, T3>, T4>, T5>& X)
    {
    char_buffer buf;
    std::string out;
    
    bool status = false;
    
    while(status == false)
      {
      int required_size = (std::snprintf)(buf.mem, size_t(buf.n_chars), X.A.A.A.A.fmt.c_str(), X.A.A.A.B, X.A.A.B, X.A.B, X.B);
      
      if(required_size < 0)  { break; }
      
      if(uword(required_size) >= buf.n_chars)
        {
        if(buf.n_chars > char_buffer::n_chars_prealloc)  { break; }
        
        buf.set_size(1 + uword(required_size));
        }
      else
        {
        status = true;
        }
      
      if(status)  { out = buf.mem; }
      }
    
    return out;
    }
  
  
  
  template<typename T2, typename T3, typename T4, typename T5, typename T6>
  inline
  std::string
  str(const basic_format< basic_format< basic_format< basic_format< basic_format< format, T2>, T3>, T4>, T5>, T6>& X)
    {
    char_buffer buf;
    std::string out;
    
    bool status = false;
    
    while(status == false)
      {
      int required_size = (std::snprintf)(buf.mem, size_t(buf.n_chars), X.A.A.A.A.A.fmt.c_str(), X.A.A.A.A.B, X.A.A.A.B, X.A.A.B, X.A.B, X.B);
      
      if(required_size < 0)  { break; }
      
      if(uword(required_size) >= buf.n_chars)
        {
        if(buf.n_chars > char_buffer::n_chars_prealloc)  { break; }
        
        buf.set_size(1 + uword(required_size));
        }
      else
        {
        status = true;
        }
      
      if(status)  { out = buf.mem; }
      }
    
    return out;
    }
  
  
  
  template<typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  inline
  std::string
  str(const basic_format< basic_format< basic_format< basic_format< basic_format< basic_format< format, T2>, T3>, T4>, T5>, T6>, T7>& X)
    {
    char_buffer buf;
    std::string out;
    
    bool status = false;
    
    while(status == false)
      {
      int required_size = (std::snprintf)(buf.mem, size_t(buf.n_chars), X.A.A.A.A.A.A.fmt.c_str(), X.A.A.A.A.A.B, X.A.A.A.A.B, X.A.A.A.B, X.A.A.B, X.A.B, X.B);
      
      if(required_size < 0)  { break; }
      
      if(uword(required_size) >= buf.n_chars)
        {
        if(buf.n_chars > char_buffer::n_chars_prealloc)  { break; }
        
        buf.set_size(1 + uword(required_size));
        }
      else
        {
        status = true;
        }
      
      if(status)  { out = buf.mem; }
      }
    
    return out;
    }
  
  
  
  template<typename T1>
  struct format_metaprog
    {
    static constexpr uword depth = 0;
    
    inline
    static  
    const std::string&
    get_fmt(const T1& X)
      {
      return X.A;
      }
    };
  
  
  
  //template<>
  template<typename T1, typename T2>
  struct format_metaprog< basic_format<T1,T2> >
    {
    static constexpr uword depth = 1 + format_metaprog<T1>::depth;
    
    inline
    static
    const std::string&
    get_fmt(const T1& X)
      {
      return format_metaprog<T1>::get_fmt(X.A);
      }
    
    };
  
  
  
  template<typename T1, typename T2>
  inline
  std::string
  str(const basic_format<T1,T2>& X)
    {
    return format_metaprog< basic_format<T1,T2> >::get_fmt(X.A);
    }
  
  
  
  template<typename T1, typename T2>
  inline
  std::ostream&
  operator<< (std::ostream& o, const basic_format<T1,T2>& X)
    {
    o << str(X);
    return o;
    }
   
  
  template<typename T> struct string_only              { };
  template<>           struct string_only<std::string> { typedef std::string result; };
  
  template<typename T> struct char_only                { };
  template<>           struct char_only<char         > { typedef char        result; };
  
  template<typename T>
  struct basic_format_only { };
  
  template<typename T1, typename T2>
  struct basic_format_only< basic_format<T1, T2> > { typedef basic_format<T1,T2> result; };
  
  
  
  template<typename T1>
  inline
  static
  const T1&
  str_wrapper(const T1& x, const typename string_only<T1>::result* junk = nullptr)
    {
    arma_ignore(junk);
    
    return x;
    }
  
  
  
  template<typename T1>
  inline
  static
  const T1*
  str_wrapper(const T1* x, const typename char_only<T1>::result* junk = nullptr)
    {
    arma_ignore(junk);
    
    return x;
    }
  
  
  
  template<typename T1>
  inline
  static
  std::string
  str_wrapper(const T1& x, const typename basic_format_only<T1>::result* junk = nullptr)
    {
    arma_ignore(junk);
    
    return str(x);
    }
  
  }


//! @}
