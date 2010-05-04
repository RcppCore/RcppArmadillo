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


//! \addtogroup format_wrap
//! @{


//! \namespace arma_boost namespace for functions and classes which partially emulate Boost functionality 
namespace arma_boost
  {
  
  #if defined(ARMA_USE_BOOST_FORMAT)

    using boost::format;
    using boost::basic_format;
    using boost::str;

  #else
  
    #if defined(ARMA_HAVE_STD_SNPRINTF)

      #define arma_snprintf std::snprintf

    #else

      // better-than-nothing emulation of C99 snprintf(),
      // with correct return value and null-terminated output string.
      // note that _snprintf() provided by MS is not a good substitute for snprintf()

      inline
      int
      arma_snprintf(char* out, size_t size, const char* fmt, ...)
        {
        size_t i;
        
        for(i=0; i<size; ++i)
          {
          out[i] = fmt[i];
          if(fmt[i] == char(0))
            break;
          }
        
        if(size > 0)
          out[size-1] = char(0);
        
        return int(i);
        }

    #endif
    
    class format
      {
      public:
    
      format(const char* in_fmt)
        : A(in_fmt)
        {
        }
    
      format(const std::string& in_fmt)
        : A(in_fmt)
        {
        }
    
    
      const std::string A;
    
      private:
      format();
      };
    
    
    
    template<typename T1, typename T2>
    class basic_format
      {
      public:
    
      basic_format(const T1& in_A, const T2& in_B)
        : A(in_A)
        , B(in_B)
        {
        }
    
      const T1& A;
      const T2& B;
    
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
      char  local_buffer[1024];
      char* buffer = local_buffer;
      
      int buffer_size   = 1024;
      int required_size = buffer_size;
   
      bool using_local_buffer = true;
      
      std::string out;
      
      do
        {
        if(using_local_buffer == false)
          {
          buffer = new char[buffer_size];
          }
        
        required_size = arma_snprintf(buffer, buffer_size, X.A.A.c_str(), X.B);
        
        if(required_size < buffer_size)
          {
          if(required_size > 0)
            {
            out = buffer;
            }
          }
        else
          {
          buffer_size *= 2;
          }
        
        if(using_local_buffer == true)
          {
          using_local_buffer = false;
          }
        else
          {
          delete[] buffer;
          }
        
        } while( (required_size >= buffer_size) );

      return out;
      }
    
    
    
    template<typename T2, typename T3>
    inline
    std::string
    str(const basic_format< basic_format< format, T2>, T3>& X)
      {
      char  local_buffer[1024];
      char* buffer = local_buffer;
      
      int buffer_size   = 1024;
      int required_size = buffer_size;
   
      bool using_local_buffer = true;
      
      std::string out;
      
      do
        {
        if(using_local_buffer == false)
          {
          buffer = new char[buffer_size];
          }
        
        required_size = arma_snprintf(buffer, buffer_size, X.A.A.A.c_str(), X.A.B, X.B);
        
        if(required_size < buffer_size)
          {
          if(required_size > 0)
            {
            out = buffer;
            }
          }
        else
          {
          buffer_size *= 2;
          }
        
        if(using_local_buffer == true)
          {
          using_local_buffer = false;
          }
        else
          {
          delete[] buffer;
          }
        
        } while( (required_size >= buffer_size) );

      return out;
      }
    
    
    
    template<typename T2, typename T3, typename T4>
    inline
    std::string
    str(const basic_format< basic_format< basic_format< format, T2>, T3>, T4>& X)
      {
      char  local_buffer[1024];
      char* buffer = local_buffer;
      
      int buffer_size   = 1024;
      int required_size = buffer_size;
   
      bool using_local_buffer = true;
      
      std::string out;
      
      do
        {
        if(using_local_buffer == false)
          {
          buffer = new char[buffer_size];
          }
        
        required_size = arma_snprintf(buffer, buffer_size, X.A.A.A.A.c_str(), X.A.A.B, X.A.B, X.B);
        
        if(required_size < buffer_size)
          {
          if(required_size > 0)
            {
            out = buffer;
            }
          }
        else
          {
          buffer_size *= 2;
          }
        
        if(using_local_buffer == true)
          {
          using_local_buffer = false;
          }
        else
          {
          delete[] buffer;
          }
        
        } while( (required_size >= buffer_size) );

      return out;
      }
    
    
    
    template<typename T2, typename T3, typename T4, typename T5>
    inline
    std::string
    str(const basic_format< basic_format< basic_format< basic_format< format, T2>, T3>, T4>, T5>& X)
      {
      char  local_buffer[1024];
      char* buffer = local_buffer;
      
      int buffer_size   = 1024;
      int required_size = buffer_size;
   
      bool using_local_buffer = true;
      
      std::string out;
      
      do
        {
        if(using_local_buffer == false)
          {
          buffer = new char[buffer_size];
          }
        
        required_size = arma_snprintf(buffer, buffer_size, X.A.A.A.A.A.c_str(), X.A.A.A.B, X.A.A.B, X.A.B, X.B);
        
        if(required_size < buffer_size)
          {
          if(required_size > 0)
            {
            out = buffer;
            }
          }
        else
          {
          buffer_size *= 2;
          }
        
        if(using_local_buffer == true)
          {
          using_local_buffer = false;
          }
        else
          {
          delete[] buffer;
          }
        
        } while( (required_size >= buffer_size) );

      return out;
      }
    
    
    
    template<typename T2, typename T3, typename T4, typename T5, typename T6>
    inline
    std::string
    str(const basic_format< basic_format< basic_format< basic_format< basic_format< format, T2>, T3>, T4>, T5>, T6>& X)
      {
      char  local_buffer[1024];
      char* buffer = local_buffer;
      
      int buffer_size   = 1024;
      int required_size = buffer_size;
   
      bool using_local_buffer = true;
      
      std::string out;
      
      do
        {
        if(using_local_buffer == false)
          {
          buffer = new char[buffer_size];
          }
        
        required_size = arma_snprintf(buffer, buffer_size, X.A.A.A.A.A.A.c_str(), X.A.A.A.A.B, X.A.A.A.B, X.A.A.B, X.A.B, X.B);
        
        if(required_size < buffer_size)
          {
          if(required_size > 0)
            {
            out = buffer;
            }
          }
        else
          {
          buffer_size *= 2;
          }
        
        if(using_local_buffer == true)
          {
          using_local_buffer = false;
          }
        else
          {
          delete[] buffer;
          }
        
        } while( (required_size >= buffer_size) );

      return out;
      }
    
    
    
    template<typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
    inline
    std::string
    str(const basic_format< basic_format< basic_format< basic_format< basic_format< basic_format< format, T2>, T3>, T4>, T5>, T6>, T7>& X)
      {
      char  local_buffer[1024];
      char* buffer = local_buffer;
      
      int buffer_size   = 1024;
      int required_size = buffer_size;
   
      bool using_local_buffer = true;
      
      std::string out;
      
      do
        {
        if(using_local_buffer == false)
          {
          buffer = new char[buffer_size];
          }
        
        required_size = arma_snprintf(buffer, buffer_size, X.A.A.A.A.A.A.A.c_str(), X.A.A.A.A.A.B, X.A.A.A.A.B, X.A.A.A.B, X.A.A.B, X.A.B, X.B);
        
        if(required_size < buffer_size)
          {
          if(required_size > 0)
            {
            out = buffer;
            }
          }
        else
          {
          buffer_size *= 2;
          }
        
        if(using_local_buffer == true)
          {
          using_local_buffer = false;
          }
        else
          {
          delete[] buffer;
          }
        
        } while( (required_size >= buffer_size) );

      return out;
      }
    
    
    
    template<typename T1>
    struct format_metaprog
      {
      static const u32 depth = 0;
    
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
      static const u32 depth = 1 + format_metaprog<T1>::depth;
    
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
        
        
  #endif
  
  }

//! @}
