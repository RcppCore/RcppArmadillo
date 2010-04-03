// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_shuffle
//! @{



template<typename T1>
inline
void
op_shuffle::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_shuffle>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  const u32 dim = in.aux_u32_a;
  const u32 N   = (dim == 0) ? X.n_rows : X.n_cols;
  
  // see "fn_sort_index.hpp" for the definition of "arma_sort_index_packet_ascend"
  // and the associated "operator<"
  std::vector< arma_sort_index_packet_ascend<int,u32> > packet_vec(N);
  
  for(u32 i=0; i<N; ++i)
    {
    packet_vec[i].val   = std::rand();
    packet_vec[i].index = i;
    }
  
  std::sort( packet_vec.begin(), packet_vec.end() );
  
  if(X.is_vec() == false)
    {
    if(&out != &X)
      {
      arma_extra_debug_print("op_shuffle::apply(): matrix");
      
      out.copy_size(X);
      
      if(dim == 0)
        {
        for(u32 i=0; i<N; ++i)
          {
          out.row(i) = X.row(packet_vec[i].index);
          }
        }
      else
        {
        for(u32 i=0; i<N; ++i)
          {
          out.col(i) = X.col(packet_vec[i].index);
          }
        }
      }
    else  // in-place shuffle
      {
      arma_extra_debug_print("op_shuffle::apply(): in-place matrix");
      
      // reuse the val member variable of packet_vec
      // to indicate whether a particular row or column
      // has already been shuffled
      
      for(u32 i=0; i<N; ++i)
        {
        packet_vec[i].val = 0;
        }
        
      if(dim == 0)
        {
        for(u32 i=0; i<N; ++i)
          {
          if(packet_vec[i].val == 0)
            {
            const u32 j = packet_vec[i].index;
            
            out.swap_rows(i, j);
            
            packet_vec[j].val = 1;
            }
          }
        }
      else
        {
        for(u32 i=0; i<N; ++i)
          {
          if(packet_vec[i].val == 0)
            {
            const u32 j = packet_vec[i].index;
            
            out.swap_cols(i, j);
            
            packet_vec[j].val = 1;
            }
          }
        }
      }
    }
  else  // we're dealing with a vector
    {
    if(&out != &X)
      {
      arma_extra_debug_print("op_shuffle::apply(): vector");
      
      out.copy_size(X);
      
      if(dim == 0)
        {
        if(X.n_rows > 1)  // i.e. column vector
          {
          for(u32 i=0; i<N; ++i)
            {
            out[i] = X[ packet_vec[i].index ];
            }
          }
        else
          {
          out = X;
          }
        }
      else
        {
        if(X.n_cols > 1)  // i.e. row vector
          {
          for(u32 i=0; i<N; ++i)
            {
            out[i] = X[ packet_vec[i].index ];
            }
          }
        else
          {
          out = X;
          }
        }
      }
    else  // in-place shuffle
      {
      arma_extra_debug_print("op_shuffle::apply(): in-place vector");
      
      // reuse the val member variable of packet_vec
      // to indicate whether a particular row or column
      // has already been shuffled
      
      for(u32 i=0; i<N; ++i)
        {
        packet_vec[i].val = 0;
        }
        
      if(dim == 0)
        {
        if(X.n_rows > 1)  // i.e. column vector
          {
          for(u32 i=0; i<N; ++i)
            {
            if(packet_vec[i].val == 0)
              {
              const u32 j = packet_vec[i].index;
              
              std::swap(out[i], out[j]);
              
              packet_vec[j].val = 1;
              }
            }
          }
        }
      else
        {
        if(X.n_cols > 1)  // i.e. row vector
          {
          for(u32 i=0; i<N; ++i)
            {
            if(packet_vec[i].val == 0)
              {
              const u32 j = packet_vec[i].index;
              
              std::swap(out[i], out[j]);
              
              packet_vec[j].val = 1;
              }
            }
          }
        }
      }
    }
  
  }


//! @}
