// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppArmadilloForward.h: Rcpp/Armadillo glue
//
// Copyright (C)  2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
//
// This file is part of RcppArmadillo.
//
// RcppArmadillo is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppArmadillo is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RcppArmadillo__RcppArmadilloForward__h
#define RcppArmadillo__RcppArmadilloForward__h

#include <RcppCommon.h>
#include <armadillo>

/* forward declarations */
namespace Rcpp {
    /* support for wrap */
    template <typename T> SEXP wrap ( const arma::Mat<T>& ) ;
    template <typename T> SEXP wrap ( const arma::Row<T>& ) ;
    template <typename T> SEXP wrap ( const arma::Col<T>& ) ;
    template <typename T> SEXP wrap ( const arma::field<T>& ) ;
    template <typename T> SEXP wrap ( const arma::Cube<T>& ) ;
    
    template <typename T1, typename T2, typename glue_type> 
    SEXP wrap(const arma::Glue<T1, T2, glue_type>& X ) ;
    
    template <typename T1, typename op_type>
    SEXP wrap(const arma::Op<T1, op_type>& X ) ;
    
    template <typename T1, typename T2, typename glue_type> 
    SEXP wrap(const arma::eGlue<T1, T2, glue_type>& X ) ;
    
    template <typename T1, typename op_type>
    SEXP wrap(const arma::eOp<T1, op_type>& X ) ;
    
    template <typename T1, typename op_type>
    SEXP wrap(const arma::OpCube<T1,op_type>& X ) ;
    
    template <typename T1, typename T2, typename glue_type>
    SEXP wrap(const arma::GlueCube<T1,T2,glue_type>& X ) ;
    
    template <typename T1, typename op_type>
    SEXP wrap(const arma::eOpCube<T1,op_type>& X ) ;
    
    template <typename T1, typename T2, typename glue_type>
    SEXP wrap(const arma::eGlueCube<T1,T2,glue_type>& X ) ;

    template<typename out_eT, typename T1, typename op_type>
    SEXP wrap( const arma::mtOp<out_eT,T1,op_type>& X ) ;

    template<typename out_eT, typename T1, typename T2, typename glue_type>
    SEXP wrap( const arma::mtGlue<out_eT,T1,T2,glue_type>& X );
    
    namespace traits {

	/* support for as */
	template <typename T> class Exporter< arma::Mat<T> > ;
	template <typename T> class Exporter< arma::Row<T> > ;
	template <typename T> class Exporter< arma::Col<T> > ;
// template <typename T> class Exporter< arma::field<T> > ;
// template <typename T> class Exporter< arma::Cube<T> > ;

    } // namespace traits 

}

#endif
