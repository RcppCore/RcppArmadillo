// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppArmadillo.cpp: Rcpp/Armadillo glue
//
// Copyright (C)  2010 Dirk Eddelbuettel and Romain Francois
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

#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma ;
	
extern "C" SEXP RcppArmadillo_wrap(){
	
    // using the _ notation
    List matrices = make_list( 
    	_["Mat<int>"]          = eye<imat>( 3,3 ), 
    	_["Mat<double>"]       = eye<mat>( 3,3 ),
    	_["Mat<float>"]        = eye<fmat>( 3, 3 ), 
    	_["Mat<unsigned int>"] = eye<umat>( 3, 3 )
    	) ;
    
    // using the Named(.) = . notation
    List cols = make_list( 
    	Named( "Col<double>" ) = zeros<mat>(5,1), 
    	Named( "Col<float>" )  = zeros<fmat>(5,1)
    	) ; 
    
    // using the Named( . , . ) notation
    List rows = make_list( 
    	Named( "Row<double>", zeros<mat>(1,5) ),
    	Named( "Row<float>" , zeros<fmat>(1,5) )
    	) ;
	
    // creating an empty list and grow it on demand
    List fields ;
    field<int> f1( 2, 2 ) ;
    f1( 0, 0 ) = 0 ; 
    f1( 1, 0 ) = 1 ; 
    f1( 0, 1 ) = 2 ; 
    f1( 1, 1 ) = 3 ; 
    fields["field<int>"] = f1 ;
    field<std::string> f2(2,2) ;
    f2( 0, 0 ) = "a" ; 
    f2( 1, 0 ) = "b" ; 
    f2( 0, 1 ) = "c" ; 
    f2( 1, 1 ) = "d" ; 
    fields["field<std::string>"] = f2 ;
    field<colvec> f3(2,2) ;
    f3(0,0) = zeros<mat>(5,1) ;
    f3(1,0) = zeros<mat>(4,1) ;
    f3(0,1) = zeros<mat>(3,1) ;
    f3(1,1) = zeros<mat>(2,1) ;
    fields["field<colvec>"] = f3 ;
	
    List output = make_list( 
    	_["matrices : Mat<T>"]  = matrices, 
    	_["rows : Row<T>"]      = rows,
    	_["columns : Col<T>"]   = cols, 
    	_["fields  : field<T>"] = fields ) ;
    
    return output ;
	
}

extern "C" SEXP RcppArmadillo_as_Mat(SEXP input_){

    List input(input_) ;
    imat m1 = input[0] ; /* implicit as */
    mat  m2 = input[1] ; /* implicit as */
    umat m3 = input[0] ; /* implicit as */
    fmat m4 = input[1] ; /* implicit as */
	
    List res = make_list( 
    	accu( m1 ),
    	accu( m2 ),
    	accu( m3 ),
    	accu( m4 ) ) ;
    
    return res ;
}

extern "C" SEXP RcppArmadillo_as_Col( SEXP input_){
    List input(input_) ;
    icolvec m1 = input[0] ; /* implicit as */
    colvec  m2 = input[1] ; /* implicit as */
    ucolvec m3 = input[0] ; /* implicit as */
    fcolvec m4 = input[1] ; /* implicit as */
	
    List res = make_list( 
    	accu( m1 ),
    	accu( m2 ),
    	accu( m3 ),
    	accu( m4 ) ) ;
    
    return res ;
}

extern "C" SEXP RcppArmadillo_as_Row(SEXP input_){
    List input(input_) ;
    irowvec m1 = input[0] ; /* implicit as */
    rowvec  m2 = input[1] ; /* implicit as */
    urowvec m3 = input[0] ; /* implicit as */
    frowvec m4 = input[1] ; /* implicit as */
	
    List res = make_list( 
    	accu( m1 ),
    	accu( m2 ),
    	accu( m3 ),
    	accu( m4 ) ) ;
    return res ;
}

extern "C" SEXP RcppArmadillo_wrap_Glue(){
	
    mat m1 = eye<mat>( 3, 3 ) ;
    mat m2 = eye<mat>( 3, 3 ) ;
	
    List res ;
    res["mat+mat"] = m1 + m2 ;
    return res ;
}

extern "C" SEXP RcppArmadillo_wrap_Op(){

    mat m1 = eye<mat>( 3, 3 ) ;
	
    List res ;
    res["- mat"] = - m1 ;
    return res ;
}

extern "C" SEXP armadillo_version(SEXP single_){
    struct arma::arma_version av;
    bool single = Rcpp::as<bool>( single_) ;
    if( single ){
	return Rcpp::wrap( 10000*av.major + 100*av.minor + av.patch ) ;
    }
    Rcpp::IntegerVector version(3); 
    version = av.major, av.minor, av.patch ;
    Rcpp::CharacterVector names(3); 
    names = "major", "minor", "patch" ;
    version.names() = names ;
    version.attr("class" ) = "armadillo_version" ;
    return version ;
}

