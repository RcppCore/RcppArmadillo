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

extern "C" SEXP RcppArmadillo_wrap(){
	
    // using the Argument(.) = . notation
    List cols = List::create( 
    	Argument( "Col<double>" ) = arma::zeros<arma::mat>(5,1), 
    	Argument( "Col<float>" )  = arma::zeros<arma::fmat>(5,1)
    	) ; 
    List rows = List::create( 
    	Argument( "Row<double>") = arma::zeros<arma::mat>(1,5),
    	Argument( "Row<float>" ) = arma::zeros<arma::fmat>(1,5)
    	) ;
	
    // using the _[.] = . notation
    List matrices = List::create( 
    	_["Mat<int>"]          = arma::eye<arma::imat>( 3,3 ), 
    	_["Mat<double>"]       = arma::eye<arma::mat>( 3,3 ),
    	_["Mat<float>"]        = arma::eye<arma::fmat>( 3, 3 ), 
    	_["Mat<unsigned int>"] = arma::eye<arma::umat>( 3, 3 )
    	) ;
    
    // creating an empty list and grow it on demand
    List fields ;
    arma::field<int> f1( 2, 2 ) ;
    f1( 0, 0 ) = 0 ; 
    f1( 1, 0 ) = 1 ; 
    f1( 0, 1 ) = 2 ; 
    f1( 1, 1 ) = 3 ; 
    fields["field<int>"] = f1 ;
    
    arma::field<std::string> f2(2,2) ;
    f2( 0, 0 ) = "a" ; 
    f2( 1, 0 ) = "b" ; 
    f2( 0, 1 ) = "c" ; 
    f2( 1, 1 ) = "d" ; 
    fields["field<std::string>"] = f2 ;
    
    arma::field<arma::colvec> f3(2,2) ;
    f3(0,0) = arma::zeros<arma::mat>(5,1) ;
    f3(1,0) = arma::zeros<arma::mat>(4,1) ;
    f3(0,1) = arma::zeros<arma::mat>(3,1) ;
    f3(1,1) = arma::zeros<arma::mat>(2,1) ;
    fields["field<colvec>"] = f3 ;
	
    List output = List::create( 
    	_["matrices : Mat<T>"]  = matrices, 
    	_["rows : Row<T>"]      = rows,
    	_["columns : Col<T>"]   = cols, 
    	_["fields  : field<T>"] = fields ) ;
    
    return output ;
	
}

extern "C" SEXP RcppArmadillo_as_Mat(SEXP input_){

    List input(input_) ;
    arma::imat m1 = input[0] ; /* implicit as */
    arma::mat  m2 = input[1] ; /* implicit as */
    arma::umat m3 = input[0] ; /* implicit as */
    arma::fmat m4 = input[1] ; /* implicit as */
	
    List res = List::create( 
    	arma::accu( m1 ),
    	arma::accu( m2 ),
    	arma::accu( m3 ),
    	arma::accu( m4 ) ) ;
    
    return res ;
}

extern "C" SEXP RcppArmadillo_as_Col( SEXP input_){
    List input(input_) ;
    arma::icolvec m1 = input[0] ; /* implicit as */
    arma::colvec  m2 = input[1] ; /* implicit as */
    arma::ucolvec m3 = input[0] ; /* implicit as */
    arma::fcolvec m4 = input[1] ; /* implicit as */
	
    List res = List::create( 
    	arma::accu( m1 ),
    	arma::accu( m2 ),
    	arma::accu( m3 ),
    	arma::accu( m4 ) ) ;
    
    return res ;
}

extern "C" SEXP RcppArmadillo_as_Row(SEXP input_){
    List input(input_) ;
    arma::irowvec m1 = input[0] ; /* implicit as */
    arma::rowvec  m2 = input[1] ; /* implicit as */
    arma::urowvec m3 = input[0] ; /* implicit as */
    arma::frowvec m4 = input[1] ; /* implicit as */
	
    List res = List::create( 
    	arma::accu( m1 ),
    	arma::accu( m2 ),
    	arma::accu( m3 ),
    	arma::accu( m4 ) ) ;
    return res ;
}

extern "C" SEXP RcppArmadillo_wrap_Glue(){
	
    arma::mat m1 = arma::eye<arma::mat>( 3, 3 ) ;
    arma::mat m2 = arma::eye<arma::mat>( 3, 3 ) ;
	                     
    List res ;                                  
    res["mat+mat"] = m1 + m2 ;
    return res ;
}

extern "C" SEXP RcppArmadillo_wrap_Op(){

	arma::mat m1 = arma::eye<arma::mat>( 3, 3 ) ;
	
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
    IntegerVector version = IntegerVector::create( 
    	_["major"] = av.major, 
    	_["minor"] = av.minor, 
    	_["patch"] = av.patch
    	) ;
    return version ;
}

