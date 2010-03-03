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
	
	List output(4); 
	List matrices(4) ;
	List rows(2) ;
	List cols(2) ;
	
	imat x1 = eye<imat>( 3,3 ) ;
	// x1.swap_rows(0,1)  ;
	matrices[0] = x1 ;
	
	mat x2 = eye<mat>( 3,3 ) ;
	// x2.swap_rows(0,1)  ;
	matrices[1] = x2 ;
	
	fmat x3 = eye<fmat>( 3, 3 );
	// x3.swap_rows(0,1)  ;
	matrices[2] = x3 ;
	
	umat x4 = eye<umat>( 3, 3 );
	// x4.swap_rows(0,1)  ;
	matrices[3] = x4 ;
	
	colvec r1 = zeros<mat>(5,1);
	cols[0] = r1  ;
	
	fcolvec r2 = zeros<fmat>(5,1);
	cols[1] = r2  ;
	
	rowvec c1 = zeros<mat>(1,5);
	rows[0] = c1  ;
	
	frowvec c2 = zeros<fmat>(1,5);
	rows[1] = c2  ;
	
	std::vector<std::string> names(4) ;
	names[0] = "Mat<int>" ;
	names[1] = "Mat<double>" ;
	names[2] = "Mat<float>" ;
	names[3] = "Mat<unsigned int>" ;
	matrices.names() = names ;
	
	names.resize(2) ;
	names[0] = "Row<double>" ;
	names[1] = "Row<float>" ;
	rows.names() = names ;
	
	names[0] = "Col<double>" ;
	names[1] = "Col<float>" ;
	cols.names() = names ;
	
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
	
	output[0] = matrices;
	output[1] = rows ;
	output[2] = cols ;
	output[3] = fields ;
	
	names[0] = "matrices : Mat<T>" ; 
	names[1] = "rows : Row<T>" ;
	names.push_back( "columns : Col<T>" );
	names.push_back( "fields  : field<T>" );
	output.names() = names ;
	
	return output ;
	
}

extern "C" SEXP RcppArmadillo_as_Mat(SEXP input_){

	List input(input_) ;
	arma::imat m1 = input[0] ; /* implicit as */
	arma::mat  m2 = input[1] ; /* implicit as */
	arma::umat m3 = input[0] ; /* implicit as */
	arma::fmat m4 = input[1] ; /* implicit as */
	
	List res(4) ;
	res[0] = arma::accu( m1 ) ;
	res[1] = arma::accu( m2 ) ;
	res[2] = arma::accu( m3 ) ;
	res[3] = arma::accu( m4 ) ;
	
	return res ;
}

extern "C" SEXP RcppArmadillo_as_Col( SEXP input_){
	List input(input_) ;
	arma::icolvec m1 = input[0] ; /* implicit as */
	arma::colvec  m2 = input[1] ; /* implicit as */
	arma::ucolvec m3 = input[0] ; /* implicit as */
	arma::fcolvec m4 = input[1] ; /* implicit as */
	
	List res(4) ;
	res[0] = arma::accu( m1 ) ;
	res[1] = arma::accu( m2 ) ;
	res[2] = arma::accu( m3 ) ;
	res[3] = arma::accu( m4 ) ;
	
	return res ;
}

extern "C" SEXP RcppArmadillo_as_Row(SEXP input_){
	List input(input_) ;
	arma::irowvec m1 = input[0] ; /* implicit as */
	arma::rowvec  m2 = input[1] ; /* implicit as */
	arma::urowvec m3 = input[0] ; /* implicit as */
	arma::frowvec m4 = input[1] ; /* implicit as */
	
	List res(4) ;
	res[0] = arma::accu( m1 ) ;
	res[1] = arma::accu( m2 ) ;
	res[2] = arma::accu( m3 ) ;
	res[3] = arma::accu( m4 ) ;
	
	return res ;
	
}

extern "C" SEXP RcppArmadillo_wrap_Glue(){
	
	arma::mat m1 = eye<mat>( 3, 3 ) ;
	arma::mat m2 = eye<mat>( 3, 3 ) ;
	
	List res ;
	res["mat+mat"] = m1 + m2 ;
	return res ;
}

extern "C" SEXP RcppArmadillo_wrap_Op(){
	arma::mat m1 = eye<mat>( 3, 3 ) ;
	
	List res ;
	res["mat+mat"] = - m1 ;
	return res ;
	
}

extern "C" SEXP fastLm(SEXP sy, SEXP sX) {

    Rcpp::NumericVector yr(sy);
    Rcpp::NumericMatrix Xr(sX);
    std::vector<int> dims = Xr.attr("dim") ;
    int n = dims[0], k = dims[1];

#if ARMA_VERSION_GE_090
    arma::mat X(Xr.begin(), n, k, false);   // reuses memory and avoids extra copy
#else
    arma::mat X(Xr.begin(), n, k);
#endif

    arma::colvec y(yr.begin(), yr.size());

    arma::colvec coef = solve(X, y);            // fit model y ~ X

    arma::colvec resid = y - X*coef; 

    double sig2 = SCALAR( trans(resid)*resid/(n-k) );
    arma::mat covmat = sig2 * arma::inv(arma::trans(X)*X);
    arma::colvec stderrest = sqrt(covmat.diag());

    Rcpp::Pairlist res(Rcpp::Named( "coef", coef),
                       Rcpp::Named( "stderr", stderrest));
    return res;
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

