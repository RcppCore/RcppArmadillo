#!/usr/bin/r -t
#
# Copyright (C) 2010	Dirk Eddelbuettel, Romain Francois and Douglas Bates
#
# This file is part of RcppArmadillo.
#
# RcppArmadillo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# RcppArmadillo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

test.wrap.R <- function(){
	
	fx <- cxxfunction( , '
	
    // using the Named(.) = . notation
    List cols = List::create( 
    	Named( "Col<double>" ) = arma::zeros<arma::mat>(5,1), 
    	Named( "Col<float>" )  = arma::zeros<arma::fmat>(5,1)
    	) ;
    
    // using the Named(., .)  notation
    List rows = List::create( 
    	Named( "Row<double>",  arma::zeros<arma::mat>(1,5)  ),
    	Named( "Row<float>" ,  arma::zeros<arma::fmat>(1,5) )
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
	' , plugin = "RcppArmadillo" )
	
    res <- fx()

    checkEquals( res[[1]][[1]], matrix(as.integer((diag(3))),nr=3), msg = "eye<imat>(3,3)" )
    checkEquals( res[[1]][[2]], diag(3), msg = "eye<mat>(3,3)" )
    checkEquals( res[[1]][[3]], diag(3), msg = "eye<fmat>(3,3)" )
    checkEquals( res[[1]][[4]], matrix(as.integer((diag(3))),nr=3), msg = "eye<umat>(3,3)" )

    checkEquals( res[[2]][[1]], matrix(0, ncol = 5, nrow=1), msg = "zeros<mat>(5,1)" )
    checkEquals( res[[2]][[2]], matrix(0, ncol = 5, nrow=1), msg = "zeros<fmat>(5,1)" )

    checkEquals( res[[3]][[1]], matrix(0, ncol = 1, nrow=5), msg = "zeros<mat>(1,5)" )
    checkEquals( res[[3]][[2]], matrix(0, ncol = 1, nrow=5), msg = "zeros<mat>(1,5)" )

    checkEquals( res[[4]][[1]], matrix(0:3, ncol = 2, nrow=2), msg = "field<int>" )
    checkEquals( res[[4]][[2]], matrix(letters[1:4], ncol = 2, nrow=2), msg = "field<std::string>" )
}

test.wrap.Glue <- function(){
    
	fx <- cxxfunction( , '
	
    arma::mat m1 = arma::eye<arma::mat>( 3, 3 ) ;
    arma::mat m2 = arma::eye<arma::mat>( 3, 3 ) ;
	                     
    List res ;                                  
    res["mat+mat"] = m1 + m2 ;
    return res ;
	
	', plugin = "RcppArmadillo" )
	
	res <- fx()
    checkEquals( res[[1]], 2*diag(3), msg = "wrap(Glue)" )
}

test.wrap.Op <- function(){
	
	fx <- cxxfunction( , '
	
	arma::mat m1 = arma::eye<arma::mat>( 3, 3 ) ;
	
    List res ;
    res["- mat"] = - m1 ;
    return res ;
	
	', plugin = "RcppArmadillo" )
    res <- fx()
    checkEquals( res[[1]], -1*diag(3), msg = "wrap(Op)" )
}

test.as.Mat <- function(){
	
	fx <- cxxfunction( signature(input_ = "list" ) , '
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
    ', plugin = "RcppArmadillo" )
	
    integer_mat <- matrix( as.integer(diag(4)), ncol = 4, nrow = 4 )
    numeric_mat <- diag(5)
    res <- fx( list( integer_mat, numeric_mat ) )
    checkEquals( unlist( res), c(4L, 5L, 4L, 5L ), msg = "as<Mat>" )
}

test.as.Col <- function(){
	fx <- cxxfunction( signature(input_ = "list" ) , '
		
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
	
	', plugin = "RcppArmadillo" )
	
    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Col>" )
}

test.as.Row <- function(){
	fx <- cxxfunction( signature(input_ = "list" ) , '

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
	
	', plugin = "RcppArmadillo" )

    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Row>" )
}

test.cxmat <- function(){
	
	fx <- cxxfunction( signature() , '

    arma::cx_mat m1  = arma::eye<arma::cx_mat> ( 3, 3 ) ;
    arma::cx_fmat m2 = arma::eye<arma::cx_fmat>( 3, 3 ) ;
    return List::create( _["double"] = m1, _["float"] = m2 ) ;
    
	', plugin = "RcppArmadillo" )
	checkEquals( fx(), 
		list( double = (1+0i)*diag(3), float = (1+0i)*diag(3) ), 
		msg = "support for complex matrices" )
	
}

test.mtOp <- function(){

	fx <- cxxfunction( signature() , '

	std::complex<double> x( 1.0, 2.0 ) ;
    arma::mat m1  = arma::eye<arma::mat> ( 3, 3 ) ;

    return wrap( x * m1 ) ;
    
	', plugin = "RcppArmadillo" )
	checkEquals( fx(), 
		(1+2i)*diag(3), 
		msg = "support for mtOp" )
	
}

test.mtGlue <- function(){

	fx <- cxxfunction( signature() , '

	arma::imat m2 = arma::eye<arma::imat> ( 3, 3 ) ;
	arma::mat m1  = arma::eye<arma::mat> ( 3, 3 ) ;

    return wrap( m1 + m2 ) ;
    
	', plugin = "RcppArmadillo" )
	checkEquals( fx(), 
		2.0 * diag(3) , 
		msg = "support for mtOp" )
	
}


test.sugar <- function(){

	fx <- cxxfunction( signature(x= "numeric") , '
	NumericVector xx(x) ;
	arma::mat m = forward( xx + xx ) ; 
    return wrap( m ) ;
    
	', plugin = "RcppArmadillo" )
	checkEquals( fx(1:10), 
		matrix( 2*(1:10), nrow = 10 ) , 
		msg = "RcppArmadillo and sugar" )
	
}

test.sugar.cplx <- function(){

	fx <- cxxfunction( signature(x= "complex") , '
	ComplexVector xx(x) ;
	arma::cx_mat m = forward( exp( xx ) ) ; 
	
    return wrap( m ) ;
    
	', plugin = "RcppArmadillo" )
	x <- 1:10*(1+1i) 
	checkEquals( fx(x), 
		matrix( exp(x), nrow = 10 ) , 
		msg = "RcppArmadillo and sugar (complex)" )
	
}


