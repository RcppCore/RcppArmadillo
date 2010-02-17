
RcppArmadilloExample <- function(){
	.Call( "RcppArmadilloExample", PACKAGE = "RcppArmadillo" )
}

RcppArmadilloExample_as_Mat <- function(){
	integer_mat <- matrix( as.integer(diag(4)), ncol = 4, nrow = 4 )
	numeric_mat <- diag(5)
	.Call( "RcppArmadilloExample_as_Mat", 
		list( integer_mat, numeric_mat ), 
		PACKAGE = "RcppArmadillo" )
}

RcppArmadilloExample_as_Col <- function(){
	.Call( "RcppArmadilloExample_as_Col", 
		list( 1:10, as.numeric(1:10) ), 
		PACKAGE = "RcppArmadillo" )
}

RcppArmadilloExample_as_Row <- function(){
	.Call( "RcppArmadilloExample_as_Row", 
		list( 1:10, as.numeric(1:10) ), 
		PACKAGE = "RcppArmadillo" )
}

