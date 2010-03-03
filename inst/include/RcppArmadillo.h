// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppArmadillo.h: Rcpp/Armadillo glue
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

#ifndef RcppArmadillo__RcppArmadillo__h
#define RcppArmadillo__RcppArmadillo__h

#include <RcppCommon.h>
#include <armadillo>
#include <RcppArmadilloDefines.h>

/* forward declarations */
namespace Rcpp {
	
    /* support for wrap */
    template <typename T> SEXP wrap ( const arma::Mat<T>& ) ;
    template <typename T> SEXP wrap ( const arma::Row<T>& ) ;
    template <typename T> SEXP wrap ( const arma::Col<T>& ) ;
    template <typename T> SEXP wrap ( const arma::field<T>& ) ;
    #if ARMA_VERSION_GE_070
    template <typename T> SEXP wrap ( const arma::Cube<T>& ) ;
    #endif

    template <typename T1, typename T2, typename glue_type> 
    SEXP wrap(const arma::Glue<T1, T2, glue_type>& X ) ;
    
    template <typename T1, typename op_type>
    SEXP wrap(const arma::Op<T1, op_type>& X ) ;
    
    #if ARMA_VERSION_GE_090
    template <typename T1, typename T2, typename glue_type> 
    SEXP wrap(const arma::eGlue<T1, T2, glue_type>& X ) ;
    
    template <typename T1, typename op_type>
    SEXP wrap(const arma::eOp<T1, op_type>& X ) ;
    #endif 
    
    namespace traits {

	/* support for as */
	template <typename T> class Exporter< arma::Mat<T> > ;
	template <typename T> class Exporter< arma::Row<T> > ;
	template <typename T> class Exporter< arma::Col<T> > ;
// template <typename T> class Exporter< arma::field<T> > ;
// #ifdef ARMA_VERSION_GE_070
// 	template <typename T> class Exporter< arma::Cube<T> > ;
// #endif

    } // namespace traits 

}

#include <Rcpp.h>

RcppExport SEXP RcppArmadilloExample() ;
RcppExport SEXP RcppArmadilloExample_as_Mat( SEXP );
RcppExport SEXP RcppArmadilloExample_as_Col( SEXP );
RcppExport SEXP RcppArmadilloExample_as_Row( SEXP );

namespace Rcpp{

    namespace RcppArmadillo{

	template <typename T>
	SEXP arma_wrap( const T& object, const ::Rcpp::Dimension& dim){
	    ::Rcpp::RObject x = ::Rcpp::wrap( object.memptr() , object.memptr() + object.n_elem ) ;
	    x.attr( "dim" ) = dim ;
	    return x; 
	}

    } /* namespace RcppArmadillo */
	
    /* wrap */

    template <typename T> SEXP wrap ( const arma::Mat<T>& data ){
	return RcppArmadillo::arma_wrap( data, Dimension( data.n_rows, data.n_cols ) ) ;
    } ;

    template <typename T> SEXP wrap( const arma::Col<T>& data ){
	return RcppArmadillo::arma_wrap( data, Dimension( data.n_elem, 1) ) ;
    } ;

    template <typename T> SEXP wrap( const arma::Row<T>& data ){
	return RcppArmadillo::arma_wrap(data, Dimension( 1, data.n_elem ) ) ;
    } ;

    #if ARMA_VERSION_GE_070
    template <typename T> SEXP wrap( const arma::Cube<T>& data ){
	return RcppArmadillo::arma_wrap(data, Dimension(  data.n_rows, data.n_cols, data.n_slices ) ) ;
    }
    #endif

    namespace RcppArmadillo {
	
	/* Importer class for field<T> */
	template <typename T> class FieldImporter {
	public:
	    typedef T r_import_type ;
	    FieldImporter( const arma::field<T>& data_ ) : data(data_){}
	    inline int size() const { return data.n_elem ; }
	    inline T get(int i) const { return data[i] ; }
	    inline SEXP wrap( int i) const { return ::Rcpp::wrap( data[i] ) ; } 
	private:
	    const arma::field<T>& data ;
	} ;

    } // namespace RcppArmadillo

    template <typename T> 
    SEXP wrap( const arma::field<T>& data){
	RObject x = wrap( RcppArmadillo::FieldImporter<T>( data ) ) ;
	x.attr("dim" ) = Dimension( data.n_rows, data.n_cols ) ;
	return x ;
    }
    
    /* TODO: maybe we could use the advanced constructor to avoid creating the 
             temporary Mat */
    template <typename T1, typename T2, typename glue_type>
    SEXP wrap(const arma::Glue<T1, T2, glue_type>& X ){
    	    return wrap( arma::Mat<typename T1::elem_type>(X) ) ;
    }

    template <typename T1, typename op_type>
    SEXP wrap(const arma::Op<T1, op_type>& X ){
    	    return wrap( arma::Mat<typename T1::elem_type>(X) ) ;
    }
    
    /* TODO: will do better when I can use 0.9.0 */
    #if ARMA_VERSION_GE_090
    namespace RcppArmadillo{
    	
    	/* we can intercept and directly build the resulting matrix using 
    	   memory allocated by R */
    	template <typename T1, typename T2, typename eglue_type>
    	SEXP wrap_eglue( const arma::eGlue<T1, T2, eglue_type>& X, ::Rcpp::traits::false_type ){
    		int n_rows = X.P1.n_rows ;
    		int n_cols = X.P1.n_cols ;
    		typedef typename ::Rcpp::Vector< ::Rcpp::traits::r_sexptype_traits< typename T1::elem_type>::rtype > VECTOR ;
    		VECTOR res(::Rcpp::Dimension( n_rows , n_cols )) ;
    		::arma::Mat<typename T1::elem_type> result( res.begin(), n_rows, n_cols, false ) ;
    		result = X ;
    		return res ;
    	}
    	
    	template <typename T1, typename T2, typename eglue_type>
    	SEXP wrap_eglue( const arma::eGlue<T1, T2, eglue_type>& X, ::Rcpp::traits::true_type ){
    		return ::Rcpp::wrap( arma::Mat<typename T1::elem_type>(X) ) ;
    	}
    	
    	template <typename T1, typename eop_type>
    	SEXP wrap_eop( const arma::eOp<T1,eop_type>& X, ::Rcpp::traits::false_type ){
    		int n_rows = X.P.n_rows ;
    		int n_cols = X.P.n_cols ;
    		typedef typename ::Rcpp::Vector< ::Rcpp::traits::r_sexptype_traits< typename T1::elem_type>::rtype > VECTOR ;
    		VECTOR res(::Rcpp::Dimension( n_rows , n_cols )) ;
    		::arma::Mat<typename T1::elem_type> result( res.begin(), n_rows, n_cols, false ) ;
    		result = X ;
    		return res ;
    	}
    	
    	template <typename T1, typename eop_type>
    	SEXP wrap_eop( const arma::eOp<T1,eop_type>& X, ::Rcpp::traits::true_type ){
    		return ::Rcpp::wrap( arma::Mat<typename T1::elem_type>(X) ) ;
    	}
    	
    } // namespace RcppArmadillo
    
    template <typename T1, typename T2, typename glue_type>
    SEXP wrap(const arma::eGlue<T1, T2, glue_type>& X ){
    	return RcppArmadillo::wrap_eglue( X, typename traits::r_sexptype_needscast<typename T1::elem_type>::type() ) ;
    }

    template <typename T1, typename op_type>
    SEXP wrap(const arma::eOp<T1, op_type>& X ){
    	    return RcppArmadillo::wrap_eop( X, typename traits::r_sexptype_needscast<typename T1::elem_type>::type() ) ;
    }
    #endif

    /* support for Rcpp::as */

    namespace traits {
	
	template <typename T> 
	class Exporter< arma::Col<T> > : public IndexingExporter< arma::Col<T>, T > {
	public: 
	    Exporter(SEXP x) : IndexingExporter< arma::Col<T>, T >(x){}
	}; 

	template <typename T> 
	class Exporter< arma::Row<T> > : public IndexingExporter< arma::Row<T>, T > {
	public:
	    Exporter(SEXP x) : IndexingExporter< arma::Row<T>, T >(x){}
	}; 

	template <typename T> 
	class Exporter< arma::Mat<T> > : public MatrixExporter< arma::Mat<T>, T > {
	public:
	    Exporter(SEXP x) : MatrixExporter< arma::Mat<T>, T >(x){}
	}; 
	
    } // namespace traits

}

#endif

