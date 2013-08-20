// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// RcppArmadilloAs.h: Rcpp/Armadillo glue, support for as
//
// Copyright (C)  2013  Dirk Eddelbuettel, Romain Francois
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

#ifndef RcppArmadillo__RcppArmadilloAs__h
#define RcppArmadillo__RcppArmadilloAs__h

namespace Rcpp{
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
		
	 
	template <typename T>
	class Exporter< arma::SpMat<T> > {
	public:
		Exporter( SEXP x ) : mat(x){}
		
		arma::SpMat<T> get(){
			IntegerVector dims = mat.slot("Dim");
			arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
			arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));     
			arma::Col<T> x     = Rcpp::as< arma::Col<T> >(mat.slot("x"));
			
			int nrow = dims[0], ncol = dims[1];
			arma::SpMat<T> res(nrow, ncol);
			
			// create space for values, and copy
			arma::access::rw(res.values) = arma::memory::acquire_chunked<T>(x.size() + 1);
			arma::arrayops::copy(arma::access::rwp(res.values), x.begin(), x.size() + 1);
			
			// create space for row_indices, and copy 
			arma::access::rw(res.row_indices) = 
			    arma::memory::acquire_chunked<arma::uword>(i.size() + 1);
			arma::arrayops::copy(arma::access::rwp(res.row_indices), i.begin(), i.size() + 1);
			
			// create space for col_ptrs, and copy 
			arma::access::rw(res.col_ptrs) = arma::memory::acquire<arma::uword>(p.size() + 2);
			arma::arrayops::copy(arma::access::rwp(res.col_ptrs), p.begin(), p.size() + 1);
			
			// important: set the sentinel as well
			arma::access::rwp(res.col_ptrs)[p.size()+1] = std::numeric_limits<arma::uword>::max();
			
			// set the number of non-zero elements
			arma::access::rw(res.n_nonzero) = x.size();
			
			return res;
		}
		
	private:
		S4 mat ;
	} ;
}	
}

#endif

