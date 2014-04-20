// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// RcppArmadilloAs.h: Rcpp/Armadillo glue, support for as
//
// Copyright (C)  2013  Dirk Eddelbuettel and Romain Francois
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
    class Exporter< const arma::Mat<T>& > {
    public:  
        typedef typename Rcpp::Matrix< Rcpp::traits::r_sexptype_traits<T>::rtype > MATRIX ;
                
        Exporter(SEXP x) : mat(x) {}
                
        inline arma::Mat<T>* get(){
            return new arma::Mat<T>( mat.begin(), mat.nrow(), mat.ncol(), false ) ;
        }
                
    private:
        MATRIX mat ;
    };
         
    template <typename T>
    class Exporter< arma::SpMat<T> > {
    public:
        Exporter( SEXP x ) : mat(x){}
                
        arma::SpMat<T> get(){
            const int  RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
        
            IntegerVector dims = mat.slot("Dim");
            IntegerVector i = mat.slot("i") ;
            IntegerVector p = mat.slot("p") ;     
            Vector<RTYPE> x = mat.slot("x") ;
                        
            arma::SpMat<T> res(dims[0], dims[1]);
                        
            // create space for values, and copy
            arma::access::rw(res.values) = arma::memory::acquire_chunked<T>(x.size() + 1);
            arma::arrayops::copy(arma::access::rwp(res.values), x.begin(), x.size() + 1);
                        
            // create space for row_indices, and copy 
            arma::access::rw(res.row_indices) = arma::memory::acquire_chunked<arma::uword>(i.size() + 1);
            std::copy( i.begin(), i.end(), arma::access::rwp(res.row_indices) ) ;
                        
            // create space for col_ptrs, and copy 
            arma::access::rw(res.col_ptrs) = arma::memory::acquire<arma::uword>(p.size() + 2);
            std::copy(p.begin(), p.end(), arma::access::rwp(res.col_ptrs) );
                        
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
        
    /* Begin Armadillo vector as support classes */
    
    template <typename T, typename MAT, typename REF, 
      typename NEEDS_CAST = typename Rcpp::traits::r_sexptype_needscast<T>::type 
    >
    class ArmaMat_InputParameter;
    
    template <typename T, typename MAT, typename REF>
    class ArmaMat_InputParameter<T, MAT, REF, Rcpp::traits::false_type> {
    public:
        ArmaMat_InputParameter( SEXP x_ ) : m(x_), mat( reinterpret_cast<T*>( m.begin() ), m.nrow(), m.ncol(), false) {} 
                        
        inline operator REF(){
            return mat ;        
        }
                        
    private:
        Rcpp::Matrix< Rcpp::traits::r_sexptype_traits<T>::rtype > m ;
        MAT mat ;
    } ;
    
    template <typename T, typename MAT, typename REF>
    class ArmaMat_InputParameter<T, MAT, REF, Rcpp::traits::true_type> {
    public:
        ArmaMat_InputParameter( SEXP x_ ): m(x_), mat( as<MAT>(m) ) {}
                        
        inline operator REF(){
            return mat ;        
        }
                        
    private:
        Rcpp::Matrix< Rcpp::traits::r_sexptype_traits<T>::rtype > m ;
        MAT mat ;
    } ;

    /* End Armadillo vector as support classes */


    /* Begin Armadillo vector as support classes */
    
    template <typename T, typename VEC, typename REF, 
      typename NEEDS_CAST = typename Rcpp::traits::r_sexptype_needscast<T>::type 
    >
    class ArmaVec_InputParameter;
    
    template <typename T, typename VEC, typename REF>
    class ArmaVec_InputParameter<T, VEC, REF, Rcpp::traits::false_type> {
    public:
        ArmaVec_InputParameter( SEXP x_ ) : v(x_), vec( reinterpret_cast<T*>( v.begin() ), v.size(), false ){}
                        
        inline operator REF(){
            return vec ;        
        }
                        
    private:
        Rcpp::Vector< Rcpp::traits::r_sexptype_traits<T>::rtype > v ;
        VEC vec ;
    } ;
    
    template <typename T, typename VEC, typename REF>
    class ArmaVec_InputParameter<T, VEC, REF, Rcpp::traits::true_type> {
    public:
        ArmaVec_InputParameter( SEXP x_ ): v(x_), vec( as<VEC>(v) ) {}
                        
        inline operator REF(){
            return vec ;        
        }
                        
    private:
        Rcpp::Vector< Rcpp::traits::r_sexptype_traits<T>::rtype > v ;
        VEC vec ;
    } ;

    /* End Armadillo vector as support classes */
    
#define MAKE_INPUT_PARAMETER(INPUT_TYPE,TYPE,REF)                       \
    template <typename T>                                               \
    class INPUT_TYPE<TYPE> : public ArmaVec_InputParameter<T, TYPE, REF >{ \
    public:                                                             \
    INPUT_TYPE( SEXP x) : ArmaVec_InputParameter<T, TYPE, REF >(x){} \
    } ;                                                                                                                  
    
    MAKE_INPUT_PARAMETER(ConstReferenceInputParameter, arma::Col<T>, const arma::Col<T>& )
    MAKE_INPUT_PARAMETER(ReferenceInputParameter     , arma::Col<T>, arma::Col<T>&       )
    MAKE_INPUT_PARAMETER(ConstInputParameter         , arma::Col<T>, const arma::Col<T>  )
    
    MAKE_INPUT_PARAMETER(ConstReferenceInputParameter, arma::Row<T>, const arma::Row<T>& )
    MAKE_INPUT_PARAMETER(ReferenceInputParameter     , arma::Row<T>, arma::Row<T>&       )
    MAKE_INPUT_PARAMETER(ConstInputParameter         , arma::Row<T>, const arma::Row<T>  )
    
#undef MAKE_INPUT_PARAMETER

    
#define MAKE_INPUT_PARAMETER(INPUT_TYPE,TYPE,REF)                       \
    template <typename T>                                               \
    class INPUT_TYPE<TYPE> : public ArmaMat_InputParameter<T, TYPE, REF >{ \
    public:                                                             \
    INPUT_TYPE( SEXP x) : ArmaMat_InputParameter<T, TYPE, REF >(x){} \
    } ;
 
    MAKE_INPUT_PARAMETER(ConstReferenceInputParameter, arma::Mat<T>, const arma::Mat<T>& )
    MAKE_INPUT_PARAMETER(ReferenceInputParameter     , arma::Mat<T>, arma::Mat<T>&       )
    MAKE_INPUT_PARAMETER(ConstInputParameter         , arma::Mat<T>, const arma::Mat<T>  )

#undef MAKE_INPUT_PARAMETER
    
}

#endif

