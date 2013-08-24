// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// sample.cpp: Rcpp/Armadillo equivalent to R's sample().  
// This is to be used in C++ functions, and should *not* be called from R.
// It should yield identical results to R in most cases
// (note that Walker's alias method is not implemented).
//
// Copyright (C)  2012 - 2013  Christian Gunning
// Copyright (C)  2013  Romain Francois
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

#ifndef RCPPARMADILLO__EXTENSIONS__SAMPLE_H
#define RCPPARMADILLO__EXTENSIONS__SAMPLE_H

#include <RcppArmadillo.h>

namespace Rcpp{

    namespace RcppArmadillo{

        //Declarations
        template <typename INDEX>
        void SampleReplace( INDEX &index, int nOrig, int size);
        
        template <typename INDEX>
        void SampleNoReplace( INDEX &index, int nOrig, int size);
        
        template <typename INDEX>
        void ProbSampleReplace(INDEX &index, int nOrig, int size, arma::vec &prob);
        
        template <typename INDEX>
        void ProbSampleNoReplace(INDEX &index, int nOrig, int size, arma::vec &prob);
        
        template <typename PROB>
        void FixProb(PROB &prob, int size, bool replace);

        template <class T> 
        T sample(const T &x, const int size, const bool replace, NumericVector prob_ = NumericVector(0) ) {
            // Templated sample -- should work on any Rcpp Vector
            int ii, jj;
            int nOrig = x.size();
            int probsize = prob_.size();
            // Create return object
            T ret(size);
            if ( size > nOrig && !replace) throw std::range_error( "Tried to sample more elements than in x without replacement" ) ;
            // Store the sample ids here, modify in-place
            IntegerVector index(size);
            if (probsize == 0) { // No probabilities given
                if (replace) {
                    SampleReplace<IntegerVector>(index, nOrig, size);
                } else {
                    SampleNoReplace<IntegerVector>(index, nOrig, size);
                }
            } else { 
                if (probsize != nOrig) throw std::range_error( "Number of probabilities must equal input vector length" ) ;
                // normalize, error-check probability vector
                // fixprob will be modified in-place
                NumericVector fixprob = clone(prob_);
                FixProb<NumericVector>(fixprob, size, replace);
                // 
                // Copy the given probabilities into an arma vector
                arma::vec prob(fixprob.begin(), fixprob.size());
                if (replace) {
                    // check for walker alias conditions 
                    int walker_test = sum( (probsize * prob) > 0.1);
                    if (walker_test < 200) { 
                        ProbSampleReplace<IntegerVector>(index, nOrig, size, prob);
                    } else {
                        throw std::range_error( "Walker Alias method not implemented. R-core sample() is likely faster for this problem.");
                        // WalkerProbSampleReplace(index, nOrig, size, prob);
                    }
                } else {
                    ProbSampleNoReplace<IntegerVector>(index, nOrig, size, prob);
                }
            }
            // copy the results into the return vector
            for (ii=0; ii<size; ii++) {
                jj = index[ii];
                ret[ii] = x[jj];
            }
            return(ret);
        }

        template <typename INDEX>
        void SampleReplace( INDEX &index, int nOrig, int size) {
            int ii;
            for (ii = 0; ii < size; ii++) {
                index[ii] = nOrig * unif_rand();
            }
        }

        template <typename INDEX>
        void SampleNoReplace( INDEX &index, int nOrig, int size) {
            int ii, jj;
            IntegerVector sub(nOrig);
            for (ii = 0; ii < nOrig; ii++) {
                sub[ii] = ii;
            }
            for (ii = 0; ii < size; ii++) {
                jj = nOrig * unif_rand();
                index[ii] = sub[jj];
                // replace sampled element with last, decrement
                sub[jj] = sub[--nOrig];
            }
        }

        template <typename PROB>
        void FixProb(PROB &prob, const int size, const bool replace) {
            // prob is modified in-place.  
            // Is this better than cloning it?
            double sum = 0.0;
            int ii, nPos = 0;
            int nn = prob.size();
            for (ii = 0; ii < nn; ii++) {
                if (!R_FINITE(prob[ii])) //does this work??
                    throw std::range_error( "NAs not allowed in probability" ) ;
                if (prob[ii] < 0.0)
                    throw std::range_error( "Negative probabilities not allowed" ) ;
                if (prob[ii] > 0.0) {
                    nPos++;
                    sum += prob[ii];
                }
            }
            if (nPos == 0 || (!replace && size > nPos)) {
                throw std::range_error("Not enough positive probabilities");
            }
            prob = prob / sum;  //sugar
        }

        // Unequal probability sampling with replacement 
        template <typename INDEX>
        void ProbSampleReplace(INDEX &index, int nOrig, int size, arma::vec &prob){
            double rU;
            int ii, jj;
            int nOrig_1 = nOrig - 1;
            arma::uvec perm = arma::sort_index(prob, 1); //descending sort of index
            prob = arma::sort(prob, 1);  // descending sort of prob
            // cumulative probabilities 
            prob = arma::cumsum(prob);
            // compute the sample 
            for (ii = 0; ii < size; ii++) {
                rU = unif_rand();
                for (jj = 0; jj < nOrig_1; jj++) {
                    if (rU <= prob[jj])
                        break;
                }
                index[ii] = perm[jj];
            }
        }

        // Unequal probability sampling without replacement 
        template <typename INDEX>
        void ProbSampleNoReplace(INDEX &index, int nOrig, int size, arma::vec &prob){
            int ii, jj, kk;
            int nOrig_1 = nOrig - 1;
            double rT, mass, totalmass = 1.0;
            arma::uvec perm = arma::sort_index(prob, 1); //descending sort of index
            prob = arma::sort(prob, 1);  // descending sort of prob
            // compute the sample 
            for (ii = 0; ii < size; ii++, nOrig_1--) {
                rT = totalmass * unif_rand();
                mass = 0;
                for (jj = 0; jj < nOrig_1; jj++) {
                    mass += prob[jj];
                    if (rT <= mass)
                        break;
                }
                index[ii] = perm[jj];
                totalmass -= prob[jj];
                for ( kk = jj; kk < nOrig_1; kk++) {
                    prob[kk] = prob[kk+1];
                    perm[kk] = perm[kk+1];
                }
            }
        }
    }
}

#endif
