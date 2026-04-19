// rmultinom.h: Rcpp/Armadillo equivalent to R's stats::rmultinom().
//
// This is intended for use in C++ functions, and should *not* be called from R.
// It should yield identical results to R.
//
// Copyright (C)  2014-2026  Christian Gunning
// Copyright (C)  2026       Dirk Eddelbuettel and R Core
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

#ifndef RCPPARMADILLO__EXTENSIONS__MULTINOM_H
#define RCPPARMADILLO__EXTENSIONS__MULTINOM_H

#include <RcppArmadillo.h>
namespace Rcpp{
    namespace RcppArmadillo{

        IntegerVector rmultinom(int size, NumericVector prob) {
            // meaning of n, size, prob as in ?rmultinom
            // opposite of sample() - n=number of draws
            int probsize = prob.size();
            IntegerVector draws(probsize);            // Return object

            if (size < 0 || size == NA_INTEGER)
                Rcpp::stop( "Invalid size");

            // Using Kahan compensated summation for platform-independent accuracy
            // avoids relying on LDOUBLE (long double) which varies across platforms.
            // Code borrowed with full credits from R.
            double p_tot = 0.0,
                p_comp = 0.0; // Kahan compensation term

            for (int ii = 0; ii < probsize; ii++) {
                double pp = prob[ii];
                if (!R_FINITE(pp) || pp < 0. || pp > 1.) {
                    Rcpp::warning("Domain issue in rmultinom");
                    draws[ii] = NA_INTEGER;
                    return draws;
                }
                // Kahan summation: p_tot += pp with compensation
                double y = pp - p_comp,
                    t = p_tot + y;
                p_comp = (t - p_tot) - y;
                p_tot = t;
                draws[ii] = 0;
            }

            if (fabs((double)(p_tot - 1.0)) > 1e-7)
                Rcpp::stop("Probabilities do not sum to 1, please use FixProb");
            if (probsize == 1 && p_tot == 0.0) 	// trivial border case: do as rbinom
                return draws;
            if (size == 0 )             		// do as rbinom
                return draws;

            // rmultinom(size, REAL(prob), k, &INTEGER(ans)[ik]);
            // generate first draws-1 obs via binomials
            // for each slot
            for (int ii = 0; ii < probsize-1; ii++) { /* (p_tot, n) are for "remaining binomial" */
                if (prob[ii] != 0.) {
                    double pp = prob[ii] / p_tot;
                    // >= 1; > 1 happens because of rounding
                    draws[ii] = ((pp < 1.) ? (int) Rf_rbinom((double) size,  pp) : size);
                    size -= draws[ii];
                } else {
                    draws[ii] = 0;
                }
                // all done
                if (size <= 0) return draws;
                // Kahan subtraction: p_tot -= prob[k] with compensation
                double y = -prob[ii] - p_comp,
                    t = p_tot + y;
                p_comp = (t - p_tot) - y;
                p_tot = t;
            }
            // the rest go here
            draws[probsize-1] = size;
            return draws;
        }
    }
}

#endif
