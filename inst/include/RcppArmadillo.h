// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// RcppArmadillo.h: Rcpp/Armadillo glue
//
// Copyright (C)  2010 - 2013  Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

#if defined(Rcpp_hpp) && !defined(COMPILING_RCPPARMADILLO)
    #error "The file 'Rcpp.h' should not be included. Please correct to include only 'RcppArmadillo.h'."
#endif

#include <RcppArmadilloForward.h>
#include <Rcpp.h>
#include <RcppArmadilloWrap.h>
#include <RcppArmadilloAs.h>
#include <RcppArmadilloSugar.h>

/* ZGESDD - compute the singular value decomposition (SVD); of a   */
/* complex M-by-N matrix A, optionally computing the left and/or	   */
/* right singular vectors.  If singular vectors are desired, it uses a */
/* divide-and-conquer algorithm.				   */
// namespace {
//     extern "C" void F77_NAME(zgesdd)(const char *jobz, const int *m, const int *n,
//                                      Rcomplex *a, const int *lda, double *s,
//                                      Rcomplex *u, const int *ldu,
//                                      Rcomplex *vt, const int *ldvt,
//                                      Rcomplex *work, const int *lwork, 
//                                      double *rwork,
//                                      int *iwork, int *info);
// }

#endif

