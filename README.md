
## RcppArmadillo: R and Armadillo via Rcpp

[![Build Status](https://travis-ci.org/RcppCore/RcppArmadillo.svg)](https://travis-ci.org/RcppCore/RcppArmadillo) [![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html) [![CRAN](https://www.r-pkg.org/badges/version/RcppArmadillo)](https://cran.r-project.org/package=RcppArmadillo) [![Dependencies](https://tinyverse.netlify.com/badge/RcppArmadillo)](https://cran.r-project.org/package=RcppArmadillo) [![Downloads](https://cranlogs.r-pkg.org/badges/RcppArmadillo?color=brightgreen)](https://www.r-pkg.org/pkg/RcppArmadillo) [![StackOverflow](https://img.shields.io/badge/stackoverflow-rcpp-orange.svg)](https://stackoverflow.com/questions/tagged/rcpp)

### Overview

[Armadillo](http://arma.sf.net) is a templated C++ linear algebra library
written by Conrad Sanderson that aims towards a good balance between speed and ease of use. Integer,
floating point and complex numbers are supported, as well as a subset of
trigonometric and statistics functions. Various matrix decompositions are
provided through optional integration with LAPACK and ATLAS libraries.

A delayed evaluation approach is employed (during compile time) to combine
several operations into one, and to reduce (or eliminate) the need for
temporaries. This is accomplished through recursive templates and template
meta-programming.

This library is useful if C++ has been decided as the language of choice
(due to speed and/or integration capabilities), rather than another language.

The RcppArmadillo package includes the header files from the templated
Armadillo library. Thus users do not need to install Armadillo itself in
order to use RcppArmadillo.

This Armadillo integration provides a nice illustration of the
capabilities of the [Rcpp](http://www.rcpp.org) package for seamless R and
C++ integration.

### Status

The package is under active development with releases to
[CRAN](https://cran.r-project.org) about once every other month, and
widely-used by other CRAN packages as can be seen from the
[CRAN package page](https://cran.r-project.org/package=RcppArmadillo). 
As of September 2019, there are 658 CRAN packages using RcppArmadillo.

### Documentation

The package contains a pdf vignette which is a pre-print of the
[paper by Eddelbuettel and Sanderson](http://dx.doi.org/10.1016/j.csda.2013.02.005)
in CSDA (2014), as well as an introductory vignette for the sparse
matrix conversions.

### Authors

Dirk Eddelbuettel, Romain Francois, Doug Bates and Binxiang Ni

### License

GPL (>= 2)
