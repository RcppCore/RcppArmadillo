
## RcppArmadillo: R and Armadillo via Rcpp

[![Build Status](https://travis-ci.org/RcppCore/RcppArmadillo.svg)](https://travis-ci.org/RcppCore/RcppArmadillo)
[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](https://www.r-pkg.org/badges/version/RcppArmadillo)](https://cran.r-project.org/package=RcppArmadillo)
[![Dependencies](https://tinyverse.netlify.com/badge/RcppArmadillo)](https://cran.r-project.org/package=RcppArmadillo)
[![Debian package](https://img.shields.io/debian/v/r-cran-rcpparmadillo/sid?color=brightgreen)](https://packages.debian.org/sid/r-cran-rcpparmadillo)
[![Last Commit](https://img.shields.io/github/last-commit/RcppCore/RcppArmadillo)](https://github.com/RcppCore/RcppArmadillo)  
[![Downloads](https://cranlogs.r-pkg.org/badges/RcppArmadillo?color=brightgreen)](https://www.r-pkg.org/pkg/RcppArmadillo)
[![CRAN use](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/CRAN_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![CRAN indirect](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/indirect_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![BioConductor use](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/BioC_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![StackOverflow](https://img.shields.io/badge/stackoverflow-rcpp-orange.svg)](https://stackoverflow.com/questions/tagged/rcpp)
[![CSDA](https://img.shields.io/badge/CSDA-10.1016%2Fj.csda.2013.02.005-brightgreen)](https://doi.org/10.1016/j.csda.2013.02.005)

### Synopsis

RcppArmadillo provides an interface from R to and from [Armadillo](http://arma.sourceforge.net) by
utilising the [Rcpp R/C++ interface library](http://dirk.eddelbuettel.com/code/rcpp.html).

### What is Armadillo?

[Armadillo](http://arma.sourceforge.net) is a high-quality linear algebra library for the C++ language,
aiming towards a good balance between speed and ease of use. It provides high-level syntax and
[functionality](http://arma.sourceforge.net/docs.html) deliberately similar to Matlab (TM).
See [its website](http://arma.sourceforge.net) more information about Armadillo.

### So give me an example!

Glad you asked. Here is a light-weight and fast implementation of linear regression:

```c++
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List fastLm(const arma::mat& X, const arma::colvec& y) {
    int n = X.n_rows, k = X.n_cols;

    arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
    arma::colvec res  = y - X*coef;           // residuals

    // std.errors of coefficients
    double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);

    arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));  

    return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                              Rcpp::Named("stderr")       = std_err,
                              Rcpp::Named("df.residual")  = n - k);
}
```

You can
[`Rcpp::sourceCpp()`](https://cran.r-project.org/package=Rcpp/vignettes/Rcpp-attributes.pdf)
the file above to compile the function.  A slightly more involved version is also included in the
package [as the `fastLm()`](https://github.com/RcppCore/RcppArmadillo/blob/master/R/fastLm.R)
function.

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

### Installation

RcppArmadillo is a [CRAN package](https://cran.r-project.org/package=RcppArmadillo), and lives
otherwise in its own habitat on [GitHub](https://github.com/RcppCore/RcppArmadillo) within the
[RcppCore](https://github.com/RcppCore) GitHub organization.

Run

```r
install.packages("RcppArmadillo")
```

to install from your nearest CRAN mirror.

### Authors

Dirk Eddelbuettel, Romain Francois, Doug Bates and Binxiang Ni

### License

GPL (>= 2)
