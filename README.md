
## RcppArmadillo: R and Armadillo via Rcpp

[![CI](https://github.com/RcppCore/RcppArmadillo/workflows/ci/badge.svg)](https://github.com/RcppCore/RcppArmadillo/actions?query=workflow%3Aci)
[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](https://www.r-pkg.org/badges/version/RcppArmadillo)](https://cran.r-project.org/package=RcppArmadillo)
[![Dependencies](https://tinyverse.netlify.com/badge/RcppArmadillo)](https://cran.r-project.org/package=RcppArmadillo)
[![Coverage Status](https://codecov.io/gh/RcppCore/RcppArmadillo/graph/badge.svg)](https://app.codecov.io/github/RcppCore/RcppArmadillo?branch=master)
[![Debian package](https://img.shields.io/debian/v/r-cran-rcpparmadillo/sid?color=brightgreen)](https://packages.debian.org/sid/r-cran-rcpparmadillo)
[![Last Commit](https://img.shields.io/github/last-commit/RcppCore/RcppArmadillo)](https://github.com/RcppCore/RcppArmadillo)
[![Downloads (monthly)](https://cranlogs.r-pkg.org/badges/RcppArmadillo?color=brightgreen)](https://www.r-pkg.org/pkg/RcppArmadillo)
[![Downloads (total)](https://cranlogs.r-pkg.org/badges/grand-total/RcppArmadillo?color=brightgreen)](https://www.r-pkg.org/pkg/RcppArmadillo)
[![CRAN use](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/CRAN_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![CRAN indirect](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/indirect_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![BioConductor use](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/BioC_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![StackOverflow](https://img.shields.io/badge/stackoverflow-rcpp-orange.svg)](https://stackoverflow.com/questions/tagged/rcpp)
[![CSDA](https://img.shields.io/badge/CSDA-10.1016%2Fj.csda.2013.02.005-brightgreen)](https://doi.org/10.1016/j.csda.2013.02.005)

### Synopsis

RcppArmadillo provides an interface from R to and from [Armadillo](https://arma.sourceforge.net) by
utilising the [Rcpp R/C++ interface library](http://dirk.eddelbuettel.com/code/rcpp.html).

### What is Armadillo?

[Armadillo](https://arma.sourceforge.net) is a high-quality linear algebra library for the C++ language,
aiming towards a good balance between speed and ease of use. It provides high-level syntax and
[functionality](https://arma.sourceforge.net/docs.html) deliberately similar to Matlab (TM).
See [its website](https://arma.sourceforge.net) more information about Armadillo.

### So give me an example!

Glad you asked. Here is a light-weight and fast implementation of linear regression:

```c++
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List fastLm(const arma::mat& X, const arma::colvec& y) {
    int n = X.n_rows, k = X.n_cols;

    arma::colvec coef = arma::solve(X, y);     // fit model y ~ X
    arma::colvec res  = y - X*coef;            // residuals
    double s2 = arma::dot(res, res) / (n - k); // std.errors of coefficients
    arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));

    return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                              Rcpp::Named("stderr")       = std_err,
                              Rcpp::Named("df.residual")  = n - k);
}
```

You can
[`Rcpp::sourceCpp()`](https://cran.r-project.org/package=Rcpp/vignettes/Rcpp-attributes.pdf)
the file above to compile the function.  A version is also included in the
package [as the `fastLm()`](https://github.com/RcppCore/RcppArmadillo/blob/master/R/fastLm.R)
function.

### Status

The package is mature yet under active development with releases to
[CRAN](https://cran.r-project.org) about once every other month, and
widely-used by other CRAN packages as can be seen from the
[CRAN package page](https://cran.r-project.org/package=RcppArmadillo).
As of January 2023, there are 1035 CRAN packages using RcppArmadillo.

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

Dirk Eddelbuettel, Romain Francois, Doug Bates, Binxiang Ni, and Conrad Sanderson

### License

GPL (>= 2)
