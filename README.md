
## RcppArmadillo: R and Armadillo via Rcpp

[![CI](https://github.com/RcppCore/RcppArmadillo/workflows/ci/badge.svg)](https://github.com/RcppCore/RcppArmadillo/actions?query=workflow%3Aci)
[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](https://www.r-pkg.org/badges/version/RcppArmadillo)](https://cran.r-project.org/package=RcppArmadillo)
[![Debian package](https://img.shields.io/debian/v/r-cran-rcpparmadillo/sid?color=brightgreen)](https://packages.debian.org/sid/r-cran-rcpparmadillo)
[![r-universe](https://rcppcore.r-universe.dev/badges/RcppArmadillo)](https://rcppcore.r-universe.dev/RcppArmadillo)
[![Dependencies](https://tinyverse.netlify.app/badge/RcppArmadillo)](https://cran.r-project.org/package=RcppArmadillo)
[![Coverage Status](https://codecov.io/gh/RcppCore/RcppArmadillo/graph/badge.svg)](https://app.codecov.io/github/RcppCore/RcppArmadillo?branch=master)
[![Last Commit](https://img.shields.io/github/last-commit/RcppCore/RcppArmadillo)](https://github.com/RcppCore/RcppArmadillo)
[![Downloads (monthly)](https://cranlogs.r-pkg.org/badges/RcppArmadillo?color=brightgreen)](https://www.r-pkg.org/pkg/RcppArmadillo)
[![Downloads (total)](https://cranlogs.r-pkg.org/badges/grand-total/RcppArmadillo?color=brightgreen)](https://www.r-pkg.org/pkg/RcppArmadillo)
[![CRAN use](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/CRAN_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![CRAN indirect](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/indirect_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![BioConductor use](https://jangorecki.gitlab.io/rdeps/RcppArmadillo/BioC_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppArmadillo)
[![CSDA](https://img.shields.io/badge/CSDA-10.1016%2Fj.csda.2013.02.005-brightgreen)](https://doi.org/10.1016/j.csda.2013.02.005)

### Synopsis

RcppArmadillo provides an interface from R to and from [Armadillo][armadillo] by utilising the [Rcpp
R/C++ interface library][rcpp].

### What is Armadillo?

[Armadillo][armadillo] is a high-quality linear algebra library for the C++ language, aiming towards
a good balance between speed and ease of use. It provides high-level syntax and
[functionality](https://arma.sourceforge.net/docs.html) deliberately similar to Matlab (TM).  See
[its website][armadillo] more information about Armadillo.

### So give me an example!

Glad you asked. Here is a light-weight and fast implementation of linear regression:

```c++
#include <RcppArmadillo/Lighter>
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

The `RcppArmadillo/Lighter` header includes [Rcpp][rcpp] via its `Rcpp/Lighter` header which
precludes some more compile-time heavy features such as 'Rcpp Modules' which we may not need. See
the [Rcpp][rcpp] docs more details about 'Light', 'Lighter' and 'Lightest'.  In the example above,
the switch saves about 15% of total compilation time.

### Status

The package is mature yet under active development with releases to [CRAN][cran] about once every
other month, and widely-used by other CRAN packages as can be seen from the [CRAN package page][cran
pkg].  As of August 2025, there are 1266 CRAN packages using RcppArmadillo.

As of Armadillo 15.0.0, the minimum compilation standard is C++14. However, as several hundred CRAN
packages still impose C++11 as their compilation standard, the RcppArmadillo package also includes
the final version allowing C++11, namely Armadillo 14.6.3, as a fallback used when C++11 compilation
is detected. Conversion to and compilation under C++14 or later is encouraged. R [defaults to
C++17](https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Portable-C-and-C_002b_002b-code-1)
since version 4.3.0. See [GitHub issue #475](https://github.com/RcppCore/RcppArmadillo/issues/475)
for more about choosing between 'legacy' Armadillo 14.6.3 or 'current' Armadillo 15.0.1 or later.

### Performance

Performance is excellent, and on par or exceeding the performance of another package (or two,
following a renaming) claiming otherwise. See [this post][benchmark] and its [repo][ldlasb2] for a
detailed debunking of that claim.

### Documentation

The package contains a pdf vignette which is a pre-print of the [paper by Eddelbuettel and
Sanderson][rcpparmapaper] in CSDA (2014), as well as an introductory vignette for the sparse matrix
conversions.

### Installation

RcppArmadillo is a [CRAN package][cran pkg], and lives otherwise in its own habitat on
[GitHub](https://github.com/RcppCore/RcppArmadillo) within the
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


[armadillo]: https://arma.sourceforge.net
[rcpp]: https://www.rcpp.org
[cran]: https://cran.r-project.org
[cran pkg]: https://cran.r-project.org/package=RcppArmadillo
[benchmark]: https://eddelbuettel.github.io/ldlasb2/benchmarks.html
[ldlasb2]: https://github.com/eddelbuettel/ldlasb2
[rcpparmapaper]: http://dx.doi.org/10.1016/j.csda.2013.02.005
