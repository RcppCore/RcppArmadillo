# Copyright (C) 2010 - 2013  Dirk Eddelbuettel, Romain Francois and Douglas Bates
# Copyright (C) 2014 - 2017  Dirk Eddelbuettel
# Earlier copyrights Gregor Gorjanc, Martin Maechler and Murray Stokely as detailed below
#
# This file is part of RcppArmadillo.
#
# RcppArmadillo is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 2 of the
# License, or (at your option) any later version.
#
# RcppArmadillo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

## doRUnit.R --- Run RUnit tests
##
## with credits to package fUtilities in RMetrics
## which credits Gregor Gojanc's example in CRAN package  'gdata'
## as per the R Wiki http://wiki.r-project.org/rwiki/doku.php?id=developers:runit
## and changed further by Martin Maechler
## and more changes by Murray Stokely in HistogramTools
## and then used adapted in RProtoBuf
## and now used in Rcpp and here
##
## Dirk Eddelbuettel, Feb - June 2014

if (requireNamespace("RUnit", quietly=TRUE) &&
    requireNamespace("RcppArmadillo", quietly=TRUE)) {

    library(RUnit)
    library(RcppArmadillo)

    ## Define tests
    testSuite <- defineTestSuite(name="RcppArmadillo Unit Tests",
                                 dirs=system.file("unitTests", package = "RcppArmadillo"),
                                 testFuncRegexp = "^[Tt]est.+")

    Sys.setenv("R_TESTS"="")	    	# without this, we get (or used to get) unit test failures

    tests <- runTestSuite(testSuite)    # run tests
    printTextProtocol(tests)		# print results

    ## Return success or failure to R CMD CHECK
    if (getErrors(tests)$nFail > 0) stop("TEST FAILED!")
    if (getErrors(tests)$nErr > 0) stop("TEST HAD ERRORS!")
    if (getErrors(tests)$nTestFunc < 1) stop("NO TEST FUNCTIONS RUN!")
}

