
if (requireNamespace("tinytest", quietly=TRUE) &&
    utils::packageVersion("tinytest") >= "1.0.0") {

    ## Set a seed to make the test deterministic
    set.seed(42)

    ## there are several more granular ways to test files in a tinytest directory,
    ## see its package vignette; tests can also run once the package is installed
    ## using the same command `test_package(pkgName)`, or by director or file
    ## parked, now:  , getOption("Ncpus", 1)
    tinytest::test_package("RcppArmadillo", side_effects=TRUE)
}
