#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP RcppArmadillo_armadillo_set_seed(SEXP);
extern SEXP RcppArmadillo_armadillo_set_seed_random();
extern SEXP RcppArmadillo_armadillo_version(SEXP);
extern SEXP RcppArmadillo_fastLm(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"RcppArmadillo_armadillo_set_seed",        (DL_FUNC) &RcppArmadillo_armadillo_set_seed,        1},
    {"RcppArmadillo_armadillo_set_seed_random", (DL_FUNC) &RcppArmadillo_armadillo_set_seed_random, 0},
    {"RcppArmadillo_armadillo_version",         (DL_FUNC) &RcppArmadillo_armadillo_version,         1},
    {"RcppArmadillo_fastLm",                    (DL_FUNC) &RcppArmadillo_fastLm,                    2},
    {NULL, NULL, 0}
};

void R_init_RcppArmadillo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
