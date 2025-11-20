
##  RcppArmadillo configure.ac
##
##  'Rcpp' Integration for the 'Armadillo' Templated Linear Algebra Library
##
##  Copyright (C) 2016 - 2025  Dirk Eddelbuettel
##
##  Licensed under GPL-2 or later

## require at least autoconf 2.69
AC_PREREQ([2.69])

## Process this file with autoconf to produce a configure script.
AC_INIT([RcppArmadillo],[15.2.2-0],[edd@debian.org])

## Set R_HOME, respecting an environment variable if one is set
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi

## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

## Is R already configured to compile things using OpenMP without
## any extra hand-holding?
#openmp_already_works="no"

## default to not even thinking about OpenMP as Armadillo wants a pragma
## variant available if and only if C++11 is used with g++ 5.4 or newer
#can_use_openmp="no"

## Ensure TMPDIR is set.
AC_MSG_CHECKING([whether we have a suitable tempdir])
TMPDIR=$("${R_HOME}/bin/R" --vanilla --slave -e "cat(dirname(tempdir()))")
AC_MSG_RESULT([${TMPDIR}])

## Check if R is configured to compile programs using OpenMP out-of-the-box.
AC_MSG_CHECKING([whether R CMD SHLIB can already compile programs using OpenMP])

## Create private directory in TMPDIR.
BUILDDIR="${TMPDIR}/rcpparmadillo-$$-$RANDOM"
mkdir -p "${BUILDDIR}"

owd=$(pwd)
cd "${BUILDDIR}"

cat <<EOF > test-omp.cpp
#include <omp.h>
int main() {
  return omp_get_num_threads();
}
EOF

## Execute R CMD SHLIB.
"${R_HOME}/bin/R" CMD SHLIB test-omp.cpp >/dev/null 2>&1
if test x"$?" = x"0"; then
    AC_MSG_RESULT([yes])
    openmp_already_works="yes"
    arma_have_openmp="#define ARMA_USE_OPENMP 1"
    can_use_openmp="yes"
else
    AC_MSG_RESULT([no])
    openmp_already_works="no"
    arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"
    can_use_openmp="no"
fi

## Go back home.
cd "${owd}"
rm -rf "${BUILDDIR}"

## Additional Apple check
AC_MSG_CHECKING([for macOS])
RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
if test x"${RSysinfoName}" = x"Darwin"; then
    AC_MSG_RESULT([found])
    AC_MSG_CHECKING([for macOS Apple compiler])
    apple_compiler=$($CXX --version 2>&1 | grep -i -c -e 'apple llvm')
    if test x"${apple_compiler}" = x"1"; then
        AC_MSG_RESULT([found])
        AC_MSG_WARN([OpenMP unavailable and turned off.])
        can_use_openmp="no"
    fi
else
    AC_MSG_RESULT([not found])
fi

if test x"${can_use_openmp}" = x"yes"; then
    AC_MSG_CHECKING([for OpenMP])
    ## if R has -fopenmp we should be good
    allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
    hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)
    if test x"${hasOpenMP}" = x""; then
        AC_MSG_RESULT([missing])
        arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"
        openmp_flag=""
    else
        AC_MSG_RESULT([found and suitable])
        arma_have_openmp="#define ARMA_USE_OPENMP 1"
        openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
    fi
fi

## now use all these
AC_SUBST([ARMA_HAVE_OPENMP], ["${arma_have_openmp}"])
AC_SUBST([OPENMP_FLAG], ["${openmp_flag}"])
AC_CONFIG_FILES([inst/include/RcppArmadillo/config/RcppArmadilloConfigGenerated.h src/Makevars])
AC_OUTPUT
