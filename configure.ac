
##  RcppArmadillo configure.ac
##
##  'Rcpp' Integration for the 'Armadillo' Templated Linear Algebra Library
##
##  Copyright (C) 2016 - 2023  Dirk Eddelbuettel
##
##  Licensed under GPL-2 or later

## require at least autoconf 2.69
AC_PREREQ([2.69])

## Process this file with autoconf to produce a configure script.
AC_INIT([RcppArmadillo],[0.12.4.0.0],[edd@debian.org])

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
openmp_already_works="no"

## default to not even thinking about OpenMP as Armadillo wants a pragma
## variant available if and only if C++11 is used with g++ 5.4 or newer
can_use_openmp="no"

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
else
    AC_MSG_RESULT([no])
fi

## Go back home.
cd "${owd}"
rm -rf "${BUILDDIR}"

## If the above checks failed, then perform other heuristics
## based on the compiler version, etc.
if test x"${openmp_already_works}" = x"no"; then

    ## Check the C++ compiler using the CXX value set

    ## If it is g++, we have GXX set so let's examine it
    if test "${GXX}" = yes; then
        AC_MSG_CHECKING([whether g++ version is sufficient])
        gxx_version=$(${CXX} -v 2>&1 | awk '/^.*g.. version/ {print $3}')
        case ${gxx_version} in
            1.*|2.*|3.*|4.0.*|4.1.*|4.2.*|4.3.*|4.4.*|4.5.*|4.6.*|4.7.0|4.7.1)
                 AC_MSG_RESULT([no])
                 AC_MSG_WARN([Only g++ version 4.7.2 or greater can be used with RcppArmadillo.])
                 AC_MSG_ERROR([Please use a different compiler.])
            ;;
            4.7.*|4.8.*|4.9.*|5.0*|5.1*|5.2*|5.3*)
                 AC_MSG_RESULT([yes, but without OpenMP as version ${gxx_version} (Armadillo constraint)])
                 ## we know this one is bad
                 can_use_openmp="no"
            ;;
            5.4*|5.5*|5.6*|5.7*|5.8*|5.9*|6.*|7.*|8.*|9.*|10.*|11.*|12.*)
                 AC_MSG_RESULT([yes, with OpenMP as version ${gxx_version}])
                 ## we know this one is good, yay
                 can_use_openmp="yes"
            ;;
            *)
                 AC_MSG_RESULT([almost])
                 AC_MSG_WARN([Compiler self-identifies as being compliant with GNUC extensions but is not g++.])
                 ## we know nothing, so no
                 can_use_openmp="no"
            ;;
        esac
    fi

    ## Check for Apple LLVM

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
        else
            AC_MSG_RESULT([not found])
            AC_MSG_CHECKING([for clang compiler])
            clang_compiler=$($CXX --version 2>&1 | grep -i -c -e 'clang ')

            if test x"${clang_compiler}" = x"1"; then
                AC_MSG_RESULT([found])
                AC_MSG_CHECKING([for OpenMP compatible version of clang])
                clang_version=$(${CXX} -v 2>&1 | awk '/^.*clang version/ {print $3}')

                case ${clang_version} in
                    4.*|5.*|6.*|7.*|8.*|9.*|10.*|11.*)
                        AC_MSG_RESULT([found and suitable])
                        can_use_openmp="yes"
                    ;;
                    *)
                        AC_MSG_RESULT([not found])
                        AC_MSG_WARN([OpenMP unavailable and turned off.])
                        can_use_openmp="no"
                    ;;
                esac
            else
                AC_MSG_RESULT([not found])
                AC_MSG_WARN([unsupported macOS build detected; if anything breaks, you keep the pieces.])
            fi
        fi
    else
        AC_MSG_RESULT([no])
    fi

fi # if test x"${openmp_already_works}" = x"no"

## Check for suitable LAPACK_LIBS
AC_MSG_CHECKING([LAPACK_LIBS])

## external LAPACK has the required function
lapack=$(${R_HOME}/bin/R CMD config LAPACK_LIBS)
hasRlapack=$(echo ${lapack} | grep lRlapack)

## in what follows below we substitute both side of the define/undef
## while this may seem a little unusual we do it to fully reproduce the
## previous bash-based implementation

if test x"${hasRlapack}" = x""; then
    ## We are using a full Lapack and can use zgbsv -- so #undef remains
    AC_MSG_RESULT([system LAPACK found])
    arma_lapack="#undef ARMA_CRIPPLED_LAPACK"
else
    ## We are using R's subset of Lapack and CANNOT use zgbsv etc, so we mark it
    AC_MSG_RESULT([R-supplied partial LAPACK found])
    AC_MSG_WARN([Some complex-valued LAPACK functions may not be available])
    arma_lapack="#define ARMA_CRIPPLED_LAPACK 1"
fi

## Default the OpenMP flag to the empty string.
## If and only if OpenMP is found, expand to $(SHLIB_OPENMP_CXXFLAGS)
openmp_flag=""

## Set the fallback, by default it is nope
arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"

if test x"${openmp_already_works}" = x"yes"; then
    arma_have_openmp="#define ARMA_USE_OPENMP 1"
    openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
fi

if test x"${can_use_openmp}" = x"yes"; then
    AC_MSG_CHECKING([for OpenMP])
    ## if R has -fopenmp we should be good
    allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
    hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)
    if test x"${hasOpenMP}" = x""; then
        AC_MSG_RESULT([missing])
        arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"
    else
        AC_MSG_RESULT([found and suitable])
        arma_have_openmp="#define ARMA_USE_OPENMP 1"
        openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
    fi
fi


## now use all these
AC_SUBST([ARMA_LAPACK],["${arma_lapack}"])
AC_SUBST([ARMA_HAVE_OPENMP], ["${arma_have_openmp}"])
AC_SUBST([OPENMP_FLAG], ["${openmp_flag}"])
AC_CONFIG_FILES([inst/include/RcppArmadillo/config/RcppArmadilloConfigGenerated.h src/Makevars])
AC_OUTPUT
