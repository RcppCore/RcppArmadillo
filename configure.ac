
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

AC_MSG_CHECKING([what system we are on])
#SYSKERNEL=$(uname -s)
#SYSMACHINE=$(uname -m)
#SYSOS=$(uname -o)
#AC_MSG_RESULT([${SYSKERNEL} ${SYSMACHINE} ${SYSOS}])
## Rely on Sys.info() from R for cross-compilation cases
SYSKERNEL=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
SYSMACHINE=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["machine"]])')
AC_MSG_RESULT([running ${SYSKERNEL} on ${SYSMACHINE}])


## Is R already configured to compile things using OpenMP without
## any extra hand-holding?
#openmp_already_works="no"

## default to not even thinking about OpenMP as Armadillo wants a pragma
## variant available if and only if C++11 is used with g++ 5.4 or newer
can_use_openmp="no"

## Ensure TMPDIR is set.
AC_MSG_CHECKING([whether we have a suitable tempdir])
TMPDIR=$("${R_HOME}/bin/R" --vanilla --slave -e "cat(dirname(tempdir()))")
AC_MSG_RESULT([${TMPDIR}])

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

if test x"${SYSKERNEL}" = x"Linux"; then

   ## Check if R is configured to compile OpenMP programs out-of-the-box.
   if test x"${can_use_openmp}" = x"no"; then
       AC_MSG_CHECKING([whether R CMD SHLIB can already compile OpenMP programs])

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
   fi

   ## If needed, check if R is configured to compile OpenMP programs using -fopenmp
   if test x"${can_use_openmp}" = x"no"; then
       AC_MSG_CHECKING([whether R CMD SHLIB can compile OpenMP via -fopenmp])

       PKG_CXXFLAGS="${PKG_CXXFLAGS} -fopenmp" PKG_LIBS="${PKG_LIBS} -fopen" "${R_HOME}/bin/R" CMD SHLIB -fopenmp test-omp.cpp >/dev/null 2>&1
       if test x"$?" = x"0"; then
           AC_MSG_RESULT([yes])
           openmp_already_works="yes"
           arma_have_openmp="#define ARMA_USE_OPENMP 1"
           can_use_openmp="yes"
           # keep any entries user may have set
           PKG_CXXFLAGS="${PKG_CXXFLAGS} -fopenmp"
           PKG_LIBS="${PKG_LIBS} -fopenmp"
       else
           AC_MSG_RESULT([no])
           openmp_already_works="no"
           arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"
           can_use_openmp="no"
       fi
   fi
fi

if test x"${SYSKERNEL}" = x"Darwin"; then

   ## Check if R is configured to compile OpenMP programs using -Xclang -fopenmp
   if test x"${can_use_openmp}" = x"no"; then
       AC_MSG_CHECKING([whether R CMD SHLIB can compile OpenMP programs using '-Xclang -fopenmp'])

       PKG_CXXFLAGS="${PKG_CXXFLAGS} -Xclang -fopenmp" PKG_LIBS="${PKG_LIBS} -lomp" "${R_HOME}/bin/R" CMD SHLIB test-omp.cpp >/dev/null 2>&1
       if test x"$?" = x"0"; then
           AC_MSG_RESULT([yes])
           openmp_already_works="yes"
           arma_have_openmp="#define ARMA_USE_OPENMP 1"
           can_use_openmp="yes"
           # keep any entries user may have set
           PKG_CXXFLAGS="${PKG_CXXFLAGS} -Xclang -fopenmp"
           PKG_LIBS="${PKG_LIBS} -lomp"
       else
           AC_MSG_RESULT([no])
           openmp_already_works="no"
           arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"
           can_use_openmp="no"
       fi
    fi
fi

## Go back home.
cd "${owd}"
rm -rf "${BUILDDIR}"

dnl ## Additional Apple check
dnl AC_MSG_CHECKING([for macOS])
dnl RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
dnl if test x"${RSysinfoName}" = x"Darwin"; then
dnl     AC_MSG_RESULT([found])
dnl     AC_MSG_CHECKING([for macOS Apple compiler])
dnl     apple_compiler=$($CXX --version 2>&1 | grep -i -c -e 'apple llvm')
dnl     if test x"${apple_compiler}" = x"1"; then
dnl         AC_MSG_RESULT([found])
dnl         AC_MSG_WARN([OpenMP unavailable and turned off.])
dnl         can_use_openmp="no"
dnl     fi
dnl else
dnl     AC_MSG_RESULT([not found])
dnl fi

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
AC_SUBST([PKG_CXXFLAGS], ["${PKG_CXXFLAGS}"])
AC_SUBST([PKG_LIBS], ["${PKG_LIBS}"])
AC_CONFIG_FILES([inst/include/RcppArmadillo/config/RcppArmadilloConfigGenerated.h src/Makevars])
AC_OUTPUT
