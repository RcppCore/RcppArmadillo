
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
AC_INIT([RcppArmadillo],[15.2.3-1],[edd@debian.org])

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

## Rely on Sys.info() from R for cross-compilation cases
AC_MSG_CHECKING([what system we are on])
SYSKERNEL=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
SYSMACHINE=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["machine"]])')
AC_MSG_RESULT([running ${SYSKERNEL} on ${SYSMACHINE}])


## Default to not assuming OpenMP, but then test a variety of setups
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

AC_MSG_CHECKING([whether on Linux])
if test x"${SYSKERNEL}" = x"Linux"; then
   AC_MSG_RESULT([yes])

   ## Check if R is configured to compile OpenMP programs out-of-the-box.
   if test x"${can_use_openmp}" = x"no"; then
       AC_MSG_CHECKING([whether R CMD SHLIB can already compile OpenMP programs])
       "${R_HOME}/bin/R" CMD SHLIB test-omp.cpp >/dev/null 2>&1
       if test x"$?" = x"0"; then
           AC_MSG_RESULT([yes])
           can_use_openmp="yes"
       else
           AC_MSG_RESULT([no])
           can_use_openmp="no"
       fi
   fi

   ## If needed, check if R is configured to compile OpenMP programs using -fopenmp
   if test x"${can_use_openmp}" = x"no"; then
       AC_MSG_CHECKING([whether R CMD SHLIB can compile OpenMP via -fopenmp])
       PKG_CXXFLAGS="${PKG_CXXFLAGS} -fopenmp" PKG_LIBS="${PKG_LIBS} -fopenmp" "${R_HOME}/bin/R" CMD SHLIB -fopenmp test-omp.cpp >/dev/null 2>&1
       if test x"$?" = x"0"; then
           AC_MSG_RESULT([yes])
           can_use_openmp="yes"
           # keep any entries user may have set
           PKG_CXXFLAGS="${PKG_CXXFLAGS} -fopenmp"
           PKG_LIBS="${PKG_LIBS} -fopenmp"
       else
           AC_MSG_RESULT([no])
           can_use_openmp="no"
       fi
   fi
else
   AC_MSG_RESULT([no])
fi # if Linux

AC_MSG_CHECKING([whether on macOS])
if test x"${SYSKERNEL}" = x"Darwin"; then
   AC_MSG_RESULT([yes])

   ## Check if R is configured to compile OpenMP programs using -Xclang -fopenmp
   if test x"${can_use_openmp}" = x"no"; then
       AC_MSG_CHECKING([whether R CMD SHLIB can compile OpenMP programs using '-Xclang -fopenmp'])
       PKG_CXXFLAGS="${PKG_CXXFLAGS} -Xclang -fopenmp" PKG_LIBS="${PKG_LIBS} -lomp" "${R_HOME}/bin/R" CMD SHLIB test-omp.cpp >/dev/null 2>&1
       if test x"$?" = x"0"; then
           AC_MSG_RESULT([yes])
           can_use_openmp="yes"
           # keep any entries user may have set
           PKG_CXXFLAGS="${PKG_CXXFLAGS} -Xclang -fopenmp"
           PKG_LIBS="${PKG_LIBS} -lomp"
       else
           AC_MSG_RESULT([no])
           can_use_openmp="no"
       fi
   fi

   ## Check if R is configured to compile OpenMP programs using -fopenmp (cf data.table #6409)
   if test x"${can_use_openmp}" = x"no"; then
       AC_MSG_CHECKING([whether R CMD SHLIB can compile OpenMP programs using '-fopenmp'])
       PKG_CXXFLAGS="${PKG_CXXFLAGS} -fopenmp" PKG_LIBS="${PKG_LIBS} -fopenmp" "${R_HOME}/bin/R" CMD SHLIB test-omp.cpp >/dev/null 2>&1
       if test x"$?" = x"0"; then
           AC_MSG_RESULT([yes])
           can_use_openmp="yes"
           # keep any entries user may have set
           PKG_CXXFLAGS="${PKG_CXXFLAGS} -fopenmp"
           PKG_LIBS="${PKG_LIBS} -fopenmp"
       else
           AC_MSG_RESULT([no])
           can_use_openmp="no"
       fi
   fi

   if test x"${can_use_openmp}" = x"no"; then
       if test x"${SYSMACHINE}" = x"arm64"; then
           HOMEBREW_PREFIX=/opt/homebrew
       else
           HOMEBREW_PREFIX=/usr/local
       fi
       if test -e "${HOMEBREW_PREFIX}/opt/libomp"; then
           AC_MSG_CHECKING([whether R CMD SHLIB can use libomp at ${HOMEBREW_PREFIX}/opt/libomp])
           LIBOMP_INCLUDE="-I${HOMEBREW_PREFIX}/opt/libomp/include -Xclang -fopenmp"
           LIBOMP_LINK="-L${HOMEBREW_PREFIX}/opt/libomp/lib -lomp"
           PKG_CXXFLAGS="${PKG_CXXFLAGS} ${LIBOMP_INCLUDE}" PKG_LIBS="${PKG_LIBS} ${LIBOMP_LINK}" "${R_HOME}/bin/R" CMD SHLIB test-omp.cpp >/dev/null 2>&1
           if test x"$?" = x"0"; then
               AC_MSG_RESULT([yes])
               can_use_openmp="yes"
               # keep any entries user may have set
               PKG_CXXFLAGS="${PKG_CXXFLAGS} ${LIBOMP_INCLUDE}"
               PKG_LIBS="${PKG_LIBS} ${LIBOMP_LINK}"
           else
               AC_MSG_RESULT([no])
               can_use_openmp="no"
           fi
       fi
   fi
else
   AC_MSG_RESULT([no])
fi # if macOS

## Go back home.
cd "${owd}"
rm -rf "${BUILDDIR}"

# Overall summary
AC_MSG_CHECKING([for OpenMP])
if test x"${can_use_openmp}" = x"yes"; then
   AC_MSG_RESULT([found and suitable])
   arma_have_openmp="#define ARMA_USE_OPENMP 1"
   openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
else
   AC_MSG_RESULT([missing so no OpenMP acceleration])
   arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"
   openmp_flag=""
fi

## Now use all these
AC_SUBST([ARMA_HAVE_OPENMP], ["${arma_have_openmp}"])
AC_SUBST([OPENMP_FLAG], ["${openmp_flag}"])
AC_SUBST([PKG_CXXFLAGS], ["${PKG_CXXFLAGS}"])
AC_SUBST([PKG_LIBS], ["${PKG_LIBS}"])
AC_CONFIG_FILES([inst/include/RcppArmadillo/config/RcppArmadilloConfigGenerated.h src/Makevars])
AC_OUTPUT
