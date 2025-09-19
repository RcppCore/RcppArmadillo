
// version/arma.h: Armadillo version selection
//
// Copyright (C)  2025-current  Dirk Eddelbuettel
//
// This file is part of RcppArmadillo.
//
// RcppArmadillo is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppArmadillo is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RcppArmadillo__version__arma__h
#define RcppArmadillo__version__arma__h


// Purpose: By explicitly defining or undefining a variable the 'defined(...) that
// is in RcppArmadilloForward.h works as expected allowing us to selectively include
// either a 'legacy' Armadillo (i.e. 14.6.3. the last version to a) allow C++11 and
// b) permit suppression of deprecation warnings) or a 'current' Armadillo (i.e. as
// of this writing 15.0.1, or any later version)
//
// See https://github.com/RcppCore/RcppArmadillo/issues/475 for more details

// Sanity check: Cannot select both current and legacy but only one
#if defined(ARMA_USE_CURRENT) && defined(ARMA_USE_LEGACY)
    #error "Do not select both 'current' and 'legacy', only one choice is possible."
#endif

// Sanity check: Cannot select current under C++11
#if defined(ARMA_USE_CURRENT) && __cplusplus <= 201103L
    #error "Do not select 'current' with C++11 (or older) as 'current' requires C++14 or newer."
#endif

// Carefully check and set if the user has -DARMA_USE_CURRENT
// This can be set in src/Makevars(.win) via PKG_CPPFLAGS (or equivalent)
#if defined(ARMA_USE_CURRENT)
#define ARMA_SELECTED_CURRENT_VERSION
#else
#undef ARMA_SELECTED_CURRENT_VERSION
#endif

// Carefully check and set if the user has -DARMA_USE_LEGACY
// This can be set in src/Makevars(.win) via PKG_CPPFLAGS (or equivalent)
#if defined(ARMA_USE_LEGACY)
#define ARMA_SELECTED_LEGACY_VERSION
#else
#undef ARMA_SELECTED_LEGACY_VERSION
#endif

// Fallback: Use 'legacy' and warn.
// This toggle can be turned to default to 'current' (allowing a legacy override) which we plan
// to do after a (sufficiently long) adjustment period
#if !defined(ARMA_SELECTED_LEGACY_VERSION) && !defined(ARMA_SELECTED_CURRENT_VERSION)
    // Show messages, adjusted for compilation standard (as we also want to move on from C++11)
    // We plan to 'at some point' flip to for C++14
    #if __cplusplus < 201402L
        #pragma message("Using fallback compilation with Armadillo 14.6.3. Please consider defining -DARMA_USE_CURRENT and also removing C++11 compilation directive. See GitHub issue #475 for more.")
        // Define selector used in RcppArmadilloForward.h
        #define ARMA_SELECTED_LEGACY_VERSION
        #undef ARMA_SELECTED_CURRENT_VERSION
    #endif
    // -- no longer automatically fall back to legacy version (unless in C++11 mode)
    // #else
    //     #pragma message("Using fallback compilation with Armadillo 14.6.3. Please consider defining -DARMA_USE_CURRENT. See GitHub issue #475 for more.")
    //     // Define selector used in RcppArmadilloForward.h
    //     // It is our intention to select current here after transition instead of legacy
    //     #define ARMA_SELECTED_LEGACY_VERSION
    //     #undef ARMA_SELECTED_CURRENT_VERSION
    // #endif
#endif


#endif
