## .SciPy2R.R: Conversion of SciPy sparse matrix to R
##
## Copyright (C) 2017 - 2020  Binxiang Ni and Dirk Eddelbuettel
##
## This file is part of RcppArmadillo.
##
## RcppArmadillo is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppArmadillo is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

.SciPy2R <- function(spmat) {           # #nocov start
    .Deprecated("This function is no longer needed as 'reticulate' now converts from 'scipy'.",
                package="reticulate")

    if (!requireNamespace("reticulate", quietly=TRUE)) {
        stop("You must install the 'reticulate' package (and have SciPy).", call.=FALSE)
    }

    type <- spmat$getformat()
    shape <- unlist(spmat$shape)
    data <- as.vector(reticulate::py_get_attr(spmat, "data"))

    # Since CSC Matrix from SciPy has been automatically converted to dgCMatrix,
    # the conversion for csc matrix is no longer needed.

    # if (type == "csc") {
    #     indices <- as.vector(reticulate::py_get_attr(spmat, "indices"))
    #     indptr <- as.vector(reticulate::py_get_attr(spmat, "indptr"))
    #     res <- new("dgCMatrix", i = indices, p = indptr, x = data, Dim = shape)
    # }
    if (type == "coo") {
        row <- as.vector(reticulate::py_get_attr(spmat, "row"))
        col <- as.vector(reticulate::py_get_attr(spmat, "col"))
        res <- new("dgTMatrix", i = row, j = col, x = data, Dim = shape)
    } else if (type == "csr") {
        indices <- as.vector(reticulate::py_get_attr(spmat, "indices"))
        indptr <- as.vector(reticulate::py_get_attr(spmat, "indptr"))
        res <- new("dgRMatrix", j = indices, p = indptr, x = data, Dim = shape)
    } else {
        stop("Only CSC, COO and CSR matrices from SciPy are supported.", call.=FALSE)
    }
    return(res)
}										# #nocov end
