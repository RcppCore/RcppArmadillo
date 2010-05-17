
## fastLm.R: Rcpp/Armadillo implementation of lm()
##
## Copyright (C)  2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

fastLmPure <- function(y, X) {

    stopifnot(is.matrix(X))
    stopifnot(nrow(y)==nrow(X))

    res <- .Call("fastLm", y, X, package="RcppArmadillo")
}

fastLm <- function(x, ...) UseMethod("fastLm")

fastLm.default <- function(x, y, ...) {

    x <- as.matrix(x)
    y <- as.numeric(y)

    res <- fastLmPure(y, x)

    res$coefficients <- res$coefficient[,1] # force into single-col vector

    names(res$coefficients) <- colnames(x)

    res$fitted.values <- as.vector(x %*% res$coefficients)
    res$residuals <- y - res$fitted.values
    res$call <- match.call()

    class(res) <- "fastLm"
    res
}

print.fastLm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients, digits=5)
}

summary.fastLm <- function(object, ...) {
    se <- object$stderr
    tval <- coef(object)/se

    TAB <- cbind(Estimate = coef(object),
                 StdErr = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval), df=object$df))

    # why do I need this here?
    rownames(TAB) <- names(object$coefficients)
    colnames(TAB) <- c("Estimate", "StdErr", "t.value", "p.value")

    res <- list(call=object$call,
                coefficients=TAB)

    class(res) <- "summary.fastLm"
    res
}

print.summary.fastLm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")

    printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

fastLm.formula <- function(formula, data=list(), ...) {
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)

    res <- fastLm.default(x, y, ...)
    res$call <- match.call()
    res$formula <- formula
    res
}

predict.fastLm <- function(object, newdata=NULL, ...) {
    if (is.null(newdata)) {
        y <- fitted(object)
    } else {
        if (!is.null(object$formula)) {
            x <- model.matrix(object$formula, newdata)
        } else {
            x <- newdata
        }
        y <- as.vector(x %*% coef(object))
    }
    y
}
