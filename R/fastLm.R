## fastLm.R: Rcpp/Armadillo implementation of lm()
##
## Copyright (C)  2010 - 2021  Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

fastLmPure <- function(X, y) {

    stopifnot(is.matrix(X), is.numeric(y), nrow(y)==nrow(X))

    .Call(`_RcppArmadillo_fastLm_impl`, X, y)
}

fastLm <- function(X, ...) UseMethod("fastLm")

fastLm.default <- function(X, y, ...) {

    X <- as.matrix(X)
    y <- as.numeric(y)

    res <- fastLmPure(X, y)

    res$coefficients <- as.vector(res$coefficient)

    names(res$coefficients) <- colnames(X)

    res$fitted.values <- as.vector(X %*% res$coefficients)
    res$residuals <- y - res$fitted.values
    res$call <- match.call()
    res$intercept <- any(apply(X, 2, function(x) all(x == x[1])))

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

    ## cf src/library/stats/R/lm.R and case with no weights and an intercept
    f <- object$fitted.values
    r <- object$residuals
    #mss <- sum((f - mean(f))^2)
    mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)
    rss <- sum(r^2)

    r.squared <- mss/(mss + rss)
    df.int <- if (object$intercept) 1L else 0L

    n <- length(f)
    rdf <- object$df
    adj.r.squared <- 1 - (1 - r.squared) * ((n - df.int)/rdf)

    res <- list(call=object$call,
                coefficients=TAB,
                r.squared=r.squared,
                adj.r.squared=adj.r.squared,
                sigma=sqrt(sum((object$residuals)^2)/rdf),
                df=object$df,
                residSum=summary(object$residuals, digits=5)[-4])

    class(res) <- "summary.fastLm"
    res
}

print.summary.fastLm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nResiduals:\n")
    print(x$residSum)
    cat("\n")

    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
    digits <- max(3, getOption("digits") - 3)
    cat("\nResidual standard error: ", formatC(x$sigma, digits=digits), " on ",
        formatC(x$df), " degrees of freedom\n", sep="")
    cat("Multiple R-squared: ", formatC(x$r.squared, digits=digits),
        ",\tAdjusted R-squared: ",formatC(x$adj.r.squared, digits=digits),
        "\n", sep="")
    invisible(x)
}

fastLm.formula <- function(formula, data=list(), ...) {
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)

    res <- fastLm.default(x, y, ...)
    res$call <- match.call()
    res$formula <- formula
    res$intercept <- attr(attr(mf, "terms"), "intercept")
    res
}

predict.fastLm <- function(object, newdata=NULL, ...) {
    if (is.null(newdata)) {
        y <- fitted(object)
    } else {
        if (!is.null(object$formula)) {
            x <- model.matrix(object$formula, newdata)
        } else {
            x <- newdata 						# #nocov
        }
        y <- as.vector(x %*% coef(object))
    }
    y
}
