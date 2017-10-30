
# R Header

###
### TO NOT EDIT THIS FILE
### 
### Edit `libcoin.w' and run `nuweb -r libcoin.w'
###

# LinStatExpCov

LinStatExpCov <- function(X, Y, ix = NULL, iy = NULL, weights = integer(0),
                          subset = integer(0), block = integer(0),
                          varonly = FALSE, nperm = 0, standardise = FALSE,
                          tol = sqrt(.Machine$double.eps))
{
    if (is.null(ix) & is.null(iy))
        return(.LinStatExpCov1d(X = X, Y = Y, weights = weights,
                                subset = subset, block = block,
                                varonly = varonly, nperm = nperm,
                                standardise = standardise, tol = tol))

    if (!is.null(ix) & !is.null(iy))
        return(.LinStatExpCov2d(X = X, Y = Y, ix = ix, iy = iy,
                                weights = weights, subset = subset,
                                block = block, varonly = varonly, nperm = nperm,
                                standardise = standardise, tol = tol))

    if (missing(X) & !is.null(ix))
        return(.LinStatExpCov1d(X = ix, Y = Y, weights = weights,
                                subset = subset, block = block,
                                varonly = varonly, nperm = nperm,
                                standardise = standardise, tol = tol))

    stop("incorrect call to LinStatExpCov")
}

# LinStatExpCov1d

.LinStatExpCov1d <- function(X, Y, weights = integer(0), subset = integer(0), block = integer(0),
                             varonly = FALSE, nperm = 0, standardise = FALSE,
                             tol = sqrt(.Machine$double.eps))
{

    if (NROW(X) != NROW(Y))
        stop("dimensions of X and Y don't match")
    N <- NROW(X)

    if (is.integer(X)) {
        if (is.null(attr(X, "levels")))
            attr(X, "levels") <- 1:max(X)
    }

    # Check weights, subset, block
    
    if (length(weights) > 0) {
        if (!((N == length(weights)) && all(weights >= 0)))
            stop("incorrect weights")
    }

    if (length(subset) > 0) {
        rs <- range(subset)
        if (!((rs[2] <= N) && (rs[1] >= 1L)))
            stop("incorrect subset")
    }

    if (length(block) > 0) {
        if (!((N == length(block)) && is.factor(block)))
            stop("incorrect block")
    }
    

    ms <- !(complete.cases(X) & complete.cases(Y))
    if (all(ms))
        stop("all observations are missing")
    if (any(ms)) {
        if (length(subset) > 0) {
            if (all(subset %in% which(ms)))
                stop("all observations are missing")
            subset <- subset[!(subset %in% which(ms))]
        } else {
            subset <- (1:N)[-which(ms)]
        }
    }

    ret <- .Call(R_ExpectationCovarianceStatistic, X, Y, weights, subset,
                 block, as.integer(varonly), as.double(tol))
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    if (nperm > 0) {
        if (standardise) {
            standardise <- ret
        } else {
            standardise <- integer(0)
        }
        ret$PermutedLinearStatistic <-
            .Call(R_PermutedLinearStatistic, X, Y, weights, subset,
                  block, as.double(nperm), standardise)
    }
    class(ret) <- c("LinStatExpCov1d", "LinStatExpCov")
    ret
}

# LinStatExpCov2d

.LinStatExpCov2d <- function(X = numeric(0), Y, ix, iy, weights = integer(0), subset = integer(0),
                             block = integer(0), varonly = FALSE, nperm = 0,
                             standardise = FALSE,
                             tol = sqrt(.Machine$double.eps))
{
    if (!((length(ix) == length(iy)) &&
          is.integer(ix) && is.integer(iy)))
        stop("incorrect ix and/or iy")
    N <- length(ix)

    if (is.null(attr(ix, "levels")))
        attr(ix, "levels") <- 1:max(ix)
    if (is.null(attr(iy, "levels")))
        attr(iy, "levels") <- 1:max(iy)

    if (length(X) > 0) {
        if (!((min(ix) >= 0 && nrow(X) == (length(attr(ix, "levels")) + 1)) &&
              all(complete.cases(X)) &&
              (nrow(X) == (length(attr(ix, "levels")) + 1))))
            stop("incorrect X")
    }

    if (!(all(complete.cases(Y))) &&
          (nrow(Y) == (length(attr(iy, "levels")) + 1)) &&
          (min(iy) >= 0L && nrow(Y) == (length(attr(iy, "levels")) + 1)))
        stop("incorrect Y")

    # Check weights, subset, block
    
    if (length(weights) > 0) {
        if (!((N == length(weights)) && all(weights >= 0)))
            stop("incorrect weights")
    }

    if (length(subset) > 0) {
        rs <- range(subset)
        if (!((rs[2] <= N) && (rs[1] >= 1L)))
            stop("incorrect subset")
    }

    if (length(block) > 0) {
        if (!((N == length(block)) && is.factor(block)))
            stop("incorrect block")
    }
    

    ret <- .Call(R_ExpectationCovarianceStatistic_2d, X, ix, Y, iy,
                 weights, subset, block, as.integer(varonly), as.double(tol))
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    if (nperm > 0) {
        if (standardise) {
            standardise <- ret
        } else {
            standardise <- integer(0)
        }
        ret$PermutedLinearStatistic <-
            .Call(R_PermutedLinearStatistic_2d, X, ix, Y, iy, weights, subset,
                  block, nperm, ret$Table, standardise)
    }
    class(ret) <- c("LinStatExpCov2d", "LinStatExpCov")
    ret
}

# vcov LinStatExpCov

vcov.LinStatExpCov <- function(object, ...) {
    if (object$varonly)
        stop("cannot extract covariance matrix")
    PQ <- prod(object$dim)
    ret <- matrix(0, nrow = PQ, ncol = PQ)
    ret[lower.tri(ret, diag = TRUE)] <- object$Covariance
    ret <- ret + t(ret)
    diag(ret) <- diag(ret) / 2
    ret
}


# doTest

### note: lower = FALSE => p-value; lower = TRUE => 1 - p-value
doTest <- function(object, teststat = c("maximum", "quadratic", "scalar"),
                     alternative = c("two.sided", "less", "greater"),
                     pvalue = TRUE, lower = FALSE, log = FALSE,
                     minbucket = 10L, ordered = TRUE, pargs = GenzBretz())
{

    ### avoid match.arg for performance reasons
    teststat <- teststat[1]
    if (!any(teststat == c("maximum", "quadratic", "scalar")))
        stop("incorrect teststat")
    alternative <- alternative[1]
    if (!any(alternative == c("two.sided", "less", "greater")))
        stop("incorrect alternative")

    if (teststat == "quadratic") {
        if (alternative != "two.sided")
            stop("incorrect alternative")
    }

    test <- which(c("maximum", "quadratic", "scalar") == teststat)
    if (test == 3) {
        if (length(object$LinearStatistic) != 1)
            stop("scalar test statistic not applicable")
        test <- 1L ### scalar is maximum internally
    }
    alt <- which(c("two.sided", "less", "greater") == alternative)

    if (!pvalue & (NCOL(object$PermutedLinearStatistic) > 0))
        object$PermutedLinearStatistic <- matrix(NA_real_, nrow = 0, ncol = 0)

    if (!object$Xfactor) {
        if (teststat == "quadratic") {
            ret <- .Call(R_QuadraticTest, object,
                         as.integer(pvalue), as.integer(lower), as.integer(log))
        } else {
            ret <- .Call(R_MaximumTest, object,
                         as.integer(alt), as.integer(pvalue), as.integer(lower),
                         as.integer(log), as.integer(pargs$maxpts),
                         as.double(pargs$releps), as.double(pargs$abseps))
            if (teststat == "scalar") {
                var <- if (object$varonly) object$Variance else object$Covariance
                ret$TestStatistic <- object$LinearStatistic - object$Expectation
                ret$TestStatistic <-
                    if (var > object$tol) ret$TestStatistic / sqrt(var) else NaN
            }
        }
    } else {
        ret <- .Call(R_MaximallySelectedTest, object, as.integer(ordered),
                     as.integer(test), as.integer(minbucket),
                     as.integer(lower), as.integer(log))
    }
    ret
}

.libcoinCall <- function(FUN, ...) .Call(FUN, ...)
