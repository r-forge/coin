
.LinStatExpCov1d <- function(X, Y, weights, subset, block, 
                             varonly = FALSE, B = 0L, standardise = FALSE, 
                             tol = sqrt(.Machine$double.eps)) 
{
    if (NROW(X) != NROW(Y))
        stop("dimensions of X and Y don't match")

    if (is.integer(X)) {
        if (is.null(attr(X, "levels")))
            attr(X, "levels") <- 1:max(X)
    }

    if (!missing(weights) && length(weights) > 0) {
        if (!((NROW(X) == length(weights)) && 
              is.integer(weights) &&
              all(weights >= 0)))
            stop("incorrect weights")
    } else {
        weights <- integer(0)
    }

    if (!missing(subset) && length(subset) > 0) {
        if (!((max(subset) <= NROW(X)) &&
              (min(subset) >= 1L) && 
              is.integer(subset)))
            stop("incorrect subset")
        if (min(subset) == 0) stop("subset has start 1 index")
        subset <- subset - 1L
    } else {
        subset <- integer(0)
    }

    if (!missing(block) && length(block) > 0) {
        if (!((NROW(X) == length(block)) &&
              is.factor(block)))
            stop("incorrect block")
    } else {
        block <- integer(0)
    }

    ms <- !(complete.cases(X) & complete.cases(Y))
    if (any(ms)) {
        if (length(subset) > 0) {
            subset <- subset[!(subset %in% (which(ms) - 1L))] 
        } else {
            subset <- (0:(NROW(X) - 1))[-which(ms)]
        }
    }
    
    ret <- .Call("R_ExpectationCovarianceStatistic", X, Y, weights, subset, 
                 block, as.integer(varonly), as.double(tol), PACKAGE = "libcoin")
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    if (B > 0)
        ret$PermutedLinearStatistic <- .Call("R_PermutedLinearStatistic", ret, X, Y, weights, 
                         subset, block, as.integer(B), as.integer(standardise), 
                         PACKAGE = "libcoin")
    ret
}

.LinStatExpCov2d <- function(X, Y, ix, iy, weights, subset, block, 
                             varonly = FALSE, B = 0,
                             standardise = FALSE, 
                             tol = sqrt(.Machine$double.eps)) 
{
    if (!((length(ix) == length(iy)) &&
          is.integer(ix) && is.integer(iy)))
        stop("incorrect ix and/or iy")

    if (is.null(attr(ix, "levels")))
        attr(ix, "levels") <- 1:max(ix)
    if (is.null(attr(iy, "levels")))
        attr(iy, "levels") <- 1:max(iy)

    if (!missing(X) && length(X) > 0) {
        if (!((min(ix) >= 0 && nrow(X) == (max(ix) + 1)) &&
              all(complete.cases(X)) &&
              (nrow(X) == (length(attr(ix, "levels")) + 1))))
            stop("incorrect X")
    } else  {
        X <- numeric(0)
    }

    if (!(all(complete.cases(Y))) &&
          (nrow(Y) == (length(attr(iy, "levels")) + 1)) &&
          (min(iy) >= 0 && nrow(Y) == (max(iy) + 1)))
        stop("incorrect Y")

    if (!missing(weights) && length(weights) > 0) {
        if (!((length(ix) == length(weights)) && 
              is.integer(weights) &&
              all(weights >= 0)))
            stop("incorrect weights")
    } else {
        weights <- integer(0)
    }

    if (!missing(subset) && length(subset) > 0) {
        if (!((max(subset) <= length(ix)) && 
              (min(subset) >= 1L) &&
              is.integer(subset)))
            stop("incorrect subset")
        if (min(subset) == 0) stop("subset has start 1 index")
        subset <- subset - 1L
    } else {
        subset <- integer(0)
    }

    if (!missing(block) && length(block) > 0) {
        if (!((length(ix) == length(block)) &&
              is.factor(block)))
            stop("incorrect block")
    } else {
        block <- integer(0)
    }

    ret <- .Call("R_ExpectationCovarianceStatistic_2d", X, ix, Y, iy, 
        weights, subset, block, as.integer(varonly), as.double(tol), 
        PACKAGE = "libcoin")
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    if (B > 0)
        ret$PermutedLinearStatistic <- .Call("R_PermutedLinearStatistic_2d", ret, X, ix, Y, iy, 
                         block, as.integer(B), as.integer(standardise), 
                         PACKAGE = "libcoin")
    ret
}

LinStatExpCov <- function(X, Y, ix = NULL, iy = NULL, weights, subset, block, 
                          varonly = FALSE, B = 0, standardise = FALSE, 
                          tol = sqrt(.Machine$double.eps)) 
{

    if (is.null(ix) & is.null(iy))
        return(.LinStatExpCov1d(X = X, Y = Y, weights = weights, 
                                subset = subset, block = block, 
                                varonly = varonly, B = B, 
                                standardise = standardise, tol = tol))

    if (!is.null(ix) & !is.null(iy))
        return(.LinStatExpCov2d(X = X, Y = Y, ix = ix, iy = iy, 
                                weights = weights, subset = subset,
                                block = block, varonly = varonly, B = B, 
                                standardise = standardise, tol = tol))

    if (missing(X) & !is.null(ix))
        return(.LinStatExpCov1d(X = ix, Y = Y, weights = weights, 
                                subset = subset, block = block, 
                                varonly = varonly, B = B, 
                                standardise = standardise, tol = tol))

    stop("incorrect call to LinStatExpCov")
}


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

    if (!pvalue & (NCOL(object$PermutedLinearStatistic) > 0)) {
        object$PermutedLinearStatistic <- matrix(nrow = 0, ncol = 0)
        storage.mode(object$PermutedLinearStatistic) <- "double"
    }

    if (!object$Xfactor) {
        if (teststat == "quadratic") {
            ret <- .Call("R_ChisqTest", object, 
                         as.integer(pvalue), as.integer(lower), as.integer(log), 
                         PACKAGE = "libcoin")
        } else {
            ret <- .Call("R_MaxtypeTest", object,  
                         as.integer(alt), as.integer(pvalue), as.integer(lower), 
                         as.integer(log), as.integer(pargs$maxpts), 
                         as.double(pargs$releps), as.double(pargs$abseps), 
                         PACKAGE = "libcoin")
            if (teststat == "scalar") {
                var <- ifelse(object$varonly, object$Variance, object$Covariance)
                ret$TestStatistic <- object$LinearStatistic - object$Expectation
                ret$TestStatistic <- ifelse(var > object$tol, ret$TestStatistic / sqrt(var), NaN)
            }
        }
    } else {
        ret <- .Call("R_MaxSelectTest", object, as.integer(ordered), 
                     as.integer(test), as.integer(minbucket), 
                     as.integer(lower), as.integer(log), PACKAGE = "libcoin")
    }
    ret
}
