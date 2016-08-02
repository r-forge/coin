
.LinStatExpCov1d <- function(X, Y, weights, subset, block, 
                             varonly = FALSE, B = 0L, standardise = FALSE, 
                             tol = sqrt(.Machine$double.eps)) 
{
    stopifnot(NROW(X) == NROW(Y))

    if (is.integer(X)) {
        if (is.null(attr(X, "levels")))
            attr(X, "levels") <- 1:max(X)
    }

    if (!missing(weights) && length(weights) > 0) {
        stopifnot(NROW(X) == length(weights))
        stopifnot(is.integer(weights))
    } else {
        weights <- integer(0)
    }
    if (!missing(subset) && length(subset) > 0) {
        stopifnot(all(subset %in% 1:NROW(X)))
        stopifnot(is.integer(subset))
        if (min(subset) == 0) stop("subset has start 1 index")
        subset <- subset - 1L
    } else {
        subset <- integer(0)
    }
    if (!missing(block) && length(block) > 0) {
        stopifnot(NROW(X) == length(block))
        stopifnot(is.factor(block))
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
    stopifnot(length(ix) == length(iy))
    stopifnot(is.integer(ix))
    stopifnot(is.integer(iy))
    if (is.null(attr(ix, "levels")))
        attr(ix, "levels") <- 1:max(ix)
    if (is.null(attr(iy, "levels")))
        attr(iy, "levels") <- 1:max(iy)

    if (!missing(X) && length(X) > 0) {
        stopifnot(min(ix) >= 0 && nrow(X) == (max(ix) + 1))
        stopifnot(all(complete.cases(X)))
        stopifnot(nrow(X) == (length(attr(ix, "levels")) + 1))
    } else  {
        X <- numeric(0)
    }
    stopifnot(all(complete.cases(Y)))
    stopifnot(nrow(Y) == (length(attr(iy, "levels")) + 1))
    stopifnot(min(iy) >= 0 && nrow(Y) == (max(iy) + 1))

    if (!missing(weights) && length(weights) > 0) {
        stopifnot(length(ix) == length(weights))
        stopifnot(is.integer(weights))
    } else {
        weights <- integer(0)
    }
    if (!missing(subset) && length(subset) > 0) {
        stopifnot(all(subset %in% 1:length(ix)))
        stopifnot(is.integer(subset))
        if (min(subset) == 0) stop("subset has start 1 index")
        subset <- subset - 1L
    } else {
        subset <- integer(0)
    }
    if (!missing(block) && length(block) > 0) {
        stopifnot(length(ix) == length(block))
        stopifnot(is.factor(block))
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


### <FIXME> add alternative argument for type = "maxstat" </FIXME>
### lower = FALSE => p-value; lower = TRUE => 1 - p-value
doTest <- function(object, teststat = c("maximum", "quadratic", "scalar"), 
                   alternative = c("two.sided", "less", "greater"),
                   pvalue = TRUE, lower = FALSE, log = FALSE,
                   minbucket = 10L, ordered = TRUE, pargs = GenzBretz()) 
{

    ### avoid match.arg for performance reasons
    teststat <- teststat[1]
    stopifnot(teststat %in% c("maximum", "quadratic", "scalar"))
    alternative <- alternative[1]
    stopifnot(alternative %in% c("two.sided", "less", "greater"))

    if (teststat == "quadratic") stopifnot(alternative == "two.sided")

    test <- which(c("maximum", "quadratic", "scalar") == teststat)
    if (test == 3) {
        stopifnot(length(object$LinearStatistic) == 1)
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
                         as.double(pargs$abseps), as.double(pargs$releps), 
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
