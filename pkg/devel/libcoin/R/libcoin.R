
# R Header

###
### TO NOT EDIT THIS FILE
### 
### Edit `libcoin.w' and run `nuweb -r libcoin.w'
###

# LinStatExpCov

LinStatExpCov <- function(X, Y, ix = NULL, iy = NULL, weights = integer(0),
                          subset = integer(0), block = integer(0),
                          checkNAs = TRUE, 
                          varonly = FALSE, nperm = B, B = 0, standardise = FALSE,
                          tol = sqrt(.Machine$double.eps))
{
    if (missing(X) & !is.null(ix) & is.null(iy)) {
        X <- ix
        ix <- NULL
    }

    if (missing(X)) X <- integer(0)

    ### <FIXME> for the time being only!!! </FIXME>
##    if (length(subset) > 0) subset <- sort(subset)
    
    if (is.null(ix) & is.null(iy))
        return(.LinStatExpCov1d(X = X, Y = Y, weights = weights,
                                subset = subset, block = block, 
                                checkNAs = checkNAs,
                                varonly = varonly, nperm = nperm,
                                standardise = standardise, tol = tol))

    if (!is.null(ix) & !is.null(iy))
        return(.LinStatExpCov2d(X = X, Y = Y, ix = ix, iy = iy,
                                weights = weights, subset = subset,
                                block = block, varonly = varonly, 
                                checkNAs = checkNAs, nperm = nperm,
                                standardise = standardise, tol = tol))

    stop("incorrect call to LinStatExpCov")
}

# LinStatExpCov1d

.LinStatExpCov1d <- function(X, Y, weights = integer(0), subset = integer(0), block = integer(0),
                             checkNAs = TRUE, varonly = FALSE, nperm = 0, standardise = FALSE,
                             tol = sqrt(.Machine$double.eps))
{

    if (NROW(X) != NROW(Y))
        stop("dimensions of X and Y don't match")
    N <- NROW(X)

    if (is.integer(X)) {
        if (is.null(attr(X, "levels")) || checkNAs) {
            rg <- range(X)
            if (any(is.na(rg)))
                stop("no missing values allowed in X") 
            stopifnot(rg[1] > 0) ### no missing values allowed here!
            if (is.null(attr(X, "levels")))
                attr(X, "levels") <- 1:rg[2]
        }
    }

    if (is.factor(X) && checkNAs)
        stopifnot(all(!is.na(X)))

    # Check weights, subset, block
    

    if (is.null(weights)) weights <- integer(0)

    if (length(weights) > 0) {
        if (!((N == length(weights)) && all(weights >= 0)))
            stop("incorrect weights")
        if (checkNAs) stopifnot(all(!is.na(weights)))
    }

    if (is.null(subset)) subset <- integer(0)

    if (length(subset) > 0 && checkNAs) {
        rs <- range(subset)
        if (any(is.na(rs))) stop("no missing values allowed in subset")
        if (!((rs[2] <= N) && (rs[1] >= 1L)))
            stop("incorrect subset")
    }

    if (is.null(block)) block <- integer(0)

    if (length(block) > 0) {
        if (!((N == length(block)) && is.factor(block)))
            stop("incorrect block")
        if (checkNAs) stopifnot(all(!is.na(block)))
    }
    

    if (checkNAs) {
        # Handle Missing Values
        
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
        
    }

    ret <- .Call(R_ExpectationCovarianceStatistic, X, Y, weights, subset,
                 block, as.integer(varonly), as.double(tol))
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    if (nperm > 0) {
        ret$PermutedLinearStatistic <-
            .Call(R_PermutedLinearStatistic, X, Y, weights, subset,
                  block, as.double(nperm))
        if (standardise)
            ret$StandardisedPermutedLinearStatistic <-
                .Call(R_StandardisePermutedLinearStatistic, ret)
    }
    class(ret) <- c("LinStatExpCov1d", "LinStatExpCov")
    ret
}

# LinStatExpCov2d

.LinStatExpCov2d <- function(X = numeric(0), Y, ix, iy, weights = integer(0), subset = integer(0),
                             block = integer(0), checkNAs = TRUE, varonly = FALSE, nperm = 0,
                             standardise = FALSE,
                             tol = sqrt(.Machine$double.eps))
{

    IF <- function(x) is.integer(x) || is.factor(x)

    if (!((length(ix) == length(iy)) && IF(ix) && IF(iy)))
        stop("incorrect ix and/or iy")
    N <- length(ix)

    # Check ix
    
    if (is.null(attr(ix, "levels"))) {
        rg <- range(ix)
        if (any(is.na(rg)))
            stop("no missing values allowed in ix") 
        stopifnot(rg[1] >= 0)
        attr(ix, "levels") <- 1:rg[2]
    } else {
        if (checkNAs) stopifnot(all(!is.na(ix)))
    }
    

    # Check iy
    
    if (is.null(attr(iy, "levels"))) {
        rg <- range(iy)
        if (any(is.na(rg)))
            stop("no missing values allowed in iy") 
        stopifnot(rg[1] >= 0)
        attr(iy, "levels") <- 1:rg[2]
    } else {
        if (checkNAs) stopifnot(all(!is.na(ix)))
    }
    

    if (length(X) > 0) {
        if (!(NROW(X) == (length(attr(ix, "levels")) + 1) &&
              all(complete.cases(X)) &&
             (NROW(X) == (length(attr(ix, "levels")) + 1))))
            stop("incorrect X")
    }

    if (!(all(complete.cases(Y))) &&
          (NROW(Y) == (length(attr(iy, "levels")) + 1)) &&
          (NROW(Y) == (length(attr(iy, "levels")) + 1)))
        stop("incorrect Y")

    # Check weights, subset, block
    

    if (is.null(weights)) weights <- integer(0)

    if (length(weights) > 0) {
        if (!((N == length(weights)) && all(weights >= 0)))
            stop("incorrect weights")
        if (checkNAs) stopifnot(all(!is.na(weights)))
    }

    if (is.null(subset)) subset <- integer(0)

    if (length(subset) > 0 && checkNAs) {
        rs <- range(subset)
        if (any(is.na(rs))) stop("no missing values allowed in subset")
        if (!((rs[2] <= N) && (rs[1] >= 1L)))
            stop("incorrect subset")
    }

    if (is.null(block)) block <- integer(0)

    if (length(block) > 0) {
        if (!((N == length(block)) && is.factor(block)))
            stop("incorrect block")
        if (checkNAs) stopifnot(all(!is.na(block)))
    }
    

    ret <- .Call(R_ExpectationCovarianceStatistic_2d, X, ix, Y, iy,
                 weights, subset, block, as.integer(varonly), as.double(tol))
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    if (nperm > 0) {
        ret$PermutedLinearStatistic <-
            .Call(R_PermutedLinearStatistic_2d, X, ix, Y, iy, block, nperm, ret$Table)
        if (standardise)
            ret$StandardisedPermutedLinearStatistic <-
                .Call(R_StandardisePermutedLinearStatistic, ret)
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
                     pvalue = TRUE, lower = FALSE, log = FALSE, PermutedStatistics = FALSE,
                     minbucket = 10L, ordered = TRUE, maxselect = object$Xfactor, 
                     pargs = GenzBretz())
{

    ### avoid match.arg for performance reasons
    teststat <- teststat[1]
    if (!any(teststat == c("maximum", "quadratic", "scalar")))
        stop("incorrect teststat")
    alternative <- alternative[1]
    if (!any(alternative == c("two.sided", "less", "greater")))
        stop("incorrect alternative")

    if (maxselect)
        stopifnot(object$Xfactor)

    if (teststat == "quadratic" || maxselect) {
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

    if (!maxselect) {
        if (teststat == "quadratic") {
            ret <- .Call(R_QuadraticTest, object,
                         as.integer(pvalue), as.integer(lower),
                         as.integer(log), as.integer(PermutedStatistics))
        } else {
            ret <- .Call(R_MaximumTest, object,
                         as.integer(alt), as.integer(pvalue), as.integer(lower),
                         as.integer(log), as.integer(PermutedStatistics),
                         as.integer(pargs$maxpts),
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
        if (!PermutedStatistics) ret$PermutedStatistics <- NULL
    }
    ret
}

# Contrasts

`%*%` <- function(x, y) UseMethod("%*%", y)
`%*%.default` <- function(x, y) base::`%*%`(x, y)
`%*%.LinStatExpCov` <- function(x, y) {
    stopifnot(!y$varonly)
    stopifnot(is.numeric(x))
    if (is.vector(x)) x <- matrix(x, nrow = 1)
    P <- y$dimension[1]
    stopifnot(ncol(x) == P)
    Q <- y$dimension[2]
    ret <- y
    xLS <- x %*% matrix(y$LinearStatistic, nrow = P)
    xExp <- x %*% matrix(y$Expectation, nrow = P)
    xExpX <- x %*% matrix(y$ExpectationX, nrow = P)
    xCov <- tcrossprod(x %*% vcov(y), x)
    if (!is.matrix(xCov)) xCov <- matrix(xCov)
    if (length(y$PermutedLinearStatistic) > 0) {
        xPS <- apply(y$PermutedLinearStatistic, 2, function(y)
                     as.vector(x %*% matrix(y, nrow = P)))
        if (!is.matrix(xPS)) xPS <- matrix(xPS, nrow = 1)
        ret$PermutedLinearStatistic <- xPS
    }
    ret$LinearStatistic <- as.vector(xLS)
    ret$Expectation <- as.vector(xExp)
    ret$ExpectationX <- as.vector(xExpX)
    ret$Covariance <- as.vector(xCov[lower.tri(xCov, diag = TRUE)])
    ret$Variance <- diag(xCov)
    ret$dimension <- c(NROW(x), Q)
    ret$Xfactor <- FALSE
    if (length(y$StandardisedPermutedLinearStatistic) > 0)
        ret$StandardisedPermutedLinearStatistic <-
            .Call(R_StandardisePermutedLinearStatistic, ret)
    ret
}

.libcoinCall <- function(FUN, ...) .Call(FUN, ...)
