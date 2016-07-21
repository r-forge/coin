
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
            subset <- subset[!(subset %in% which(ms))] 
        } else {
            subset <- (0:(NROW(X) - 1))[-which(ms)]
        }
    }
    
    ret <- .Call("R_ExpectationCovarianceStatistic", X, Y, weights, subset, 
                 block, as.integer(varonly), PACKAGE = "libcoin")
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    ret$sim <- double(0)
    ret$tol <- tol
    if (B > 0)
        ret$sim <- .Call("R_PermutedLinearStatistic", ret, X, Y, weights, 
                         subset, block, as.integer(B), as.integer(standardise), 
                         as.double(tol), PACKAGE = "libcoin")
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

    if (missing(X)) X <- numeric(0)
    stopifnot(all(complete.cases(X)))
    stopifnot(all(complete.cases(Y)))
    stopifnot(nrow(X) == length(attr(ix, "levels")) + 1)
    stopifnot(nrow(Y) == length(attr(iy, "levels")) + 1)

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
        weights, subset, block, as.integer(varonly), PACKAGE = "libcoin")
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    ret$sim <- double(0)
    ret$tol <- tol
    if (B > 0)
        ret$sim <- .Call("R_PermutedLinearStatistic_2d", ret, X, ix, Y, iy, 
                         block, as.integer(B), as.integer(standardise), 
                         as.double(tol), PACKAGE = "libcoin")
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
doTest <- function(object, type = c("maxstat", "quadform"), 
                   alternative = c("two.sided", "less", "greater"),
                   lower = FALSE, log = FALSE,
                   minbucket = 10L, ordered = TRUE, pargs = GenzBretz()) 
{

    type <- match.arg(type)
    alternative <- match.arg(alternative)
    if (type == "quadform") stopifnot(alternative == "two.sided")
    alt <- which(c("two.sided", "less", "greater") == alternative)
    if (!object$Xfactor) {
        if (type == "quadform") {
            ret <- .Call("R_ChisqTest", object, object$sim, object$tol, 
                         as.integer(lower), as.integer(log), 
                         PACKAGE = "libcoin")
        } else {
            ret <- .Call("R_MaxtypeTest", object, object$sim, object$tol, 
                         as.integer(alt), as.integer(lower), 
                         as.integer(log), as.integer(pargs$maxpts), 
                         as.double(pargs$abseps), as.double(pargs$releps), 
                         PACKAGE = "libcoin")
        }
    } else {
        type <- as.integer(which(c("maxstat", "quadform") == type))
        ret <- .Call("R_MaxSelectTest", object, as.integer(ordered), 
                     object$sim, type, object$tol, as.integer(minbucket), 
                     as.integer(lower), as.integer(log), PACKAGE = "libcoin")
    }
    ret
}
