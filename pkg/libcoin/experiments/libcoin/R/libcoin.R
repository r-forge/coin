
LinStatExpCov <- function(X, Y, weights, subset, block, varonly = FALSE) {

    stopifnot(NROW(X) == NROW(Y))

    if (!missing(weights)) {
        stopifnot(NROW(X) == length(weights))
        stopifnot(is.integer(weights))
    } else {
        weights <- integer(0)
    }
    if (!missing(subset)) {
        stopifnot(all(subset %in% 1:NROW(X)))
        stopifnot(is.integer(subset))
        subset <- subset - 1L
    } else {
        subset <- integer(0)
    }
    if (!missing(block)) {
        stopifnot(NROW(X) == length(block))
        stopifnot(is.factor(block))
    } else {
        block <- integer(0)
    }

    ms <- !(complete.cases(X) & complete.cases(Y))
    if (any(ms)) {
        OK <- FALSE
        if (length(subset) > 0) OK <- !any(which(ms) %in% subset)
        stopifnot(OK)
    }
    
    ret <- .Call("R_ExpectationCovarianceStatistic", X, Y, weights, subset, 
                 block, as.integer(varonly), PACKAGE = "libcoin")
    names(ret) <- c("LinStat", "Expectation", "Covariance")
    if (varonly)
        names(ret)[names(ret) == "Covariance"] <- "Variance"
    ret
}

LinStatExpCov2d <- function(X, Y, ix, iy, weights, subset, block, varonly = FALSE) {

    stopifnot(length(ix) == length(iy))
    stopifnot(is.integer(ix))
    stopifnot(is.integer(iy))
    attr(ix, "levels") <- 1:max(ix)
    attr(iy, "levels") <- 1:max(iy)

    if (missing(X)) X <- numeric(0)
    stopifnot(all(complete.cases(X)))
    stopifnot(all(complete.cases(Y)))
    stopifnot(nrow(X) == max(attr(ix, "levels")) + 1)
    stopifnot(nrow(Y) == max(attr(iy, "levels")) + 1)

    if (!missing(weights)) {
        stopifnot(length(ix) == length(weights))
        stopifnot(is.integer(weights))
    } else {
        weights <- integer(0)
    }
    if (!missing(subset)) {
        stopifnot(all(subset %in% 1:length(ix)))
        stopifnot(is.integer(subset))
        subset <- subset - 1L
    } else {
        subset <- integer(0)
    }
    if (!missing(block)) {
        stopifnot(length(ix) == length(block))
        stopifnot(is.factor(block))
    } else {
        block <- integer(0)
    }

    ret <- .Call("R_ExpectationCovarianceStatistic_2d", X, ix, Y, iy, 
        weights, subset, block, as.integer(varonly), PACKAGE = "libcoin")
    names(ret) <- c("LinStat", "Expectation", "Covariance")
    if (varonly)
        names(ret)[names(ret) == "Covariance"] <- "Variance"
    ret
}

ChisqStat <- function(object, tol = sqrt(.Machine$double.eps)) {

    stopifnot(!is.null(object$Covariance))
    tmp <- .Call("R_MPinv_sym", object$Covariance, tol)
    object$MPinv <- tmp[[1]]
    object$rank <- tmp[[2]]
    object
}

ChisqTest <- function(object, log = FALSE) {

    .Call("R_ChisqTest", object$LinStat, object$Expect, object$MPinv,
          object$rank, log = log)
}

MaxtypeStat <- function(object, tol = sqrt(.Machine$double.eps), 
                        alternative = c("two.sided", "less", "greater")) {

    alternative <- match.arg(alternative)
    object$alternative <- alternative
}
