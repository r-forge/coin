
LinStatExpCov <- function(X, Y, weights, subset, block, 
                          varonly = FALSE, B = 0L, standardise = FALSE, 
                          tol = sqrt(.Machine$double.eps)) {

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
    ret$sim <- double(0);
    if (B > 0)
        ret$sim <- .Call("R_PermutedLinearStatistic", ret, X, Y, weights, subset, 
                         block, as.integer(B), as.integer(standardise), as.double(tol),
                         PACKAGE = "libcoin")
    ret
}

LinStatExpCov2d <- function(X, Y, ix, iy, weights, subset, block, varonly = FALSE, B = 0,
                            standardise = FALSE, tol = sqrt(.Machine$double.eps)) {

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
    ret$sim <- double(0);
    if (B > 0)
        ret$sim <- .Call("R_PermutedLinearStatistic_2d", ret, X, ix, Y, iy, 
                         block, as.integer(B), as.integer(standardise), as.double(tol),
                         PACKAGE = "libcoin")
    ret
}

ChisqStat <- function(object, tol = sqrt(.Machine$double.eps)) {

    stopifnot(!is.null(object$Covariance))
    tmp <- .Call("R_MPinv_sym", object$Covariance, tol)
    object$MPinv <- tmp[[1]]
    object$rank <- tmp[[2]]
    object
}

ChisqTest <- function(object, tol = sqrt(.Machine$double.eps), log = FALSE) {
    .Call("R_ChisqTest", object, object$sim, sqrt(.Machine$double.eps), as.integer(log))
}

MaxabsstatTest <- function(object, tol = sqrt(.Machine$double.eps), log = FALSE) {
    .Call("R_MaxabsstatTest", object, object$sim, sqrt(.Machine$double.eps), as.integer(log), 
          10000L, .0001, .0001)
}

MaxstatTest <- function(object, tol = sqrt(.Machine$double.eps), minbucket = 10L, 
                        log = FALSE, ordered = TRUE, type = c("maxabs", "quadform")) {
    type <- match.arg(type)
    type <- as.integer(which(c("maxabs", "quadform") == type))
    if (ordered) 
        return(.Call("R_MaxstatTest_ordered", object, object$sim, type, tol, 
                     as.integer(minbucket), as.integer(log), PACKAGE = "libcoin"))
    return(.Call("R_MaxstatTest_unordered", object, object$sim, type, tol, 
                 as.integer(minbucket), as.integer(log), PACKAGE = "libcoin"))
}


