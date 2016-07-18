
LinStatExpCov <- function(X, Y, weights, subset, block, 
                          varonly = FALSE, B = 0L, standardise = FALSE, 
                          tol = sqrt(.Machine$double.eps)) 
{
    stopifnot(NROW(X) == NROW(Y))

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
    ret$sim <- double(0);
    if (B > 0)
        ret$sim <- .Call("R_PermutedLinearStatistic", ret, X, Y, weights, subset, 
                         block, as.integer(B), as.integer(standardise), as.double(tol),
                         PACKAGE = "libcoin")
    ret
}

LinStatExpCov2d <- function(X, Y, ix, iy, weights, subset, block, varonly = FALSE, B = 0,
                            standardise = FALSE, tol = sqrt(.Machine$double.eps)) 
{
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
    ret$sim <- double(0);
    if (B > 0)
        ret$sim <- .Call("R_PermutedLinearStatistic_2d", ret, X, ix, Y, iy, 
                         block, as.integer(B), as.integer(standardise), as.double(tol),
                         PACKAGE = "libcoin")
    ret
}

### <FIXME> add alternative argument for type = "maxstat" </FIXME>
Test <- function(object, tol = sqrt(.Machine$double.eps), lower = FALSE, log = FALSE,
                 type = c("maxstat", "quadform"), xtrafo = c("id", "maxstat"),
                 minbucket = 10L, ordered = TRUE) 
{
    type <- match.arg(type)
    xtrafo <- match.arg(xtrafo)
    if (xtrafo == "id") {
        if (type == "quadform") {
            ret <- .Call("R_ChisqTest", object, object$sim, tol, 
                         as.integer(lower), as.integer(log), PACKAGE = "libcoin")
        } else {
            ret <- .Call("R_MaxabsstatTest", object, object$sim, tol, as.integer(lower), 
                         as.integer(log), 10000L, .0001, .0001, PACKAGE = "libcoin")
        }
    } else {
        type <- as.integer(which(c("maxstat", "quadform") == type))
        if (ordered) {
            ret <- .Call("R_MaxstatTest_ordered", object, object$sim, type, tol, 
                        as.integer(minbucket), as.integer(lower), as.integer(log), PACKAGE = "libcoin")
        } else {
            ret <- .Call("R_MaxstatTest_unordered", object, object$sim, type, tol, 
                     as.integer(minbucket), as.integer(lower), as.integer(log), PACKAGE = "libcoin")
        }
    }
    if (length(ret) == 2)
        names(ret) <- c("statistic", "p.value")
    if (length(ret) == 3)
        names(ret) <- c("statistic", "p.value", "index")
    ret
}
