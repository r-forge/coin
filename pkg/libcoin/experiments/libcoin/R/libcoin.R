
LinStatExpCov <- function(X, Y, weights, subset, block, varonly = FALSE) {

    stopifnot(NROW(X) == NROW(Y))
    if (!missing(weights)) {
        stopifnot(NROW(X) == length(weights))
        stopifnot(is.integer(weights))
    }
    if (!missing(subset)) {
        stopifnot(all(subset %in% 1:NROW(X)))
        stopifnot(is.integer(subset))
    }
    if (!missing(block)) {
        stopifnot(NROW(X) == length(block))
        stopifnot(is.factor(block))
    }

    ret <- list()
    if (missing(weights) & missing(subset)) case <- "vanilla"
    if (!missing(weights) & missing(subset)) case <- "weights"
    if (missing(weights) & !missing(subset)) case <- "subset"
    if (!missing(weights) & !missing(subset)) case <- "weights_subset"
    ret$case <- case

    cc <- complete.cases(X) & complete.cases(Y)
    if (case == "vanilla") {
        stopifnot(all(cc))
        sw <- NROW(X)
    }
    if (case == "weights") {
        stopifnot(all(cc[weights > 0]))
        sw <- sum(weights)
    }
    if (case == "subset") {
        stopifnot(all(cc[subset]))
        sw <- length(subset)
    }
    if (case == "weights_subset") {
        stopifnot(all(cc[(weights > 0) & (1:length(cc) %in% subset)]))
        sw <- sum(weights[subset])
    }

    ret$LinStat <- switch(case, 
        "vanilla" = .Call("R_LinearStatistic", X, Y, PACKAGE = "libcoin"),
        "weights" = .Call("R_LinearStatistic_weights", X, Y, weights, PACKAGE = "libcoin"),
        "subset" = .Call("R_LinearStatistic_subset", X, Y, subset - 1L, PACKAGE = "libcoin"),
        "weights_subset" = .Call("R_LinearStatistic_weights_subset", X, Y, weights, subset - 1L, PACKAGE = "libcoin"))
    ExpY <- switch(case, 
        "vanilla" = .Call("R_ExpectationInfluence", Y, PACKAGE = "libcoin"),
        "weights" = .Call("R_ExpectationInfluence_weights", Y, weights, PACKAGE = "libcoin"),
        "subset" = .Call("R_ExpectationInfluence_subset", Y, subset - 1L, PACKAGE = "libcoin"),
        "weights_subset" = .Call("R_ExpectationInfluence_weights_subset", Y, weights, subset - 1L, PACKAGE = "libcoin"))
    ExpX <- switch(case, 
        "vanilla" = .Call("R_ExpectationX", X, PACKAGE = "libcoin"),
        "weights" = .Call("R_ExpectationX_weights", X, weights, PACKAGE = "libcoin"),
        "subset" = .Call("R_ExpectationX_subset", X, subset - 1L, PACKAGE = "libcoin"),
        "weights_subset" = .Call("R_ExpectationX_weights_subset", X, weights, subset - 1L, PACKAGE = "libcoin"))
    if (varonly) {
        CovY <- switch(case, 
            "vanilla" = .Call("R_VarianceInfluence", Y, PACKAGE = "libcoin"),
            "weights" = .Call("R_VarianceInfluence_weights", Y, weights, PACKAGE = "libcoin"),
            "subset" = .Call("R_VarianceInfluence_subset", Y, subset - 1L, PACKAGE = "libcoin"),
            "weights_subset" = .Call("R_VarianceInfluence_weights_subset", Y, weights, subset - 1L, PACKAGE = "libcoin"))
        CovX <- switch(case, 
            "vanilla" = .Call("R_VarianceX", X, PACKAGE = "libcoin"),
            "weights" = .Call("R_VarianceX_weights", X, weights, PACKAGE = "libcoin"),
            "subset" = .Call("R_VarianceX_subset", X, subset - 1L, PACKAGE = "libcoin"),
            "weights_subset" = .Call("R_VarianceX_weights_subset", X, weights, subset - 1L, PACKAGE = "libcoin"))
    } else {
        CovY <- switch(case, 
            "vanilla" = .Call("R_CovarianceInfluence", Y, PACKAGE = "libcoin"),
            "weights" = .Call("R_CovarianceInfluence_weights", Y, weights, PACKAGE = "libcoin"),
            "subset" = .Call("R_CovarianceInfluence_subset", Y, subset - 1L, PACKAGE = "libcoin"),
            "weights_subset" = .Call("R_CovarianceInfluence_weights_subset", Y, weights, subset - 1L, PACKAGE = "libcoin"))
        CovX <- switch(case, 
            "vanilla" = .Call("R_CovarianceX", X, PACKAGE = "libcoin"),
            "weights" = .Call("R_CovarianceX_weights", X, weights, PACKAGE = "libcoin"),
            "subset" = .Call("R_CovarianceX_subset", X, subset - 1L, PACKAGE = "libcoin"),
            "weights_subset" = .Call("R_CovarianceX_weights_subset", X, weights, subset - 1L, PACKAGE = "libcoin"))
    }

    if (!missing(block)) {
       lev <- levels(block) 
       Cov <- 0
       Exp <- 0
       for (l in lev) {
           if (case == "vanilla")
               tmp <- LinStatExpCov(X, Y, subset = which(block == l), varonly = varonly)
           if (case == "weights")
               tmp <- LinStatExpCov(X, Y, weights = weights, subset = which(block == l), varonly = varonly)
           if (case == "subset")
               tmp <- LinStatExpCov(X, Y, subset = subset[which(block[subset] == l)], varonly = varonly)
           if (case == "weights_subset")
               tmp <- LinStatExpCov(X, Y, weights = weights, subset = subset[which(block[subset] == l)], varonly = varonly)
           if (varonly) {
               Cov <- Cov + tmp$Variance
           } else {
               Cov <- Cov + tmp$Covariance
           }
           Exp <- Exp + tmp$Expectation
       }
       ret$Expectation <- Exp
       ret$Covariance <- Cov
    } else {
        ret$Expectation <- .Call("R_ExpectationLinearStatistic", ExpY, ExpX)
        if (!varonly) {
            ret$Covariance <- .Call("R_CovarianceLinearStatistic", CovY, ExpX, CovX, as.integer(sw))
        } else {
            ret$Covariance <- .Call("R_VarianceLinearStatistic", CovY, ExpX, CovX, as.integer(sw))
        }
    }
    if (varonly) names(ret) <- gsub("Covariance", "Variance", names(ret))
    ret
}


LinStatExpCov2d <- function(X, Y, ix, iy, weights, subset, block, varonly = FALSE) {

    stopifnot(length(ix) == length(iy))

    stopifnot(is.integer(ix))
    stopifnot(is.integer(iy))

    stopifnot(all(complete.cases(X)))
    stopifnot(all(complete.cases(Y)))

    attr(ix, "levels") <- 1:max(ix)
    attr(iy, "levels") <- 1:max(iy)

    if (!missing(weights)) {
        stopifnot(length(ix) == length(weights))
        stopifnot(is.integer(weights))
    }
    if (!missing(subset)) {
        stopifnot(all(subset <= length(ix)))
        stopifnot(is.integer(subset))
    }
    if (!missing(block)) {
        stopifnot(length(ix) == length(block))
        stopifnot(is.factor(block))
    }

    ret <- list()
    if (missing(weights) & missing(subset)) case <- "vanilla"
    if (!missing(weights) & missing(subset)) case <- "weights"
    if (missing(weights) & !missing(subset)) case <- "subset"
    if (!missing(weights) & !missing(subset)) case <- "weights_subset"

    ret$case <- case
    if (missing(block)) {
        tab_ixiy <- switch(case, "vanilla" = .Call("R_2dtable", ix, iy, PACKAGE = "libcoin"),
                "weights" = .Call("R_2dtable_weights", ix, iy, weights, PACKAGE = "libcoin"),
                "subset" = .Call("R_2dtable_subset", ix, iy, subset - 1L, PACKAGE = "libcoin"),
                "weights_subset" = .Call("R_2dtable_weights_subset", ix, iy, weights, subset - 1L, PACKAGE = "libcoin"))
    } else {
        tab_ixiy <- switch(case, "vanilla" = .Call("R_2dtable_block", ix, iy, block, PACKAGE = "libcoin"),
                "weights" = .Call("R_2dtable_weights_block", ix, iy, weights, block, PACKAGE = "libcoin"),
                "subset" = .Call("R_2dtable_subset_block", ix, iy, subset - 1L, block, PACKAGE = "libcoin"),
                "weights_subset" = .Call("R_2dtable_weights_subset_block", ix, iy, weights, subset - 1L, block, PACKAGE = "libcoin"))
        tab_ixiy <- array(tab_ixiy, dim = c(max(ix) + 1, max(iy) + 1, nlevels(block)))
    }

    fun <- function(tab) {
        ret <- list()
        sw <- sum(tab)
        ret$LinStat <- .Call("R_LinearStatistic_2d", X, Y, tab, PACKAGE = "libcoin")
        tab_iy <- as.integer(colSums(tab))
        tab_ix <- as.integer(rowSums(tab))
        ExpY <- .Call("R_ExpectationInfluence_weights", Y, tab_iy, PACKAGE = "libcoin")
        ExpX <- .Call("R_ExpectationX_weights", X, tab_ix, PACKAGE = "libcoin")
        if (varonly) {
            CovY <- .Call("R_VarianceInfluence_weights", Y, tab_iy, PACKAGE = "libcoin")
            CovX <- .Call("R_VarianceX_weights", X, tab_ix, PACKAGE = "libcoin")
        } else {
            CovY <- .Call("R_CovarianceInfluence_weights", Y, tab_iy, PACKAGE = "libcoin")
            CovX <- .Call("R_CovarianceX_weights", X, tab_ix, PACKAGE = "libcoin")
        }
        ret$Expectation <- .Call("R_ExpectationLinearStatistic", ExpY, ExpX)
        if (!varonly) {
            ret$Covariance <- .Call("R_CovarianceLinearStatistic", CovY, ExpX, CovX, as.integer(sw))
        } else {
            ret$Covariance <- .Call("R_VarianceLinearStatistic", CovY, ExpX, CovX, as.integer(sw))
        }
        ret
    }

    if (missing(block)) {
        ret <- fun(tab_ixiy)
    } else {
        tmp <- apply(tab_ixiy, 3, fun)
        ret$LinStat <- do.call("+", lapply(tmp, function(o) o$LinStat))
        ret$Expectation <- do.call("+", lapply(tmp, function(o) o$Expectation))
        ret$Covariance <- do.call("+", lapply(tmp, function(o) o$Covariance))
    }
    if (varonly) names(ret) <- gsub("Covariance", "Variance", names(ret))
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
