
LinStatExpCov <- function(X, Y, weights, subset, block) {

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

    if (!missing(block)) {
       lev <- levels(block) 
       Cov <- 0
       Exp <- 0
       for (l in lev) {
           if (case == "vanilla")
               tmp <- LinStatExpCov(X, Y, subset = which(block == l))
           if (case == "weights")
               tmp <- LinStatExpCov(X, Y, weights = weights, subset = which(block == l))
           if (case == "subset")
               tmp <- LinStatExpCov(X, Y, subset = subset[which(block[subset] == l)])
           if (case == "weights_subset")
               tmp <- LinStatExpCov(X, Y, weights = weights, subset = subset[which(block[subset] == l)])
           Cov <- Cov + tmp$Covariance
           Exp <- Exp + tmp$Expectation
       }
       ret$Expectation <- Exp
       ret$Covariance <- Cov
    } else {
        ret$Expectation <- .Call("R_ExpectationLinearStatistic", ExpY, ExpX)
        ret$Covariance <- .Call("R_CovarianceLinearStatistic", CovY, ExpX, CovX, as.integer(sw))
    }
    ret
}
