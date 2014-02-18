### Streitberg-Roehmel algorithm for independent two samples
SR_shift_2sample <- function(object, fact) {

    if (!extends(class(object), "ScalarIndependenceTestStatistic"))
        stop("Argument ", sQuote("object"), " is not of class ",
             sQuote("ScalarIndependenceTestStatistic"))

    if (!is_2sample(object))
        stop(sQuote("object"),
             " does not represent an independent two-sample problem")

    ytrans <- object@ytrans[, 1L]
    xtrans <- object@xtrans[, 1L]
    block <- object@block

    ## expand observations if weights are non-unity
    if (!is_unity(object@weights)) {
        idx <- rep.int(seq_along(object@weights), object@weights)
        ytrans <- ytrans[idx]
        xtrans <- xtrans[idx]
        block <- block[idx]
    }

    T <- 0
    Prob <- 1
    for (lev in levels(block)) {

        thisblock <- (block == lev)

        ## compute distribution of scores in this block
        scores <- ytrans[thisblock]
        m <- sum(xtrans[thisblock] == 1L)

        if (m == 0L) next;
        if (m == length(scores))
            dens <- list(T = sum(scores), Prob = 1)
        if (m < length(scores))
            dens <- cSR_shift_2sample(scores, m, fact = fact)

        ## update distribution of statistic over all blocks
        T <- as.vector(outer(dens$T, T, "+"))
        Prob <- drop(kronecker(Prob, dens$Prob))
    }

    T <- (T - expectation(object)) / sqrt(variance(object))

    ## T may not be distinct and ordered if blocks are present
    if (nlevels(block) > 1L) {
        n <- length(T)
        o <- order(T)
        T <- T[o]
        idx <- c(which(T[-1L] - T[-n] > eps()), n)
        T <- T[idx]
        Prob <-
            vapply(split(Prob[o], rep.int(seq_along(idx), diff(c(0L, idx)))),
                   sum, NA_real_, USE.NAMES = FALSE)
    }

    new("ExactNullDistribution",
        p = function(q) sum(Prob[LE(T, q)]),
        q = function(p) {
            idx <- which(cumsum(Prob) < p)
            if (length(idx) == 0L)
                T[1L]
            else if (length(idx) == length(Prob))
                T[max(idx)]
            else
                T[max(idx) + 1L]
        },
        d = function(x) Prob[T == x],
        pvalue = function(q) {
            switch(object@alternative,
                "less"      = sum(Prob[LE(T, q)]),
                "greater"   = sum(Prob[GE(T, q)]),
                "two.sided" = {
                    if (q == 0)
                        1
                    else
                        sum(Prob[LE(T, ifelse(q > 0, -q, q))]) +
                          sum(Prob[GE(T, ifelse(q >= 0, q, -q))])
                }
            )
        },
        support = function() T,
        name = paste0("exact independent two-sample distribution",
                      " (via Streitberg-Roehmel algorithm)"))
}

cSR_shift_2sample <- function(scores, m, fact) {

    if (m < 1L || m == length(scores))
        stop("not a two sample problem")
    n <- length(scores)
    ones <- rep.int(1L, n)

    ## integer scores with sum(scores) minimal
    scores <- scores * fact
    add <- min(scores - 1)
    scores <- scores - add
    m_b <- sum(sort(scores)[(n + 1L - m):n])

    Prob <- .Call("R_cpermdist2",
                  score_a = as.integer(ones), score_b = as.integer(scores),
                  m_a = as.integer(m), m_b = as.integer(m_b),
                  retProb = as.logical(TRUE),
                  PACKAGE = "coin")
    T <- which(Prob != 0)

    list(T = (T + add * m) / fact, Prob = Prob[T])
}


### Streitberg-Roehmel algorithm for paired samples
SR_shift_1sample <- function(object, fact) {

    if (!extends(class(object), "ScalarIndependenceTestStatistic"))
        stop("Argument ", sQuote("object"), " is not of class ",
             sQuote("ScalarIndependenceTestStatistic"))

    if (!is_2sample(object))
        stop(sQuote("object"),
             " does not represent an independent two-sample problem")

    scores <- object@ytrans[, 1L]
    if (any(scores < 0))
        stop("cannot compute exact distribution with negative scores")
    block <- object@block

    ## expand observations if weights are non-unity
    if (!is_unity(object@weights)) {
        idx <- rep.int(seq_along(object@weights), object@weights)
        scores <- scores[idx]
        block <- block[idx]
    }

    ##  table(object@block, scores == 0) checken
    scores <- vapply(unique(object@block), function(i) {
        s <- round(scores * fact)[object@block == i]
        s[s != 0] # remove zeros
    }, NA_real_)
    storage.mode(scores) <- "integer"
    Prob <- .Call("R_cpermdist1", scores, PACKAGE = "coin")
    T <- which(Prob != 0)
    Prob <- Prob[T]
    ## 0 is possible
    T <- (T - 1) / fact

    T <- (T - expectation(object)) / sqrt(variance(object))

    new("ExactNullDistribution",
        p = function(q) sum(Prob[LE(T, q)]),
        q = function(p) {
            idx <- which(cumsum(Prob) < p)
            if (length(idx) == 0L)
                T[1L]
            else if (length(idx) == length(Prob))
                T[max(idx)]
            else
                T[max(idx) + 1L]
        },
        d = function(x) Prob[T == x],
        pvalue = function(q) {
            switch(object@alternative,
                "less"      = sum(Prob[LE(T, q)]),
                "greater"   = sum(Prob[GE(T, q)]),
                "two.sided" = {
                    if (q == 0)
                        1
                    else
                        sum(Prob[LE(T, ifelse(q > 0, -q, q))]) +
                          sum(Prob[GE(T, ifelse(q >= 0, q, -q))])
                }
            )
        },
        support = function() T,
        name = paste0("exact paired two-sample distribution",
                      " (via Streitberg-Roehmel algorithm)"))
}


### van de Wiel split-up algorithm for independent two samples
vdW_split_up_2sample <- function(object) {

    ## <FIXME> on.exit(ex <- .C("FreeW", PACKAGE = "coin")) </FIXME>

    if (!extends(class(object), "ScalarIndependenceTestStatistic"))
        stop("Argument ", sQuote("object"), " is not of class ",
             sQuote("ScalarIndependenceTestStatistic"))

    if (!is_2sample(object))
        stop(sQuote("object"),
             " does not represent an independent two-sample problem")

    if (nlevels(object@block) != 1L)
        stop("cannot compute exact p-values with blocks")

    scores <- object@ytrans[, 1L]
    xtrans <- object@xtrans[, 1L]

    ## expand observations if weights are non-unity
    if (!is_unity(object@weights)) {
        idx <- rep.int(seq_along(object@weights), object@weights)
        scores <- scores[idx]
        xtrans <- xtrans[idx]
    }

    storage.mode(scores) <- "double"
    m <- sum(xtrans)
    storage.mode(m) <- "integer"
    tol <- eps()

    CumProb <- function(q) {
        T <- q * sqrt(variance(object)) + expectation(object)
        .Call("R_split_up_2sample", scores, m, T, tol, PACKAGE = "coin")
    }

    new("ExactNullDistribution",
        p = CumProb,
        q = function(p) {
            f <- function(x) CumProb(x) - p
            rr <- if (p <= 0.5)
                      uniroot(f, interval = c(-10, 1), tol = tol)
                  else
                      uniroot(f, interval = c(-1, 10), tol = tol)
            ## make sure quantile leads to pdf >= p
            if (rr$f.root < 0)
                rr$root <- rr$root + tol
            ## pdf is constant here
            if (rr$estim.prec > tol) {
                r1 <- rr$root
                d <- min(diff(sort(scores[!duplicated(scores)]))) /
                       sqrt(variance(object))
                while (d > tol) {
                    if (f(r1 - d) >= 0)
                        r1 <- r1 - d
                    else
                        d <- d / 2
                }
                rr$root <- r1
            }
            rr$root
        },
        d = function(x) NA,
        pvalue = function(q) {
            switch(object@alternative,
                "less"      = CumProb(q),
                "greater"   = 1 - CumProb(q - 10 * tol),
                "two.sided" = {
                    if (q == 0)
                        1
                    else if (q > 0)
                        CumProb(-q) + (1 - CumProb(q - 10 * tol))
                    else
                        CumProb(q) + (1 - CumProb(-q - 10 * tol))
                }
            )
        },
        support = function(p = 1e-5) NA,
        name = paste0("exact independent two-sample distribution",
                      " (via van de Wiel split-up algorithm)"))
}
