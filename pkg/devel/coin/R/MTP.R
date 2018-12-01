### joint distribution-based max-T multiple testing procedures
### (Westfall and Young, 1993)
setGeneric("joint",
    function(object1, object2, ...) {
        standardGeneric("joint")
    }
)

setMethod("joint",
    signature = list("MaxTypeIndependenceTestStatistic", "NullDistribution"),
    definition = function(object1, object2, stepdown, ...) {
        ## reorder test statistics to ensure consistency with "global"/"step-down"
        switch(object1@alternative,
            "less" = {
                z <- statistic(object1, type = "standardized")
                o <- order(z)                    # smallest z first
                pq <- length(z)
                if (stepdown) {
                    upper <- rep.int(Inf, pq)
                    lower <- z[o]
                }
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                o <- order(z, decreasing = TRUE) # largest z first
                pq <- length(z)
                if (stepdown) {
                    upper <- z[o]
                    lower <- rep.int(-Inf, pq)
                }
            },
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                o <- order(z, decreasing = TRUE) # abs. largest z first
                pq <- length(z)
                if (stepdown) {
                    upper <- z[o]
                    lower <- -upper
                }
            }
        )

        if (!stepdown) {
            ## iterate over unique test statistics only and remap
            RET <- z[o]
            idx <- c(which(RET[-1L] %NE% RET[-pq]), pq) # unique z
            RET <- pvalue(object2, RET[idx], ...)
            RET <- rep.int(RET, diff(c(0L, idx)))  # remapping
        } else {
            ## free step-down based on multivariate normality
            Rho <- cov2cor(covariance(object1))
            RET <- numeric(pq)
            RET[1] <- pmvn(lower = lower[1], upper = upper[1],
                           mean = rep.int(0, pq), corr = Rho,
                           conf.int = FALSE)
            if (pq > 1) {
                oo <- o
                for (i in 2:pq) {
                    j <- rank(oo)[1] # reindexing needed in each step
                    Rho <- Rho[-j, -j]
                    oo <- oo[-1]
                    RET[i] <- min(RET[i - 1],
                                  pmvn(lower = lower[i], upper = upper[i],
                                       mean = rep.int(0, length(oo)), corr = Rho,
                                       conf.int = FALSE))
                }
            }
            RET <- 1 - RET
        }

        RET <- matrix(RET[order(o)], nrow = nrow(z), ncol = ncol(z),
                      dimnames = dimnames(z))
        class(RET) <- "pvalue"
        RET
    }
)

setMethod("joint",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, stepdown, ...) {
        if (!stepdown) {
            RET <- callNextMethod(object1, object2, stepdown, ...)
        } else {
            ## free step-down based on the resampling distribution
            ## (Westfall & Young, 1993, p. 66-67, Algorithm 2.8)
            ## using standardized statistics instead of p-values
            mu <- expectation(object1)
            sigma <- sqrt(variance(object1))
            switch(object1@alternative,
                "less" = {
                    z <- statistic(object1, type = "standardized")
                    o <- order(z, decreasing = TRUE) # largest z first
                    RET <- (support(object2, raw = TRUE) - mu) / sigma
                    RET <- rowMeans(colCummins(RET[o, ]) %LE% z[o])
                },
                "greater" = {
                    z <- statistic(object1, type = "standardized")
                    o <- order(z)                    # smallest z first
                    RET <- (support(object2, raw = TRUE) - mu) / sigma
                    RET <- rowMeans(colCummaxs(RET[o, ]) %GE% z[o])
                },
                "two.sided" = {
                    z <- abs(statistic(object1, type = "standardized"))
                    o <- order(z)                    # abs. smallest z first
                    RET <- abs(support(object2, raw = TRUE) - mu) / sigma
                    RET <- rowMeans(colCummaxs(RET[o, ]) %GE% z[o])
                }
            )
            RET <- rev(cummax(rev(RET))) # enforce monotonicity

            RET <- matrix(RET[order(o)], nrow = nrow(z), ncol = ncol(z),
                          dimnames = dimnames(z))
            class(RET) <- "pvalue"
        }

        attr(RET, "nresample") <- object2@nresample
        RET
    }
)

setMethod("joint",
    signature = list("MaxTypeIndependenceTest", "missing"),
    definition = function(object1, object2, stepdown, ...) {
        callGeneric(object1@statistic, object1@distribution, stepdown, ...)
    }
)


### marginal distribution-based max-T multiple testing procedures
### (Westfall and Wolfinger, 1997; Westfall and Troendle, 2008)
setGeneric("marginal",
    function(object1, object2, ...) {
        standardGeneric("marginal")
    }
)

setMethod("marginal",
    signature = list("MaxTypeIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, stepdown, bonferroni, ...) {
        ## unadjusted p-values
        z <- statistic(object1, type = "standardized")
        RET <- switch(object1@alternative,
                   "less"      = pnorm(z),
                   "greater"   = 1 - pnorm(z),
                   "two.sided" = 2 * pmin.int(pnorm(z), 1 - pnorm(z))
               )

        ## adjustment
        RET <- if (!stepdown) {
                   if (bonferroni) pmin.int(1, length(RET) * RET) # Bonferroni
                   else 1 - (1 - RET)^length(RET)                 # Sidak
               } else {
                   n <- length(RET)
                   o <- order(RET)
                   if (bonferroni) # Bonferroni-Holm
                       pmin.int(1, cummax((n - seq_len(n) + 1L) * RET[o])[order(o)])
                   else            # Sidak-Holm
                       cummax(1 - (1 - RET[o])^(n - seq_len(n) + 1L))[order(o)]
               }

        RET <- matrix(RET, nrow = nrow(z), ncol = ncol(z),
                      dimnames = dimnames(z))
        class(RET) <- "pvalue"
        RET
    }
)

setMethod("marginal",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, stepdown, bonferroni, ...) {
        ## standardized observed and permuted test statistics
        mu <- expectation(object1)
        sigma <- sqrt(variance(object1))
        switch(object1@alternative,
            "less" = {
                z <- -statistic(object1, type = "standardized")
                zp <- -t((support(object2, raw = TRUE) - mu) / sigma)
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                zp <- t((support(object2, raw = TRUE) - mu) / sigma)
            },
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                zp <- abs(t((support(object2, raw = TRUE) - mu) / sigma))
            }
        )

        ## reorder simulations using the (decreasing) test statistics
        o <- order(z, decreasing = TRUE) # largest z first
        zp <- zp[, o, drop = FALSE]

        ## unadjusted p-values
        pu <- rowMeans(t(zp) %GE% z[o])

        ## permutation distribution
        p <- lapply(seq_len(ncol(zp)), function(i) {
            zp_i <- zp[, i]
            vapply(unique(zp_i), function(x, t) mean(x %GE% t), NA_real_,
                   x = zp_i)
        })

        ## discreteness adjustment (see Westfall and Wolfinger, 1997)
        RET <- rep.int(1 - bonferroni, length(z)) # zeros (ones) for Bonferroni (Sidak)
        for (i in 1:length(pu)) {
            qq <- if (stepdown) i else 1 # 'i' => successively smaller subsets
            for (q in qq:length(p)) {
                x <- p[[q]][p[[q]] <= pu[i]] # below eq. 2
                if (length(x) > 0) {
                    RET[i] <- if (bonferroni) RET[i] + max(x) # eq. 4
                              else RET[i] * (1 - max(x)) # eq. 2
                }
            }
        }
        RET <- if (!bonferroni) pmin.int(1 - RET, 1) else pmin.int(RET, 1)
        for (i in 2:length(RET))
            RET[i] <- max(RET[i - 1], RET[i]) # enforce monotonicity

        RET <- matrix(RET[order(o)], nrow = nrow(z), ncol = ncol(z),
                      dimnames = dimnames(z))
        class(RET) <- "pvalue"
        attr(RET, "nresample") <- object2@nresample
        RET
    }
)

setMethod("marginal",
    signature = list("MaxTypeIndependenceTest", "missing"),
    definition = function(object1, object2, stepdown, bonferroni, ...) {
        callGeneric(object1@statistic, object1@distribution, stepdown, bonferroni, ...)
    }
)


### unadjusted p-values
setGeneric("unadjusted",
    function(object1, object2, ...) {
        standardGeneric("unadjusted")
    }
)

setMethod("unadjusted",
    signature = list("MaxTypeIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, ...) {
        z <- statistic(object1, type = "standardized")
        RET <- switch(object1@alternative,
                   "less"      = pnorm(z),
                   "greater"   = 1 - pnorm(z),
                   "two.sided" = 2 * pmin.int(pnorm(z), 1 - pnorm(z))
               )

        RET <- matrix(RET, nrow = nrow(z), ncol = ncol(z),
                      dimnames = dimnames(z))
        class(RET) <- "pvalue"
        RET
    }
)

setMethod("unadjusted",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, ...) {
        ## standardized observed and permuted test statistics
        mu <- expectation(object1)
        sigma <- sqrt(variance(object1))
        switch(object1@alternative,
            "less" = {
                z <- statistic(object1, type = "standardized")
                RET <- (support(object2, raw = TRUE) - mu) / sigma
                RET <- rowMeans(RET %LE% as.vector(z))
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                RET <- (support(object2, raw = TRUE) - mu) / sigma
                RET <- rowMeans(RET %GE% as.vector(z))
            },
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                RET <- abs(support(object2, raw = TRUE) - mu) / sigma
                RET <- rowMeans(RET %GE% as.vector(z))
            }
        )

        RET <- matrix(RET, nrow = nrow(z), ncol = ncol(z),
                      dimnames = dimnames(z))
        class(RET) <- "pvalue"
        attr(RET, "nresample") <- object2@nresample
        RET
    }
)

setMethod("unadjusted",
    signature = list("MaxTypeIndependenceTest", "missing"),
    definition = function(object1, object2, ...) {
        callGeneric(object1@statistic, object1@distribution, ...)
    }
)


### compute p-values under subset pivotality (Westfall, 1997)
npmcp <- function(object) {

    ## extract from object
    y <- object@statistic@y[[1]]
    x <- object@statistic@x[[1]]
    ytrafo <- object@statistic@ytrafo
    alternative <- object@statistic@alternative

    ## <FIXME> it is currently hard to ask a distribution object
    ## for its type (and arguments). It's a design bug.
    distribution <- object@call$distribution
    ## </FIXME>
    z <- switch(alternative,
             "less"      = statistic(object, type = "standardized"),
             "greater"   = -statistic(object, type = "standardized"),
             "two.sided" = -abs(statistic(object, type = "standardized"))
         )

    ## get contrast matrix from xtrans
    C <- attr(object@statistic@xtrans, "contrast")
    stopifnot(inherits(C, "matrix"))

    ## order test statistics, most "extreme" one comes first
    Corder <- C[order(z), , drop = FALSE]

    ## compute allowed subsets of hypotheses
    ## returns list consisting of lists (one for each rejection step of H0)
    ms <- multcomp:::maxsets(Corder)

    ## make sure 'object' isn't serialized along with 'foo'
    ## (otherwise parallel operation using snow clusters will be very slow)
    object <- NULL # was rm(object)
    ## alternatively we could pass all relevant objects to 'foo' and then
    ## associate it with the global environment instead:
    ## foo <- function(s, y, x, ytrafo, distribution, alternative) { ... }
    ## environment(foo) <- .GlobalEnv
    ## or simply define 'foo' out of 'npmcp'

    foo <- function(s) {
        Ctmp <- Corder[s, , drop = FALSE] # current allowed subset
        ## x levels in current subset
        xlev <- apply(Ctmp, MARGIN = 2, function(col) any(col != 0))

        it <- independence_test(y ~ x,
                                subset = x %in% names(xlev)[xlev], # relevant data subset
                                xtrafo = mcp_trafo(x = Ctmp),
                                ytrafo = ytrafo,
                                distribution = distribution,
                                alternative = alternative)
        pvalue(it)
    }

    RET <- vapply(ms, function(sub) # for every list of allowed subsets
        max(vapply(sub, foo, NA_real_)), NA_real_) # for every subset

    for (i in 2:length(RET))
        RET[i] <- max(RET[i - 1], RET[i]) # enforce monotonicity

    matrix(RET[rank(z)], dimnames = dimnames(z))
}
