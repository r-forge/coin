### single-step maxT multiple testing procedure
setGeneric("singlestep",
    function(object1, object2, ...) {
        standardGeneric("singlestep")
    }
)

setMethod("singlestep",
    signature = list("MaxTypeIndependenceTestStatistic", "NullDistribution"),
    definition = function(object1, object2, ...) {
        ## reorder test statistics to ensure consistency with "global"/"step-down"
        switch(object1@alternative,
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                o <- order(z, decreasing = TRUE) # abs. largest z first
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                o <- order(z, decreasing = TRUE) # largest z first
            },
            "less" = {
                z <- statistic(object1, type = "standardized")
                o <- order(z)                    # smallest z first
            }
        )

        ## iterate over unique test statistics only and remap
        pq <- length(z)
        RET <- z[o]
        idx <- c(which(RET[-1L] %NE% RET[-pq]), pq) # unique z
        RET <- pvalue(object2, RET[idx], ...)

        matrix(rep.int(RET, diff(c(0L, idx)))[order(o)], # remapping
               nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
    }
)

setMethod("singlestep",
    signature = list("MaxTypeIndependenceTest", "missing"),
    definition = function(object1, object2, ...) {
        callGeneric(object1@statistic, object1@distribution, ...)
    }
)


### step-down maxT multiple testing procedure
setGeneric("stepdown",
    function(object1, object2, ...) {
        standardGeneric("stepdown")
    }
)

setMethod("stepdown",
    signature = list("MaxTypeIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, ...) {
        ## reorder upper and/or lower limits using test statistics
        switch(object1@alternative,
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                o <- order(z, decreasing = TRUE) # abs. largest z first
                pq <- length(z)
                upper <- z[o]
                lower <- -upper
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                o <- order(z, decreasing = TRUE) # largest z first
                pq <- length(z)
                upper <- z[o]
                lower <- rep.int(-Inf, pq)
            },
            "less" = {
                z <- statistic(object1, type = "standardized")
                o <- order(z)                    # smallest z first
                pq <- length(z)
                upper <- rep.int(Inf, pq)
                lower <- z[o]
            }
        )

        ## step-down based on multivariate normality
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

        matrix(1 - RET[order(o)], nrow = nrow(z), ncol = ncol(z),
               dimnames = dimnames(z))
    }
)

setMethod("stepdown",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, ...) {
        ## standardized observed and permuted test statistics
        mu <- expectation(object1)
        sigma <- sqrt(variance(object1))
        switch(object1@alternative,
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                zp <- abs(t((support(object2, raw = TRUE) - mu) / sigma))
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                zp <- t((support(object2, raw = TRUE) - mu) / sigma)
            },
            "less" = {
                z <- -statistic(object1, type = "standardized")
                zp <- -t((support(object2, raw = TRUE) - mu) / sigma)
            }
        )

        ## reorder simulations using (increasing) test statistics
        o <- order(z) # smallest z first
        zp <- zp[, o, drop = FALSE]

        ## algorithm 2.8 (Free Step-Down Resampling Method) in
        ## Westfall & Young (1993), page 66 _using standardized
        ## statistics instead of p-values_!
        if (ncol(zp) > 1) {
            for (j in 2:ncol(zp))
                zp[, j] <- pmax.int(zp[, j], zp[, j - 1])
        }
        RET <- rowMeans(t(zp) %GE% z[o])
        for (i in (length(RET) - 1):1)
            RET[i] <- max(RET[i], RET[i + 1]) # enforce monotonicity, page 67

        matrix(RET[order(o)], nrow = nrow(z), ncol = ncol(z),
               dimnames = dimnames(z))
    }
)

setMethod("stepdown",
    signature = list("MaxTypeIndependenceTest", "missing"),
    definition = function(object1, object2, ...) {
        callGeneric(object1@statistic, object1@distribution, ...)
    }
)


### Adjusted marginal p-values taking discreteness into account in the
### permutation case (Westfall and Wolfinger, 1997)
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
                   "two.sided" = 2 * pmin.int(pnorm(z), 1 - pnorm(z)),
                   "greater"   = 1 - pnorm(z),
                   "less"      = pnorm(z)
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

        matrix(RET, nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
    }
)

setMethod("marginal",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, stepdown, bonferroni, ...) {
        ## standardized observed and permuted test statistics
        mu <- expectation(object1)
        sigma <- sqrt(variance(object1))
        switch(object1@alternative,
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                zp <- abs(t((support(object2, raw = TRUE) - mu) / sigma))
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                zp <- t((support(object2, raw = TRUE) - mu) / sigma)
            },
            "less" = {
                z <- -statistic(object1, type = "standardized")
                zp <- -t((support(object2, raw = TRUE) - mu) / sigma)
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

        ## discreteness adjustment
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

        matrix(RET[order(o)], nrow = nrow(z), ncol = ncol(z),
               dimnames = dimnames(z))
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
                   "two.sided" = 2 * pmin.int(pnorm(z), 1 - pnorm(z)),
                   "greater"   = 1 - pnorm(z),
                   "less"      = pnorm(z)
               )

        matrix(RET, nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
    }
)

setMethod("unadjusted",
    signature = list("MaxTypeIndependenceTestStatistic", "ApproxNullDistribution"),
    definition = function(object1, object2, ...) {
        ## standardized observed and permuted test statistics
        mu <- expectation(object1)
        sigma <- sqrt(variance(object1))
        switch(object1@alternative,
            "two.sided" = {
                z <- abs(statistic(object1, type = "standardized"))
                zp <- abs(support(object2, raw = TRUE) - mu) / sigma
            },
            "greater" = {
                z <- statistic(object1, type = "standardized")
                zp <- (support(object2, raw = TRUE) - mu) / sigma
            },
            "less" = {
                z <- -statistic(object1, type = "standardized")
                zp <- -(support(object2, raw = TRUE) - mu) / sigma
            }
        )

        ## unadjusted p-values
        matrix(rowMeans(zp %GE% as.vector(z)),
               nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
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
             "two.sided" = -abs(statistic(object, type = "standardized")),
             "greater" = -statistic(object, type = "standardized"),
             "less" = statistic(object, type = "standardized")
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
