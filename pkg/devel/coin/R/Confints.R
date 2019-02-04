setGeneric("confint_location",
    function(object1, object2, ...) {
        standardGeneric("confint_location")
    }
)

setMethod("confint_location",
    signature = list("ScalarIndependenceTestStatistic", "NullDistribution"),
    definition = function(object1, object2, level, ...) {
        ## <FIXME> drop unused levels!
        if (!is_2sample(object1))
            warning(sQuote("object1"), " does not represent a two-sample problem")
        ## </FIXME>

        if (nlevels(object1@block) != 1L || !is_unity(object1@weights))
            stop("cannot compute confidence intervals with blocks or weights")

        alternative <- object1@alternative

        if(!((length(level) == 1L)
            && is.finite(level)
            && (level > 0)
            && (level < 1)))
            stop("level must be a single number between 0 and 1")

        scores <- object1@y[[1L]]
        groups <- object1@xtrans[, 1L]

        ## raw data
        x <- sort(scores[groups > 0])
        y <- sort(scores[groups < 1])
        alpha <- 1 - level

        foo <- function(x, d) x - d

        ## explicitly compute all possible steps
        steps <- c()
        for (lev in levels(object1@block)) {
            thisblock <- (object1@block == lev)
            ytmp <- sort(split(scores[thisblock], groups[thisblock])[[1L]])
            xtmp <- sort(split(scores[thisblock], groups[thisblock])[[2L]])
            steps <- c(steps, as.vector(outer(xtmp, ytmp, foo)))
        }
        steps <- sort(unique(steps))

        ## computes the statistic under the alternative 'd'
        fse <- function(d)
            sum(object1@ytrafo(data.frame(c(foo(x, d), y)))[seq_along(x)])

        ## we need to compute the statistics just to the right of
        ## each step
        ds <- diff(steps)
        justright <- min(abs(ds[abs(ds) > sqrt_eps])) / 2
        jumps <- vapply(steps + justright, fse, NA_real_)

        ## determine if the statistics are in- or decreasing
        ## jumpsdiffs <- diff(jumps)
        increasing <- all(diff(jumps[c(1L, length(jumps))]) > 0)
        decreasing <- all(diff(jumps[c(1L, length(jumps))]) < 0)

        ## this is safe
        if (!(increasing || decreasing))
            stop("cannot compute confidence intervals:",
                 "the step function is not monotone")

        cci <- function(alpha) {
            ## the quantiles:
            ## we reject iff
            ##
            ##   STATISTIC <  qlower OR
            ##   STATISTIC >= qupper
            ##
            qlower <- drop(qperm(object2, alpha / 2) *
                             sqrt(variance(object1)) + expectation(object1))
            qupper <- drop(qperm(object2, 1 - alpha / 2) *
                             sqrt(variance(object1)) + expectation(object1))

            ## Check if the statistic exceeds both quantiles first.
            if (qlower < min(jumps) || qupper > max(jumps)) {
                warning("Cannot compute confidence intervals")
                return(c(NA, NA))
            }

            if (increasing) {
                ##
                ##  We do NOT reject for all steps with
                ##
                ##     STATISTICS >= qlower AND
                ##     STATISTICS < qupper
                ##
                ##  but the open right interval ends with the
                ##  step with STATISTIC == qupper
                ##
                ci <- c(min(steps[qlower %LE% jumps]),
                        min(steps[jumps > qupper]))
            } else {
                ##
                ##  We do NOT reject for all steps with
                ##
                ##     STATISTICS >= qlower AND
                ##     STATISTICS < qupper
                ##
                ##  but the open left interval ends with the
                ##  step with STATISTIC == qupper
                ##
                ci <- c(min(steps[jumps %LE% qupper]),
                        min(steps[jumps < qlower]))
            }
            ci
        }

        cint <- switch(alternative,
                    "two.sided" = cci(alpha),
                    "greater"   = c(cci(alpha * 2)[1L], Inf),
                    "less"      = c(-Inf, cci(alpha * 2)[2L])
                )
        attr(cint, "conf.level") <- level

        ## was: median(steps) which will not work for blocks etc.
        u <- jumps - expectation(object1)
        sgr <- ifelse(decreasing, min(steps[u %LE% 0]), max(steps[u %LE% 0]))
        sle <- ifelse(decreasing, min(steps[u < 0]), min(steps[u > 0]))

        ESTIMATE <- mean(c(sle, sgr), na.rm = TRUE)
        names(ESTIMATE) <- "difference in location"

        list(conf.int = cint, estimate = ESTIMATE)
    }
)

setMethod("confint_location",
    signature = list("ScalarIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, level, ...) {
        ## <FIXME> drop unused levels!
        if (!is_2sample(object1))
            warning(sQuote("object1"), " does not represent a two-sample problem")
        ## </FIXME>

        if (nlevels(object1@block) != 1L || !is_unity(object1@weights))
            stop("cannot compute confidence intervals with blocks or weights")

        alternative <- object1@alternative

        if(!((length(level) == 1L)
            && is.finite(level)
            && (level > 0)
            && (level < 1)))
            stop("level must be a single number between 0 and 1")

        scores <- object1@y[[1L]]
        groups <- object1@xtrans[, 1L]

        ## raw data
        x <- sort(scores[groups > 0])
        y <- sort(scores[groups < 1])
        alpha <- 1 - level

        foo <- function(x, d) x - d

        ## approximate the steps
        ## Here we search the root of the function 'fsa' on the set
        ## c(mumin, mumax).
        ##
        ## This returns a value from c(mumin, mumax) for which
        ## the standardized statistic is equal to the
        ## quantile zq.  This means that the statistic is not
        ## within the critical region, and that implies that '
        ## is a confidence limit for the median.

        fsa <- function(d, zq) {
           STAT <- sum(object1@ytrafo(data.frame(c(foo(x, d), y)))[seq_along(x)])
           (STAT - expectation(object1)) / sqrt(variance(object1)) - zq
        }

        mumin <- min(x) - max(y)
        mumax <- max(x) - min(y)

        ccia <- function(alpha) {
            ## Check if the statistic exceeds both quantiles
            ## first: otherwise 'uniroot' won't work anyway
            statu <- fsa(mumin, zq = qperm(object2, alpha / 2))
            statl <- fsa(mumax, zq = qperm(object2, 1 - alpha / 2))
            if (sign(statu) == sign(statl)) {
                warning(paste("Samples differ in location:",
                              "Cannot compute confidence set,",
                              "returning NA"))
                return(c(NA, NA))
            }
            u <- uniroot(fsa, c(mumin, mumax),
                         zq = qperm(object2, alpha / 2), ...)$root
            l <- uniroot(fsa, c(mumin, mumax),
                         zq = qperm(object2, 1 - alpha / 2), ...)$root
            ## The process of the statistics does not need to be
            ## increasing: sort is ok here.
            sort(c(u, l))
        }

        cint <- switch(alternative,
                    "two.sided" = ccia(alpha),
                    "greater"   = c(ccia(alpha * 2)[1L], Inf),
                    "less"      = c(-Inf, ccia(alpha * 2)[2L])
                )
        attr(cint, "conf.level") <- level

        ## Check if the statistic exceeds both quantiles first.
        statu <- fsa(mumin, zq = 0)
        statl <- fsa(mumax, zq = 0)
        if (sign(statu) == sign(statl)) {
            ESTIMATE <- NA
            warning("Cannot compute estimate, returning NA")
        } else
            ESTIMATE <- uniroot(fsa, c(mumin, mumax), zq = 0, ...)$root
        names(ESTIMATE) <- "difference in location"

        list(conf.int = cint, estimate = ESTIMATE)
    }
)


setGeneric("confint_scale",
    function(object1, object2, ...) {
        standardGeneric("confint_scale")
    }
)

setMethod("confint_scale",
    signature = list("ScalarIndependenceTestStatistic", "NullDistribution"),
    definition = function(object1, object2, level, ...) {
        ## <FIXME> drop unused levels!
        if (!is_2sample(object1))
            warning(sQuote("object1"), " does not represent a two-sample problem")
        ## </FIXME>

        if (nlevels(object1@block) != 1L || !is_unity(object1@weights))
            stop("cannot compute confidence intervals with blocks or weights")

        alternative <- object1@alternative

        if(!((length(level) == 1L)
            && is.finite(level)
            && (level > 0)
            && (level < 1)))
            stop("level must be a single number between 0 and 1")

        scores <- object1@y[[1L]]
        groups <- object1@xtrans[, 1L]

        ## raw data
        x <- sort(scores[groups > 0])
        y <- sort(scores[groups < 1])
        alpha <- 1 - level

        foo <- function(x, d) x / d

        ## explicitly compute all possible steps
        steps <- c()
        for (lev in levels(object1@block)) {
            thisblock <- (object1@block == lev)
            ytmp <- sort(split(scores[thisblock], groups[thisblock])[[1L]])
            xtmp <- sort(split(scores[thisblock], groups[thisblock])[[2L]])
            ratio <-  outer(xtmp, ytmp, "/")
            aratio <- ratio[ratio >= 0]
            steps <- c(steps, aratio)
        }
        steps <- sort(unique(steps))

        ## computes the statistic under the alternative 'd'
        fse <- function(d)
            sum(object1@ytrafo(data.frame(c(foo(x, d), y)))[seq_along(x)])

        ## we need to compute the statistics just to the right of
        ## each step
        ds <- diff(steps)
        justright <- min(abs(ds[abs(ds) > sqrt_eps])) / 2
        jumps <- vapply(steps + justright, fse, NA_real_)

        ## determine if the statistics are in- or decreasing
        ## jumpsdiffs <- diff(jumps)
        increasing <- all(diff(jumps[c(1L, length(jumps))]) > 0)
        decreasing <- all(diff(jumps[c(1L, length(jumps))]) < 0)

        ## this is safe
        if (!(increasing || decreasing))
            stop("cannot compute confidence intervals:",
                 "the step function is not monotone")

        cci <- function(alpha) {
            ## the quantiles:
            ## we reject iff
            ##
            ##   STATISTIC <  qlower OR
            ##   STATISTIC >= qupper
            ##
            qlower <- drop(qperm(object2, alpha / 2) *
                             sqrt(variance(object1)) + expectation(object1))
            qupper <- drop(qperm(object2, 1 - alpha / 2) *
                             sqrt(variance(object1)) + expectation(object1))

            ## Check if the statistic exceeds both quantiles first.
            if (qlower < min(jumps) || qupper > max(jumps)) {
                warning("Cannot compute confidence intervals")
                return(c(NA, NA))
            }

            if (increasing) {
                ##
                ##  We do NOT reject for all steps with
                ##
                ##     STATISTICS >= qlower AND
                ##     STATISTICS < qupper
                ##
                ##  but the open right interval ends with the
                ##  step with STATISTIC == qupper
                ##
                ci <- c(min(steps[qlower %LE% jumps]),
                        min(steps[jumps > qupper]))
            } else {
                ##
                ##  We do NOT reject for all steps with
                ##
                ##     STATISTICS >= qlower AND
                ##     STATISTICS < qupper
                ##
                ##  but the open left interval ends with the
                ##  step with STATISTIC == qupper
                ##
                ci <- c(min(steps[jumps %LE% qupper]),
                        min(steps[jumps < qlower]))
            }
            ci
        }

        cint <- switch(alternative,
                    "two.sided" = cci(alpha),
                    "greater"   = c(cci(alpha * 2)[1L], Inf),
                    "less"      = c(0, cci(alpha * 2)[2L])
                )
        attr(cint, "conf.level") <- level

        u <- jumps - expectation(object1)
        sgr <- ifelse(decreasing, min(steps[u %LE% 0]), max(steps[u %LE% 0]))
        sle <- ifelse(decreasing, min(steps[u < 0]), min(steps[u > 0]))

        ESTIMATE <- mean(c(sle, sgr), na.rm = TRUE)
        names(ESTIMATE) <- "ratio of scales"

        list(conf.int = cint, estimate = ESTIMATE)
    }
)

setMethod("confint_scale",
    signature = list("ScalarIndependenceTestStatistic", "AsymptNullDistribution"),
    definition = function(object1, object2, level, ...) {
        ## <FIXME> drop unused levels!
        if (!is_2sample(object1))
            warning(sQuote("object1"), " does not represent a two-sample problem")
        ## </FIXME>

        if (nlevels(object1@block) != 1L || !is_unity(object1@weights))
            stop("cannot compute confidence intervals with blocks or weights")

        alternative <- object1@alternative

        if(!((length(level) == 1L)
            && is.finite(level)
            && (level > 0)
            && (level < 1)))
            stop("level must be a single number between 0 and 1")

        scores <- object1@y[[1L]]
        groups <- object1@xtrans[, 1L]

        ## raw data
        x <- sort(scores[groups > 0])
        y <- sort(scores[groups < 1])
        alpha <- 1 - level

        foo <- function(x, d) x / d

        ## approximate the steps
        ## Here we search the root of the function 'fsa' on the set
        ## c(mumin, mumax).
        ##
        ## This returns a value from c(mumin, mumax) for which
        ## the standardized statistic is equal to the
        ## quantile zq.  This means that the statistic is not
        ## within the critical region, and that implies that '
        ## is a confidence limit for the median.

        fsa <- function(d, zq) {
           STAT <- sum(object1@ytrafo(data.frame(c(foo(x, d), y)))[seq_along(x)])
           (STAT - expectation(object1)) / sqrt(variance(object1)) - zq
        }

        srangepos <- NULL
        srangeneg <- NULL
        if (any(x > 0) && any(y > 0))
            srangepos <-
                c(min(x[x > 0], na.rm = TRUE) / max(y[y > 0], na.rm = TRUE),
                  max(x[x > 0], na.rm = TRUE) / min(y[y > 0], na.rm = TRUE))
        if (any(x %LE% 0) && any(y < 0))
            srangeneg <-
                c(min(x[x %LE% 0], na.rm = TRUE) / max(y[y < 0], na.rm = TRUE),
                  max(x[x %LE% 0], na.rm = TRUE) / min(y[y < 0], na.rm = TRUE))
        if (any(is.infinite(c(srangepos, srangeneg)))) {
            stop(paste("Cannot compute asymptotic confidence",
                       "set or estimator"))
        }

        mumin <- range(c(srangepos, srangeneg), na.rm = FALSE)[1L]
        mumax <- range(c(srangepos, srangeneg), na.rm = FALSE)[2L]

        ccia <- function(alpha) {
            ## Check if the statistic exceeds both quantiles
            ## first: otherwise 'uniroot' won't work anyway
            statu <- fsa(mumin, zq = qperm(object2, alpha / 2))
            statl <- fsa(mumax, zq = qperm(object2, 1 - alpha / 2))
            if (sign(statu) == sign(statl)) {
                warning(paste("Samples differ in location:",
                              "Cannot compute confidence set,",
                              "returning NA"))
                return(c(NA, NA))
            }
            u <- uniroot(fsa, c(mumin, mumax),
                         zq = qperm(object2, alpha / 2), ...)$root
            l <- uniroot(fsa, c(mumin, mumax),
                         zq = qperm(object2, 1 - alpha / 2), ...)$root
            ## The process of the statistics does not need to be
            ## increasing: sort is ok here.
            sort(c(u, l))
        }

        cint <- switch(alternative,
                    "two.sided" = ccia(alpha),
                    "greater"   = c(ccia(alpha*2)[1L], Inf),
                    "less"      = c(0, ccia(alpha*2)[2L])
                )
        attr(cint, "conf.level") <- level

        ## Check if the statistic exceeds both quantiles first.
        statu <- fsa(mumin, zq = 0)
        statl <- fsa(mumax, zq = 0)
        if (sign(statu) == sign(statl)) {
            ESTIMATE <- NA
            warning("Cannot compute estimate, returning NA")
        } else
            ESTIMATE <- uniroot(fsa, c(mumin, mumax), zq = 0, ...)$root
        names(ESTIMATE) <- "ratio of scales"

        list(conf.int = cint, estimate = ESTIMATE)
    }
)


### CI for a binomial parameter
confint_binom <-
    function(x, n, level = 0.95, method = c("exact", "mid-p"), tol = eps)
{
    method <- match.arg(method)

    RET <- if (x >= 0 && x <= n) {
               alpha <- 1 - level
               if (method == "exact") {
                   ## exact Clopper-Pearson interval
                   c(if (x == 0) 0 else qbeta(    alpha / 2, x    , n - x + 1),
                     if (x == n) 1 else qbeta(1 - alpha / 2, x + 1, n - x    ))
               } else {
                   ## mid-p interval (see Berry and Armitage, 1995)
                   if (x == 0) {
                       c(0, 1 - alpha^(1 / n))
                   } else if (x == n) {
                       c(alpha^(1 / n), 1)
                   } else {
                       f <- function(p, a)
                           ## 0.5 * dbinom(...) + pbinom(..., lower.tail = FALSE)
                           mean(pbinom(c(x, x - 1), n, p, lower.tail = FALSE)) - a
                       c(uniroot(f, c(0, 1), a =     alpha / 2, tol = tol)$root,
                         uniroot(f, c(0, 1), a = 1 - alpha / 2, tol = tol)$root)
                   }
               }
           } else {
               stop(sQuote("x"), " must be larger or equal to 0 and",
                    " smaller or equal to ", sQuote("n"))
           }
    attr(RET, "conf.level") <- level
    RET
}


###
### Currently unused
###
simconfint_location <- function(object, level = 0.95,
    approx = FALSE, ...) {

    if (!(is_Ksample(object@statistic) &&
        inherits(object, "MaxTypeIndependenceTest")))
        stop(sQuote("object"), " is not an object of class ",
             sQuote("MaxTypeIndependenceTest"),
             " representing a K sample problem")

    xtrans <- object@statistic@xtrans
    if (!all(apply(xtrans, 2L, function(x) all(x %in% c(-1, 0, 1)))))
        stop("Only differences are allowed as contrasts")

    estimate <- c()
    lower <- c()
    upper <- c()

    ## transform max(abs(x))-type distribution into a
    ## distribution symmetric around zero
    nnd <- object@distribution
    nnd@q <- function(p) {
        pp <- p
        if (p > 0.5)
            pp <- 1 - pp
        q <- qperm(object@distribution, 1 - pp)
        if (p < 0.5)
            -q
        else
            q
    }

    for (i in seq_len(ncol(xtrans))) {
        thisset <- abs(xtrans[, i]) > 0
        ip <- new("IndependenceProblem",
                  object@statistic@x[thisset, , drop = FALSE],
                  object@statistic@y[thisset, , drop = FALSE],
                  object@statistic@block[thisset])

        itp <- independence_test(ip, teststat = "scalar",
            distribution = "asymptotic", alternative = "two.sided",
            yfun = object@statistic@ytrafo, ...)

        ci <- confint_location(itp@statistic, nnd,
                               level = level, approx =approx, ...)
        estimate <- c(estimate, ci$estimate)
        lower <- c(lower, ci$conf.int[1L])
        upper <- c(upper, ci$conf.int[2L])
    }
    RET <- data.frame(Estimate = estimate, lower = lower, upper = upper)
    colnames(RET)[2L:3L] <-
        paste(c((1 - level) / 2, 1 - (1 - level) / 2) * 100, "%")
    rownames(RET) <- colnames(object@statistic@xtrans)
    attr(RET, "conf.level") <- level
    RET
}
