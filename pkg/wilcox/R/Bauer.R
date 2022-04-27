
.Bauer_two_exact <- function(x, y, conf.level, digits.rank, alternative) {

    ## Exact confidence interval for the location parameter
    ## mean(x) - mean(y) in the two-sample case (cf. the
    ## one-sample case).

    n.x <- length(x)
    n.y <- length(y)
    r <- c(x, y)
    r <- rank(if(is.finite(digits.rank)) signif(r, digits.rank) else r)

    alpha <- 1 - conf.level
    diffs <- sort(outer(x, y, "-"))
    cint <-
        switch(alternative,
            "two.sided" = {
                qu <- qwilcox(alpha/2, n.x, n.y)
                if(qu == 0) qu <- 1
                ql <- n.x*n.y - qu
                achieved.alpha <- 2*pwilcox(trunc(qu)-1,n.x,n.y)
                c(diffs[qu], diffs[ql + 1])
             },
             "greater" = {
                 qu <- qwilcox(alpha, n.x, n.y)
                 if(qu == 0) qu <- 1
                 achieved.alpha <- pwilcox(trunc(qu)-1,n.x,n.y)
                 c(diffs[qu], +Inf)
             },
             "less" = {
                 qu <- qwilcox(alpha, n.x, n.y)
                 if(qu == 0) qu <- 1
                 ql <- n.x*n.y - qu
                 achieved.alpha <- pwilcox(trunc(qu)-1,n.x,n.y)
                 c(-Inf, diffs[ql + 1])
         })
    if (achieved.alpha-alpha > alpha/2) {
        warning("Requested conf.level not achievable")
        conf.level <- 1 - achieved.alpha
    }
    attr(cint, "conf.level") <- conf.level
    ESTIMATE <- c("difference in location" = median(diffs))
    list(cint, ESTIMATE)
}

.Bauer_two_asym <- function(x, y, conf.level, digits.rank, 
                            alternative, correct, tol.root) {

    ## Asymptotic confidence interval for the location
    ## parameter mean(x) - mean(y) in the two-sample case
    ## (cf. one-sample case).
    ##
    ## Algorithm not published, for a documentation see the
    ## one-sample case.

    n.x <- length(x)
    n.y <- length(y)

    alpha <- 1 - conf.level
    mumin <- min(x) - max(y)
    mumax <- max(x) - min(y)
    W <- function(d) { ## also fn (x, y, n.x, n.y, correct, alternative)
        dr <- c(x - d, y)
        dr <- rank(if(is.finite(digits.rank)) signif(dr, digits.rank) else dr)
        NTIES.CI <- table(dr)
        dz <- sum(dr[seq_along(x)]) - n.x * (n.x + 1) / 2 - n.x * n.y / 2
        CORRECTION.CI <-
            if(correct) {
                switch(alternative,
                    "two.sided" = sign(dz) * 0.5,
                    "greater" = 0.5,
                    "less" = -0.5)
                } else 0
        SIGMA.CI <- sqrt((n.x * n.y / 12) * ((n.x + n.y + 1)
                         - sum(NTIES.CI^3 - NTIES.CI)
                             / ((n.x + n.y) * (n.x + n.y - 1))))
        if (SIGMA.CI == 0)
            warning(
                "cannot compute confidence interval when all observations are tied",
                call.=FALSE)
        (dz - CORRECTION.CI) / SIGMA.CI
    }
 
    wdiff <- function(d, zq) W(d) - zq
    Wmumin <- W(mumin)
    Wmumax <- W(mumax)
    root <- function(zq) {
        ## in extreme cases we need to return endpoints,
        ## e.g.  wilcox.test(1, 2:60, conf.int=TRUE)
        f.lower <- Wmumin - zq
        if(f.lower <= 0) return(mumin)
        f.upper <- Wmumax - zq
        if(f.upper >= 0) return(mumax)
            uniroot(wdiff, lower=mumin, upper=mumax,
                    f.lower = f.lower, f.upper = f.upper,
                    tol = tol.root, zq = zq)$root
        }
    cint <- switch(alternative,
                   "two.sided" = {
                       l <- root(zq = qnorm(alpha/2, lower.tail = FALSE))
                       u <- root(zq = qnorm(alpha/2))
                       c(l, u)
                   },
                   "greater" = {
                       l <- root(zq = qnorm(alpha, lower.tail = FALSE))
                       c(l, +Inf)
                   },
                   "less" = {
                       u <- root(zq = qnorm(alpha))
                       c(-Inf, u)
                   })
    attr(cint, "conf.level") <- conf.level
    correct <- FALSE # for W(): no continuity correction for estimate
    ESTIMATE <- c("difference in location" =
                   uniroot(W, lower=mumin, upper=mumax,
                           tol = tol.root)$root)
    list(cint, ESTIMATE)
}
