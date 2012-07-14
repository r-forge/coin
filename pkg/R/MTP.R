

### single step maxT multiple testing procedure
singlestep <- function(object, ...) {

    ### reorder test statistics to ensure consistency with "global"/"step-down"
    switch(object@statistic@alternative,
           "two.sided" = {
               ts <- abs(statistic(object, "standardized"))
               ots <- rank(-ts) # original order
               rts <- order(ts, decreasing = TRUE)}, # abs. largest ts first
           "greater" = {
               ts <- statistic(object, "standardized")
               ots <- rank(-ts) # original order
               rts <- order(ts, decreasing = TRUE)}, # largest ts first
           "less" = {
               ts <- statistic(object, "standardized")
               ots <- rank(ts) # original order
               rts <- order(ts) # smallest ts first
           })

    ### iterate over unique test statistics only and remap
    pq <- length(ts)
    tsrts <- ts[rts]
    idx <- c(which(tsrts[-1L] != tsrts[-pq]), pq)
    uts <- tsrts[idx] # unique ts
    ret <- 1 - sapply(uts, pperm, object = object, ...)

    ret <- matrix(ret[ots], nrow = nrow(ts), ncol = ncol(ts))
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}

### algorithm 2.8 (Free Step-Down Resampling Method) in
### Westfall & Young (1993), page 66 _using standardized
### statistics instead of p-values_!
### <FIXME>
rsdmaxT <- function(pls, ts) {

    ### reorder simulations using (increasing) test statistics
    ots <- rank(ts) # original order
    rts <- order(ts) # smallest ts first
    q <- pls[, rts, drop = FALSE]

    ### algorithm 2.8 (Free Step-Down Resampling Method) in
    ### Westfall & Young (1993), page 66 _using standardized
    ### statistics instead of p-values_!
    if (ncol(q) > 1) {
        for (j in 2:ncol(q))
            q[,j] <- pmax(q[,j], q[,j-1])
    }
    ret <- rowMeans(GE(t(q), ts[rts]))
    for (i in (length(ret) - 1):1)
        ret[i] <- max(ret[i], ret[i + 1]) # enforce monotonicity, page 67

    ret <- matrix(ret[ots], nrow = nrow(ts), ncol = ncol(ts))
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}

### step-down using the asymptotic distribution
asdmaxT <- function(object) {

    ### reorder upper and/or lower limits using test statistics
    switch(object@statistic@alternative,
           "two.sided" = {
               ts <- abs(statistic(object, "standardized"))
               pq <- length(ts)
               ots <- rank(-ts) # original order
               rts <- order(ts, decreasing = TRUE) # abs. largest ts first
               upper <- ts[rts]
               lower <- -upper},
           "greater" = {
               ts <- statistic(object, "standardized")
               pq <- length(ts)
               ots <- rank(-ts) # original order
               rts <- order(ts, decreasing = TRUE) # largest ts first
               upper <- ts[rts]
               lower <- rep(-Inf, pq)},
           "less" = {
               ts <- statistic(object, "standardized")
               pq <- length(ts)
               ots <- rank(ts) # original order
               rts <- order(ts) # smallest ts first
               upper <- rep(Inf, pq)
               lower <- ts[rts]})

    ### correlation matrix
    corr <- cov2cor(covariance(object))

    ### step-down based on multivariate normality
    ret <- numeric(pq)
    ret[1] <- pmv(lower = lower[1], upper = upper[1],
                  mean = rep(0, pq), corr = corr)
    if (pq > 1) {
        for (i in 2:pq) {
            j <- rank(rts)[1] # reindexing needed in each step
            corr <- corr[-j, -j]
            rts <- rts[-1]
            ret[i] <- min(ret[i - 1],
                          pmv(lower = lower[i], upper = upper[i],
                              mean = rep(0, length(rts)), corr = corr))
        }
    }

    ret <- matrix(1 - ret[ots], nrow = nrow(ts), ncol = ncol(ts))
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}

### stepdown maxT multiple testing procedure
stepdown <- function(object, ...) {

    if (!(extends(class(object), "MaxTypeIndependenceTest")))
        stop(sQuote("object"), " is not of class ",
             sQuote("MaxTypeIndependenceTest"))

    if (extends(class(object@distribution), "AsymptNullDistribution"))
        ret <- asdmaxT(object)
    else {
        ## raw simulation results, scores have been handled already
        pls <- support(object, raw = TRUE)

        ## standardize
        dcov <- sqrt(variance(object))
        expect <- expectation(object)
        switch(object@statistic@alternative,
               "two.sided" = {
                   pls <- abs(t((pls - expect) / dcov))
                   ts <- abs(statistic(object, "standardized"))},
               "greater" = {
                   pls <- t((pls - expect) / dcov)
                   ts <- statistic(object, "standardized")},
               "less" = {
                   pls <- -t((pls - expect) / dcov)
                   ts <- -(statistic(object, "standardized"))})

        ret <- rsdmaxT(pls, ts)
    }

    ret
}

### Discrete permutation method (Westfall & Wolfinger, 1997, AmStat 51, 3-8)
discrete <- function(object, method = c("Bonferroni", "Sidak",
                                        "Bonferroni-Holm", "Sidak-Holm"), ...) {

    ### <FIXME> this should be possible when the _exact_ marginal
    ### distributions are available
    ### </FIXME>

    if (!(extends(class(object), "MaxTypeIndependenceTest") &&
          extends(class(object@distribution), "ApproxNullDistribution")))
        stop(sQuote("object"), " is not of class ",
             sQuote("MaxTypeIndependenceTest"),
             " or distribution was not approximated via Monte-Carlo")

    method <- match.arg(method)

    bonferroni <- pmatch(method, "Bonferroni-Holm", nomatch = 0)
    stepdown <- method %in% c("Bonferroni-Holm", "Sidak-Holm")

    ### raw simulation results, scores have been handled already
    pls <- support(object, raw = TRUE)

    ### standardize
    dcov <- sqrt(variance(object))
    expect <- expectation(object)
    switch(object@statistic@alternative,
           "two.sided" = {
               pls <- abs(t((pls - expect) / dcov))
               ts <- abs(statistic(object, "standardized"))},
           "greater" = {
               pls <- t((pls - expect) / dcov)
               ts <- statistic(object, "standardized")},
           "less" = {
               pls <- -t((pls - expect) / dcov)
               ts <- -(statistic(object, "standardized"))})

    ### reorder simulations using the (decreasing) test statistics
    ots <- rank(-ts) # original order
    rts <- order(ts, decreasing = TRUE) # largest ts first
    pls <- pls[, rts, drop = FALSE]

    ### unadjusted p-values
    pvals <- rowMeans(GE(t(pls), ts[rts]))

    ### permutation distribution
    foo <- function(x, t) mean(GE(x, t))
    p <- vector(mode = "list", length = ncol(pls))
    for (i in 1:ncol(pls)) {
        ux <- unique(pls[, i])
        p[[i]] <- sapply(ux, foo, x = pls[, i])
    }

    ### discrete adjustment
    ret <- rep(1 - bonferroni, length(ts)) # zeros (ones) for Bonferroni (Sidak)
    for (i in 1:length(pvals)) {
        qq <- if (stepdown) i else 1 # 'i' => successively smaller subsets
        for (q in qq:length(p)) {
            x <- p[[q]][p[[q]] <= pvals[i]] # below eq. 2
            if (length(x) > 0) {
                ret[i] <- if(bonferroni) ret[i] + max(x) # eq. 4
                           else ret[i] * (1 - max(x)) # eq. 2
            }
        }
    }
    if (!bonferroni)
        ret <- 1 - ret
    ret <- pmin(ret, 1)
    for (i in 2:length(ret))
        ret[i] <- max(ret[i - 1], ret[i]) # enforce monotonicity

    ret <- matrix(ret[ots], nrow = nrow(ts), ncol = ncol(ts))
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}

#####################################
## Westfall (1997) method in coin ###
#####################################
## Basic code for npmcp() taken from pqfunctions.R in package
## multcomp

### cf. mcp(x = "Tukey") in multcomp
mcp_trafo <- function(...) {

    args <- list(...)
    stopifnot(length(args) == 1)

    ret <- function(data) {

        x <- data[[names(args)]]
        stopifnot(is.factor(x))
        C <- args[[1]]
        if (is.character(C)) {
            C <- contrMat(table(x), C)
        } else {
            stopifnot(is.matrix(C))
            stopifnot(ncol(C) == nlevels(x))
            if (is.null(colnames(C)))
                colnames(C) <- levels(x)
            attr(C, "type") <- "User-defined"
            class(C) <- c("contrMat", "matrix")
        }
        ret <- trafo(data,
                     factor_trafo = function(x) model.matrix(~ x - 1) %*% t(C))
        attr(ret, "contrast") <- C
        ret
    }
    ret
}

### compute p-values under subset pivotality
npmcp <- function(object) {

    ### extract from object
    y <- object@statistic@y[[1]]
    x <- object@statistic@x[[1]]
    ytrafo <- object@statistic@ytrafo
    alternative <- object@statistic@alternative

    ### <FIXME> it is currently hard to ask a distribution object
    ### for its type (and arguments). Its a design bug.
    distribution <- object@call$distribution
    ### use default value
    if (is.null(distribution))
        distribution <- eval(formals(eval(object@call[[1]]))$distribution)[1]
    ### </FIXME>
    stand_tstat <- statistic(object, type = "standardized")
    tstat <- switch(alternative,
                    "less" = stand_tstat,
                    "greater" = -stand_tstat,
                    "two.sided" = -abs(stand_tstat))

    # get contrast matrix from xtrans
    C <- attr(object@statistic@xtrans, "contrast")
    stopifnot(inherits(C, "matrix"))

    # order test statistics, most "extreme" one comes first
    Corder <- C[order(tstat), , drop = FALSE]

    # compute allowed subsets of hypotheses
    # returns list consisting of lists (one for each rejection step of H0)
    ms <- multcomp:::maxsets(Corder)

    foo <- function(s) {
        Ctmp <- Corder[s, , drop = FALSE] # current allowed subset
        # x levels in current subset
        xlev <- apply(Ctmp, MARGIN = 2, function(col) any(col != 0))

        it <- independence_test(y ~ x,
                                subset = x %in% names(xlev)[xlev], # relevant data subset
                                xtrafo = mcp_trafo(x = Ctmp),
                                ytrafo = ytrafo,
                                distribution = distribution,
                                alternative = alternative)
        pvalue(it)
    }

    p <- sapply(ms, function(sub) # for every list of allowed subsets
        max(sapply(sub, foo))) # for every subset

    for (i in 2:length(p))
        p[i] <- max(p[i-1], p[i]) # forces pvalue monotonicity

    ret <- matrix(p[rank(tstat)])
    attr(ret, "dimnames") <- attr(tstat, "dimnames")
    return(ret)
}

### unadjusted p-values
unadjusted <- function(object, ...) {

    if (extends(class(object@distribution), "AsymptNullDistribution")) {
        ts <- statistic(object, "standardized")
        ret <- switch(object@statistic@alternative,
                      "two.sided" = 2 * pmin(pnorm(ts), 1 - pnorm(ts)),
                      "greater"   = 1 - pnorm(ts),
                      "less"      = pnorm(ts))
    } else {
        ## raw simulation results, scores have been handled already
        pls <- support(object, raw = TRUE)

        ## standardize
        dcov <- sqrt(variance(object))
        expect <- expectation(object)
        switch(object@statistic@alternative,
               "two.sided" = {
                   pls <- abs((pls - expect) / dcov)
                   ts <- abs(statistic(object, "standardized"))},
               "greater" = {
                   pls <- (pls - expect) / dcov
                   ts <- statistic(object, "standardized")},
               "less" = {
                   pls <- -(pls - expect) / dcov
                   ts <- -(statistic(object, "standardized"))})

        ## unadjusted p-values
        ret <- matrix(rowMeans(GE(pls, as.vector(ts))),
                      nrow = nrow(ts), ncol = ncol(ts))
        rownames(ret) <- rownames(ts)
        colnames(ret) <- colnames(ts)
    }

    ret
}
