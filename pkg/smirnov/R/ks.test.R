
.Smirnov_sim <- function(tr, tc, B, twosided)
    .Call("Smirnov_sim", as.integer(tr), as.integer(tc), 
          as.integer(B), as.integer(twosided), package = "smirnov")

psmirnov_exact <- function(q, sizes, z = NULL, two.sided = TRUE,
                           lower.tail = TRUE) {
    if(!is.null(z)) {
        z <- (diff(sort(z)) != 0)
        z <- if(any(z))
            c(0L, z, 1L)
        else
            NULL
    }
    .Call("psmirnov_exact", q, sizes[1L], sizes[2L], z,
          two.sided, lower.tail, package = "smirnov")
}

psmirnov <- function(q, m, n = length(z) - m, z = NULL, 
                     two.sided = TRUE, exact = TRUE, 
                     simulate = FALSE, B = 2000,
                     lower.tail = TRUE, log.p = FALSE) {

    ##
    ## Distribution function Prob(D < q) for the two-sample Smirnov test statistic
    ##
    ## D = sup_c | ECDF_x(c) - ECDF_y(c) | 	(two.sided)
    ##
    ## D = sup_c ( ECDF_x(c) - ECDF_y(c) ) 	(!two.sided)
    ##
    ## Implementation translated from APL code in Appendix C.2.3 of
    ##
    ##     Gunar Schröer, Computergestützte statistische Inferenz am Beispiel der
    ##     Kolmogorov-Smirnov Tests, Diplomarbeit Universität Osnabrück, 1991
    ##
    ##     see also
    ##     
    ##     Gunar Schröer and Dietrich Trenkler (1995), Exact and Randomization
    ##     Distributions of Kolmogorov-Smirnov Tests for Two or Three
    ##     Samples, Computational Statistics & Data Analysis, 20, 185--202
    ##
    ## Original APL code (slightly adapted for recent dyalog interpreter)
    ##
    ## ⍝ Get a free trial version of dyalog from https://www.dyalog.com/
    ## ⍝ Tested with version 18.0.40684
    ##      ∇R←SP DMN_TIES c;m;n;TIES;k;diag;u;v;diag_bit;d
    ## [1]  (m n)←⍴¨SP
    ## [2]  SP←SP[⍋SP←∊SP]
    ## [3]  TIES←(¯1↓SP≠1⌽SP),1
    ## [4]  k←1
    ## [5]  diag←1
    ## [6]  u←0
    ## [7]  ⍝ LOOP:
    ## [8]  :Repeat
    ## [9]    u←u,1+¯1↑u
    ## [10]   v←k-u
    ## [11]   diag_bit←(u≤m)∧(v≤n)∧(u≥0)∧v≥0
    ## [12]   ⍝(u v)←diag_bit/¨u v
    ## [13]   u←diag_bit/u
    ## [14]   v←diag_bit/v
    ## [15]   d←|(u÷m)-v÷n
    ## [16]   diag←diag_bit/(diag,0)+0,diag
    ## [17]   diag←∊(1 0=TIES[k])/(diag×c>d)(diag)
    ## [18]   k←k+1
    ## [19] :Until  (m+n)<k
    ## [20] ⍝ →LOOP×~(m+n)≥k←k+1
    ## [21] R←1-diag÷m!m+n
    ## [22] ∇
    ##      S←(1 2 3)(4 5 6 7)
    ##      prob←S DMN_TIES 1÷2
    ##      prob
    ## 0.6571428571
    ##      S←(1 2 3 4 5)(6 7 8 9 10 11 12)
    ##      prob←S DMN_TIES 3÷7
    ##      prob
    ## 0.5454545455
    ##      S←(1 2 2 3 3)(1 2 3 3 4 5 6)
    ##      prob←S DMN_TIES 3÷7
    ##      prob
    ## 0.2424242424

    if (is.numeric(q)) 
        q <- as.double(q)
    else stop("argument 'q' must be numeric")
    ret <- rep(0, length(q))
    ret[is.na(q) | q < -1 | q > 1] <- NA
    IND <- which(!is.na(ret))
    if (!length(IND)) return(ret)
    if (m < 1) stop("not enough 'x' data")
    if (n < 1) stop("not enough 'y' data")
    n.x <- floor(m)
    n.y <- floor(n)
    N <- n.x + n.y
    n <- n.x * n.y / (n.x + n.y)

    if (!exact) {
        if (simulate) {
            if (is.null(z)) {
                rt <- rep(1, length = N)
            } else {
                rt <- table(z)
            }
            Dsim <- .Smirnov_sim(rt, c(n.x, n.y), B, two.sided)
            ### need P(D < q)
            ret[IND] <- ecdf(Dsim)(q - sqrt(.Machine$double.eps)) 
            if (log.p & lower.tail)
                return(log(ret))
            if (!log.p & lower.tail)
                return(ret)
            if (log.p & !lower.tail)
                return(log1p(-ret))
            if (!log.p & !lower.tail)
                return(1 - ret)
        }
        ## <FIXME> let m and n be the min and max of the sample
        ## sizes, respectively.  Then, according to Kim and Jennrich
        ## (1973), if m < n/10, we should use the
        ## * Kolmogorov approximation with c.c. -1/(2*n) if 1 < m < 80;
        ## * Smirnov approximation with c.c. 1/(2*sqrt(n)) if m >= 80.
        if (two.sided) {
            ret[IND] <- .Call(stats:::C_pKS2, p = sqrt(n) * q[IND], tol = 1e-6)
            ### note: C_pKS2(0) = NA but Prob(D < 0) = 0
            ret[q < .Machine$double.eps] <- 0
        } else {
            ret[IND] <- 1 - exp(- 2 * n * q^2)
        }
        if (log.p & lower.tail)
            return(log(ret))
        if (!log.p & lower.tail)
            return(ret)
        if (log.p & !lower.tail)
            return(log1p(-ret))
        if (!log.p & !lower.tail)
            return(1 - ret)
    }

    pfun <- function(q)
        psmirnov_exact(q, sizes = c(n.x, n.y), z = z, two.sided = two.sided,
                       lower.tail = lower.tail)

    ret[IND] <- sapply(q[IND], pfun)
    if (any(is.na(ret[IND]))) {
        warning("computation of exact probability failed, returning Monte Carlo approximation")
        return(psmirnov(q = q, m = n.x, n = n.y, z = z, 
                        two.sided = two.sided, exact = FALSE, 
                        simulate = TRUE, B = B,
                        lower.tail = lower.tail, log.p = log.p))
    }

    if (log.p) ret <- log(ret)
    return(ret)
}

qsmirnov <- function(p, m, n = length(z) - m, z = NULL, two.sided = TRUE, ...) {

    n.x <- floor(m)
    n.y <- floor(n)
    if (n.x * n.y < 1e4) {
        ### note: The support is also OK in case of ties
        stat <- sort(unique(c(outer(0:n.x/n.x, 0:n.y/n.y, "-"))))
    } else {
        stat <- (-1e4):1e4 / (1e4 + 1)
    }
    if (two.sided) stat <- abs(stat)
    prb <- psmirnov(stat, m = n.x, n = n.y, z = z, 
                    two.sided = two.sided, log.p = FALSE, ...)
    if (is.null(p)) return(list(stat = stat, prob = prb))
    if (is.numeric(p)) 
        p <- as.double(p)
    else stop("argument 'p' must be numeric")
    ret <- rep(0, length(p))
    ret[is.na(p) | p < 0 | p > 1] <- NA
    IND <- which(!is.na(ret))
    if (!length(IND)) return(ret)
    ret[IND] <- sapply(p[IND], function(u) min(stat[prb >= u]))
    ret
}

ks.test <- function(x, ...) UseMethod("ks.test")

ks.test.default <-
    function(x, y, ..., alternative = c("two.sided", "less", "greater"),
             exact = NULL, simulate.p.value = FALSE, B = 2000)
{
    alternative <- match.arg(alternative)
    DNAME <- deparse1(substitute(x))
    x <- x[!is.na(x)]
    n <- length(x)
    if(n < 1L)
        stop("not enough 'x' data")
    PVAL <- NULL

    ### ordered variables can be treated as numeric ones as ties are handled
    ### now
    if (is.ordered(y)) y <- unclass(y)

    if(is.numeric(y)) { ## two-sample case
        args <- list(...)
        if (length(args) > 0L) 
            warning("Parameter(s) ", paste(names(args), collapse = ", "), " ignored")
        DNAME <- paste(DNAME, "and", deparse1(substitute(y)))
        y <- y[!is.na(y)]
        n.x <- as.double(n)             # to avoid integer overflow
        n.y <- length(y)
        if(n.y < 1L)
            stop("not enough 'y' data")
        if(is.null(exact))
            exact <- (n.x * n.y < 10000)
        if (!simulate.p.value) {
            METHOD <- paste(c("Asymptotic", "Exact")[exact + 1L], 
                            "two-sample Smirnov test")
        } else {
            METHOD <- "Monte-Carlo two-sample Smirnov test"
        }
        TIES <- FALSE
        n <- n.x * n.y / (n.x + n.y)
        w <- c(x, y)
        z <- cumsum(ifelse(order(w) <= n.x, 1 / n.x, - 1 / n.y))
        if(length(unique(w)) < (n.x + n.y)) {
            z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
            TIES <- TRUE
            if (!exact & !simulate.p.value)
                warning("p-value will be approximate in the presence of ties")
        }
        STATISTIC <- switch(alternative,
                            "two.sided" = max(abs(z)),
                            "greater" = max(z),
                            "less" = - min(z))
        nm_alternative <- switch(alternative,
                                 "two.sided" = "two-sided",
                                 "less" = "the CDF of x lies below that of y",
                                 "greater" = "the CDF of x lies above that of y")
        z <- NULL
        if (TIES)
            z <- w
        PVAL <- switch(alternative,
            "two.sided" = psmirnov(STATISTIC, m = n.x, n = n.y, z = w, exact = exact, 
                                   simulate = simulate.p.value, B = B,
                                   lower.tail = FALSE),
            "less" = psmirnov(STATISTIC, m = n.x, n = n.y, z = w,  exact = exact,
                              simulate = simulate.p.value, B = B,
                              two.sided = FALSE, lower.tail = FALSE),
            "greater" = psmirnov(STATISTIC, m = n.x, n = n.y, z = w, exact = exact,
                                 simulate = simulate.p.value, B = B,
                                 two.sided = FALSE, lower.tail = FALSE))
    } else { ## one-sample case
        if(is.character(y)) # avoid matching anything in this function
            y <- get(y, mode = "function", envir = parent.frame())
        if(!is.function(y))
            stop("'y' must be numeric or a function or a string naming a valid function")
        TIES <- FALSE
        if(length(unique(x)) < n) {
            warning("ties should not be present for the Kolmogorov-Smirnov test")
            TIES <- TRUE
        }
        if(is.null(exact)) exact <- (n < 100) && !TIES
        METHOD <- paste(c("Asymptotic", "Exact")[exact + 1L], 
                        "one-sample Kolmogorov-Smirnov test")
        x <- y(sort(x), ...) - (0 : (n-1)) / n
        STATISTIC <- switch(alternative,
                            "two.sided" = max(c(x, 1/n - x)),
                            "greater" = max(1/n - x),
                            "less" = max(x))
        if(exact) {
            PVAL <- 1 - if(alternative == "two.sided")
                .Call(stats:::C_pKolmogorov2x, STATISTIC, n)
            else {
                pkolmogorov1x <- function(x, n) {
                    ## Probability function for the one-sided
                    ## one-sample Kolmogorov statistics, based on the
                    ## formula of Birnbaum & Tingey (1951).
                    if(x <= 0) return(0)
                    if(x >= 1) return(1)
                    j <- seq.int(from = 0, to = floor(n * (1 - x)))
                    1 - x * sum(exp(lchoose(n, j)
                                    + (n - j) * log(1 - x - j / n)
                                    + (j - 1) * log(x + j / n)))
                }
                pkolmogorov1x(STATISTIC, n)
            }
        } else {
            PVAL <- if(alternative == "two.sided")
                        1 - .Call(stats:::C_pKS2, sqrt(n) * STATISTIC, 
                                                  tol = 1e-6)
                    else exp(- 2 * n * STATISTIC^2)
        }
        nm_alternative <-
            switch(alternative,
                   "two.sided" = "two-sided",
                   "less" = "the CDF of x lies below the null hypothesis",
                   "greater" = "the CDF of x lies above the null hypothesis")
    }

    names(STATISTIC) <- switch(alternative,
                               "two.sided" = "D",
                               "greater" = "D^+",
                               "less" = "D^-")

    ## fix up possible overshoot (PR#14671)
    PVAL <- min(1.0, max(0.0, PVAL))
    RVAL <- list(statistic = STATISTIC,
                 p.value = PVAL,
                 alternative = nm_alternative,
                 method = METHOD,
                 data.name = DNAME,
                 data = list(x = x, y = y), # for confband.ks.test
                 exact = exact)             # same
    class(RVAL) <- c("ks.test", "htest")
    return(RVAL)
}

ks.test.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    oneSample <- FALSE
    if (length(attr(terms(formula[-2L]), "term.labels")) != 1L)
        if (formula[[3]] == 1L)
            oneSample <- TRUE
        else
            stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    rname <- names(mf)[1L]
    DNAME <- paste(names(mf), collapse = " by ") # works in all cases
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    if (! oneSample) { # 2-sample Smirnov test
        g <- factor(mf[[-response]])
        if(nlevels(g) != 2L)
            stop("grouping factor must have exactly 2 levels")
        DATA <- setNames(split(mf[[response]], g), c("x", "y"))
        y <- do.call("ks.test", c(DATA, list(...)))
        y$alternative <- gsub("x", levels(g)[1L], y$alternative)
        y$alternative <- gsub("y", levels(g)[2L], y$alternative)
        y$response <- rname   # for plot.confband.ks.test
        y$groups <- levels(g) # for plot.confband.ks.test
    }
    else { # 1-sample Kolmogorov-Smirnov test
        respVar <- mf[[response]]
        DATA <- list(x = respVar)
        y <- do.call("ks.test", c(DATA, list(...)))
        y$alternative <- gsub("x", DNAME, y$alternative)
    }
    y$data.name <- DNAME
    y
}
