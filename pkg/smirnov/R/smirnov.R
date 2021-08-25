

### 2-sample two-sided KS distribution
pks2 <- function(q, n.x, n.y, exact = TRUE, tol = 1e-06) {
    if (q < 0) return(0)
    if (q > 1) return(1)
    if (is.numeric(q)) 
        q <- as.double(q)
    else stop("argument 'q' must be numeric")
    p <- rep(0, length(q))
    p[is.na(q)] <- NA
    IND <- which(!is.na(q))
    if (!length(IND)) return(p)
    if (n.x < 1) stop("not enough 'x' data")
    if (n.y < 1) stop("not enough 'y' data")
    if (exact) {
        p[IND] <- sapply(q[IND], function(s) .Call(stats:::C_pSmirnov2x, s, floor(n.x), floor(n.y)))
    } else {
        n <- n.x * n.y/(n.x + n.y)
        p[IND] <- .Call(stats:::C_pKS2, p = sqrt(n) * x[IND], tol)
    }
    return(p)
}


psmirnov <- function(q, n.x, n.y, obs = NULL, two.sided = TRUE,
                     lower.tail = TRUE, log.p = FALSE) {

    ###
    ### Distribution function Prob(D < q) for the Smirnov test statistic
    ###
    ### D = sup_c | ECDF_x(c) - ECDF_y(c) | 	(two.sided)
    ###
    ### D = sup_c ( ECDF_x(c) - ECDF_y(c) ) 	(!two.sided)
    ###

    if (length(q) > 1)
        return(sapply(q, psmirnov, n.x = n.x, n.y = n.y, obs = obs, 
                                   lower.tail = lower.tail, log.p = log.p))

    n.x <- as.integer(n.x)
    n.y <- as.integer(n.y)
    N <- n.x + n.y
    if (length(n.x) > 1 || length(n.y) > 1 || n.x < 2 || n.x < 2)
        stop(sQuote("n.x"), "and/or", sQuote("n.y"), "misspecified")
    umat <- matrix(0, nrow = n.x + 1L, ncol = n.y + 1L)
    ix <- 0:n.x
    jy <- 0:n.y

    ### see stats/src/ks.c line 103ff
    q <- (0.5 + floor(as.double(q) * n.x * n.y - 1e-7)) / (n.x * n.y);

    abs2 <- if (two.sided) abs else function(x) x

    if (!is.null(obs)) {
        if (length(obs) != N)
            stop(sQuote("length(obs)"), "is not equal to", sQuote("n.x + n.y"))
        z <- sort(obs)
        r <- rank(z, ties.method = "max")
        eqz <- diff(z) < sqrt(.Machine$double.eps)
        for (i in ix) {
            for (j in jy) {
                cmat <- 0L
                if (abs2(i / n.x - j / n.y) < q) 
                    cmat <- 1L
                k <- i + j
                if (r[j + 1] > n.y + i) {
                    cat("i: ", i, "j: ", j, "\n")
                    break;
                }
                if (k > 0 && k < N && eqz[k])
                    cmat <- 1L
                fct <- 1L
                if (i > 0 && j > 0) 
                    fct <- umat[i + 1, j] + umat[i, j + 1]
                umat[i + 1, j + 1] <- cmat * fct
            }
        }            
    } else {
        for (i in ix) {
            for (j in jy) {
                cmat <- 0L
                if (abs2(i / n.x - j / n.y) < q) 
                    cmat <- 1L
                fct <- 1    
                if (i > 0 && j > 0) 
                    fct <- umat[i + 1, j] + umat[i, j + 1]
                umat[i + 1, j + 1] <- cmat * fct
            }
        }            
    }

print(umat)

    term <- lgamma(n.x + n.y + 1) - lgamma(n.x + 1) - lgamma(n.y + 1)
    ret <- umat[n.x + 1, n.y + 1]
    if (lower.tail && log.p)
        return(log(ret) - term)
    if (lower.tail && !log.p)
        return(ret / exp(term))
    if (!lower.tail && !log.p)
        return(1 - ret / exp(term))
    return(log1p(-ret / exp(term)))
}

### Schröer & Trenkler (1994)
1 - psmirnov(3 / 7, 5, 7)

### Schröer & Trenkler
1 - psmirnov(3 / 7, 5, 7, obs = c(1, 2, 2, 3, 3, 1, 2, 3, 3, 4, 5, 6))

### correct
1 - 600 / choose(12, 5)
### incorrect
1 - 603 / choose(12, 5)

