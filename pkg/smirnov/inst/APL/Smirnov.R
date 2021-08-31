
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

psmirnov <- function(q, n.x, n.y = length(obs) - n.x, obs = NULL, 
                     two.sided = TRUE, lower.tail = TRUE, log.p = FALSE) {

    if (is.numeric(q)) 
        q <- as.double(q)
    else stop("argument 'q' must be numeric")
    ret <- rep(0, length(q))
    ret[is.na(q) | q < 0 | q > 1] <- NA
    IND <- which(!is.na(ret))
    if (!length(IND)) return(p)
    if (n.x < 1) stop("not enough 'x' data")
    if (n.y < 1) stop("not enough 'y' data")
    n.x <- floor(n.x)
    n.y <- floor(n.y)
    N <- n.x + n.y

    TIES <- if (!is.null(obs))
        c(diff(sort(obs)) != 0, TRUE)
    else
        rep(TRUE, N)

    ### see stats/src/ks.c line 103ff
    q <- (0.5 + floor(as.double(q) * n.x * n.y - 1e-7)) / (n.x * n.y);

#    cat("Ties =", TIES, "\n")
    pfun <- function(q) {
        k <- diag <- 1
        u <- 0
#    cat("m n =", n.x, n.y, "\n")
        repeat {
#      cat("k =", k, "\n")
            u <- c(u, 1 + u[length(u)])
            v <- k - u
            diag_bit <- (u <= n.x) & (v <= n.y) & (u >= 0) & (v >= 0)
            u <- u[diag_bit]
            v <- v[diag_bit]
#      cat("u =", u, "\n")
#      cat("v =", v, "\n")
            d <- u / n.x - v / n.y
            if (two.sided) d <- abs(d)
            diag <- (c(diag, 0) + c(0, diag))[diag_bit]
#      cat("diag =", diag, "\n")
            if (TIES[k])
                diag <- diag * (q > d)
#      cat("diag =", diag, "\n")
            k <- k + 1
        if (N < k) break
        }
        diag
    }
    ret[IND] <- sapply(q[IND], pfun)
#    cat("final diag =", ret, "\n")
    logdenom <- lgamma(N + 1) - lgamma(n.x + 1) - lgamma(n.y + 1)
    if (log.p & lower.tail)
        return(log(ret) - logdenom)
    if (!log.p & lower.tail)
        return(ret / exp(logdenom))
    if (log.p & !lower.tail)
        return(log1p(-ret / exp(logdenom)))
    if (!log.p & !lower.tail)
        return(1 - ret / exp(logdenom))
}

psmirnov(n.x = 3, n.y = 4, q = 1/2, lower.tail = FALSE)
1 - pks2(n.x = 3, n.y = 4, q = 1/2)
psmirnov(n.x = 3, n.y = 4, q = 1/6, lower.tail = FALSE, two.sided = FALSE)
psmirnov(n.x = 4, n.y = 3, q = 1/2, lower.tail = FALSE, two.sided = FALSE)
psmirnov(n.x = 5, n.y = 7, q = 3 / 7, lower.tail = FALSE)
1 - pks2(n.x = 5, n.y = 7, q = 3 / 7)


psmirnov(obs = c(1, 2, 2, 3, 3, 1, 2, 3, 3, 4, 5, 6), n.x = 5, q = 3 / 7, lower.tail = FALSE)

psmirnov(n.x = 300, n.y = 520, q = 1/8, lower.tail = FALSE)
1 - pks2(n.x = 300, n.y = 520, q = 1/8)
psmirnov(obs = round(rnorm(820)), n.x = 300, n.y = 520, q = 1/8, lower.tail = FALSE)
