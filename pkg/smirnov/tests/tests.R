
library("smirnov")

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

qks2 <- function(p, n.x, n.y, ...) {
    stat <- cumsum(c(0, rep(2 / (n.x + n.y), ceiling(n.x + n.y) / 2)))
    prb <- sapply(stat, pks2, n.x = n.x, n.y = n.y, ...)
    if (is.null(p)) return(list(stat = stat, prob = prb))
    sapply(p, function(u) min(stat[prb > u]))
}

dks2 <- function(q, n.x, n.y, ...) {
    cdf <- qks2(NULL, n.x, n.y, ...)
    pdf <- diff(c(0, cdf$prob, 1))
    if (is.null(q)) {
        attr(pdf, "support") <- cdf$stat
        return(pdf)
    }
    i <- match(cdf$stat, q)
    pdf[i]
}

### Schröer & Trenkler
psmirnov(3 / 7, 7, 5)
pks2(3 / 7, 7, 5)

system.time(p1 <- sapply(1:999/1000, function(q) pks2(q, 10, 20)))
system.time(p2 <- sapply(1:999/1000, function(q) psmirnov(q, 10, 20)))

all.equal(p1, p2)

library("coin")
gr <- rep(gl(2, 1), c(5, 7))
x <- c(1, 2, 2, 3, 3, 1, 2, 3, 3, 4, 5, 6)
mt <- maxstat_test(gr ~ x, minprob = 0, maxprob = 1)
pls <- coin:::MonteCarlo(x = mt@statistic@xtrans, y = mt@statistic@ytrans, 
                         weights = mt@statistic@weights, 
                         block = mt@statistic@block, nresample = 100000,
                         parallel = "no")

expectation(mt)

TS <- colSums(mt@statistic@xtrans)
N <- table(gr)
KS <- abs(TS / N[2] - pls * sum(N) / prod(N))
mKS <- apply(KS, 2, max)

mean(mKS >= 3 / 7 - sqrt(.Machine$double.eps))
