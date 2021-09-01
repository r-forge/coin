
library("smirnov")
set.seed(290875)

### 2-sample two-sided KS distribution from stats, with w/o ties only
pks2 <- function(q, n.x, n.y, exact = TRUE, tol = 1e-06) {
    if (any(q < 0)) return(0)
    if (any(q > 1)) return(1)
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
        p[IND] <- sapply(q[IND], function(s) 
                         .Call(stats:::C_pSmirnov2x, s, floor(n.x), floor(n.y)))
    } else {
        n <- n.x * n.y/(n.x + n.y)
        p[IND] <- .Call(stats:::C_pKS2, p = sqrt(n) * x[IND], tol)
    }
    return(p)
}

all.equal(psmirnov(1:999/1000, 10, 20), pks2(1:999/1000, 10, 20))

### examples by Schröer & Trenkler (1995) without ties
all.equal(psmirnov(3 / 7, 7, 5), pks2(3 / 7, 7, 5))

all.equal(psmirnov(n.x = 3, n.y = 4, q = 1/2, lower.tail = FALSE),
          1 - pks2(n.x = 3, n.y = 4, q = 1/2))
psmirnov(n.x = 3, n.y = 4, q = 1/6, lower.tail = FALSE, two.sided = FALSE)
psmirnov(n.x = 4, n.y = 3, q = 1/2, lower.tail = FALSE, two.sided = FALSE)

all.equal(psmirnov(n.x = 5, n.y = 7, q = 3 / 7, lower.tail = FALSE),
          1 - pks2(n.x = 5, n.y = 7, q = 3 / 7))

all.equal(psmirnov(n.x = 300, n.y = 520, q = 1/8, lower.tail = FALSE),
          1 - pks2(n.x = 300, n.y = 520, q = 1/8))
psmirnov(obs = round(rnorm(820)), n.x = 300, n.y = 520, 
         q = 1/8, lower.tail = FALSE)

library("coin")
### tied example by Schröer & Trenkler (1995); use coin to 
### approximate conditional p-value
gr <- rep(gl(2, 1), c(5, 7))
x <- c(1, 2, 2, 3, 3, 1, 2, 3, 3, 4, 5, 6)
mt <- maxstat_test(gr ~ x, minprob = 0, maxprob = 1)
pls <- coin:::MonteCarlo(x = mt@statistic@xtrans, y = mt@statistic@ytrans, 
                         weights = mt@statistic@weights, 
                         block = mt@statistic@block, nresample = 1e7,
                         parallel = "no")
TS <- colSums(mt@statistic@xtrans)
N <- table(gr)
KS <- abs(TS / N[2] - pls * sum(N) / prod(N))
mKS <- apply(KS, 2, max)
mean(mKS >= 3 / 7 - sqrt(.Machine$double.eps))
psmirnov(obs = x, n.x = 5, q = 3 / 7, lower.tail = FALSE)

### check quantiles
### Kim & Jennrich (1973) in Selected Tables in Mathematical Statistics, Vol 1
### (ed. Harter) Institute of Mathematical Statistics, page 129
all.equal(qsmirnov(1 - .05, 8, 6) * 8 * 6, 34)
all.equal(qsmirnov(1 - .01, 8, 6) * 8 * 6, 40)
all.equal(qsmirnov(1 - .001, 8, 6) * 8 * 6, 48)
all.equal(qsmirnov(1 - .05, 14, 10) * 140, 74)
all.equal(qsmirnov(1 - .01, 14, 10) * 140, 90)
all.equal(qsmirnov(1 - .001, 14, 10) * 140, 106)

### without ties
q <- qsmirnov(1:9/10, 5, 7)
p <- psmirnov(q, 5, 7)
all.equal(qsmirnov(p, 5, 7), q)

### with ties
obs <- c(1, 2, 2, 3, 3, 1, 2, 3, 3, 4, 5, 6)
q <- qsmirnov(1:9/10, 5, 7, obs = obs)
p <- psmirnov(q, 5, 7, obs = obs)
all.equal(qsmirnov(p, 5, 7, obs = obs), q)

### without ties
q <- qsmirnov(1:9/10, 5, 7, two.sided = FALSE)
p <- psmirnov(q, 5, 7, two.sided = FALSE)
all.equal(qsmirnov(p, 5, 7, two.sided = FALSE), q)

### with ties
obs <- c(1, 2, 2, 3, 3, 1, 2, 3, 3, 4, 5, 6)
q <- qsmirnov(1:9/10, 5, 7, obs = obs, two.sided = FALSE)
p <- psmirnov(q, 5, 7, obs = obs, two.sided = FALSE)
all.equal(qsmirnov(p, 5, 7, obs = obs, two.sided = FALSE), q)

