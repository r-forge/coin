
library("smirnov")
set.seed(290875)

### note: psmirnov(, obs = NULL) uses stats:::C_pSmirnov2x, obs = something
### trigger Schröer-Trenkler algo
all.equal(psmirnov(1:999/1000, 10, 20, obs = 1:30), psmirnov(1:999/1000, 10, 20))

### examples by Schröer & Trenkler (1995) without ties
all.equal(psmirnov(3 / 7, 7, 5, obs = 1:12), psmirnov(3 / 7, 7, 5))

all.equal(psmirnov(n.x = 3, n.y = 4, q = 1/2, obs = 1:7, lower.tail = FALSE),
          psmirnov(n.x = 3, n.y = 4, q = 1/2, lower.tail = FALSE))
psmirnov(n.x = 3, n.y = 4, q = 1/6, obs = 1:7, lower.tail = FALSE, two.sided = FALSE)
psmirnov(n.x = 4, n.y = 3, q = 1/2, obs = 1:7, lower.tail = FALSE, two.sided = FALSE)

all.equal(psmirnov(n.x = 5, n.y = 7, q = 3 / 7, obs = 1:12, lower.tail = FALSE),
          psmirnov(n.x = 5, n.y = 7, q = 3 / 7, obs = 1:12, lower.tail = FALSE))

all.equal(psmirnov(n.x = 300, n.y = 520, q = 1/8, obs = 1:820, lower.tail = FALSE),
          psmirnov(n.x = 300, n.y = 520, q = 1/8, lower.tail = FALSE))
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
all.equal(qsmirnov(1 - .05, n.x = 8, n.y = 6) * 8 * 6, 34)
all.equal(qsmirnov(1 - .01, n.x = 8, n.y = 6) * 8 * 6, 40)
all.equal(qsmirnov(1 - .001, n.x = 8, n.y = 6) * 8 * 6, 48)
all.equal(qsmirnov(1 - .05, n.x = 14, n.y = 10) * 140, 74)
all.equal(qsmirnov(1 - .01, n.x = 14, n.y = 10) * 140, 90)
all.equal(qsmirnov(1 - .001, n.x = 14, n.y = 10) * 140, 106)
### Schröer algo
all.equal(qsmirnov(1 - .05, n.x = 8, n.y = 6, obs = 1:14) * 8 * 6, 34)
all.equal(qsmirnov(1 - .01, n.x = 8, n.y = 6, obs = 1:14) * 8 * 6, 40)
all.equal(qsmirnov(1 - .001, n.x = 8, n.y = 6, obs = 1:14) * 8 * 6, 48)
all.equal(qsmirnov(1 - .05, n.x = 14, n.y = 10, obs = 1:24) * 140, 74)
all.equal(qsmirnov(1 - .01, n.x = 14, n.y = 10, obs = 1:24) * 140, 90)
all.equal(qsmirnov(1 - .001, n.x = 14, n.y = 10, obs = 1:24) * 140, 106)

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

### compare speed
n.x <- 1:5 * 20
n.y <- 1:4 * 20
q <- 2:9 / 10
x <- expand.grid(n.x = n.x, n.y = n.y, q = q)
st <- sapply(1:nrow(x), function(i) {
    s1 <- system.time(lp1 <- psmirnov(q = x[i, "q"], 
        n.x = x[i, "n.x"], n.y = x[i, "n.y"]))["user.self"]
    s2 <- system.time(lp2 <- psmirnov(q = x[i, "q"], 
        n.x = x[i, "n.x"], n.y = x[i, "n.y"], obs = 1:(x[i, "n.x"] + x[i, "n.y"])))["user.self"]
    ret <- c(isTRUE(all.equal(lp1, lp2)) + 0L, s1 = s1, s2 = s2)
})

summary(t(st))
