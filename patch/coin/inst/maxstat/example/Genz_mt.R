
library("maxstat")
library("coin")
set.seed(290875)

pmt <- function(mt, tstat = NULL) {

    xt <- mt@statistic@xtrans
    if(!all(order(colSums(xt)) == 1:ncol(xt)))
        stop(sQuote("mt"), " does not have an ordered covariate")
    R <- cov2cor(covariance(mt))
    R1 <- solve(R)

    if (is.null(tstat))
        tstat <- statistic(mt)
    g <- seq(from = -tstat, to = tstat, length = 512)
    g2 <- g^2
    dg <- (g[2] - g[1]) / sqrt(2 * pi)
    ggt <- g %*% t(g)
    r <- diag(R1) / 2
    R1m <- R1 * (-1)

    phi <- function(k)
        exp(R1m[k,k-1] * ggt - r[k] * g2)

    f <- rep(1, length(g))
    for (i in nrow(R):2)
        f <- colSums(phi(i) * f) * dg

    1 - (1 / sqrt(det(R))) * sum(f * exp(-r[1] * g2)) * dg
}

worsley <- function(mt, tstat = NULL) {

    cr <- cov2cor(covariance(mt))
    if (is.null(tstat))
        tstat <- statistic(mt)
    k <- nrow(cr)
    if (tstat < 0) stop("tstat must be greater zero")
    up <- 2 * pnorm(-tstat)

    bp <- numeric(k-1)
    for (i in 2:k)
        bp[i-1] <-
            2 * pmvnorm(c(tstat, tstat), c(Inf, Inf), corr = cr[c(i-1,i), c(i-1,i)]) +
            2 * pmvnorm(c(-Inf, tstat), c(-tstat, Inf), corr = cr[c(i-1,i), c(i-1,i)])
    k * up - sum(bp)
}



x <- ordered(cut(1:100, breaks = seq(from = 0, to = 100, by = 5)))
y <- rank(rnorm(100))
mt <- maxstat_test(y ~ x)

system.time(p1 <- pmt(mt))

system.time(p2 <- worsley(mt))

system.time(p3 <- pvalue(mt))

system.time(mta <- maxstat_test(y ~ x, distribution = approximate(30000)))

g <- seq(from = 2.3, to = 4, by = 0.1)
pex <- sapply(g, function(i) 1 - pmt(mt, i))
pap<- sapply(g, function(i) 1 - worsley(mt, i))
ex <- sapply(g, function(i) pperm(mta, i))
pl <- sapply(g, function(i) 1- pLausen92(i))

save(g, pex, pap, ex, pl, file = "approx.rda")

plot(g, ex, type = "l", lty = 1, xlab = expression(c), ylab = expression(P(T[max] <= c)))
lines(g, pex, lty = 2)
lines(g, pap, lty = 3)
lines(g, pl, lty = 4)
legend("bottomright", lty = 1:4, 
       legend = c("Monte-Carlo", "Asymptotic", "Improved Bonferroni", "Brownian Motion"), 
       bty = "n")
