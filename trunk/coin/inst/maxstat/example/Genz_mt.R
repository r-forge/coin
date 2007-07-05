
library("coin")
set.seed(290875)
x <- sort(runif(50))
y <- rnorm(length(x))
mt <- maxstat_test(y ~ x)

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


system.time(p1 <- pmt(mt))

system.time(p2 <- worsley(mt))

system.time(p3 <- pvalue(mt))

x <- ordered(cut(1:10000, breaks = seq(from = 0, to = 10000, by = 100)))
y <- rnorm(10000)
mt <- maxstat_test(y ~ x)

mta <- maxstat_test(y ~ x, distribution = approximate(20000))


g <- seq(from = 2.3, to = 4, by = 0.1)
pex <- sapply(g, function(i) 1 - pmt(mt, i))
pap<- sapply(g, function(i) 1 - worsley(mt, i))
ex <- sapply(g, function(i) pperm(mta, i))

save(g, pex, pap, ex, file = "approx.rda")

plot(g, pex, type = "l", lty = 1, xlab = expression(c), ylab = expression(P(T[max] <= c)))
lines(g, pap, lty = 2)
lines(g, pap, lty = 3)
legend("bottomright", lty = c(1, 2, 3), 
       legend = c("Exact Asymptotic", "Improved Bonferroni", "Exact"), bty = "n")
