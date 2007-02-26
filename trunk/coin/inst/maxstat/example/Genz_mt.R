
library("coin")
x <- sort(runif(50))
y <- rnorm(length(x))
mt <- maxstat_test(y ~ x)

pmt <- function(mt) {

    xt <- mt@statistic@xtrans
    if(!all(order(colSums(xt)) == 1:ncol(xt)))
        stop(sQuote("mt"), " does not have an ordered covariate")
    R <- cov2cor(covariance(mt))
    R1 <- solve(R)

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

system.time(p1 <- pmt(mt))

system.time(p2 <- pvalue(mt))

