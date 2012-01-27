
library("coin")
library("multcomp")
set.seed(290875)

### Jonckhere-Terpsra test mit joint ranking.

fun <- function(nsim, k = 0) {

    p <- numeric(nsim)

    for (i in 1:nsim) {

        x <- gl(4, 10)
        y <- rnorm(length(x), mean = as.numeric(x) * k)

        ff <- function(x) {
            K <- contrMat(table(x), "Tukey")[,x]
            as.vector(rep(1, nrow(K)) %*%K)
        }
        p[i] <- pvalue(independence_test(y ~ x, ytrafo = rank,
           xtrafo = function(data) trafo(data, factor_trafo = ff)))
    }
    mean(p <= 0.05)
}

power <- sapply(k <- seq(from = 0, to = 1, by = 0.1), fun, nsim = 100)
plot(k, power)


