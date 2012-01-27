
library("coin")
library("multtest")

### data generating process: two groups
dgp <- function(n = 100, mu = rep(0, 2)) {

    m <- length(mu)
    group <- gl(2, n/2)
    x <- matrix(rnorm(n * m), ncol = m)
    mu <- rbind(matrix(0, nrow = n/2, ncol = m), 
                matrix(rep(mu, rep(n/2, m)), ncol = m))
    X <- x + mu
    data.frame(group = group, X)
}

### a simple two sample problem
mydf <- dgp(n = 200, mu = 0.25)
it <- independence_test(X ~ group, data = mydf, distribution = approximate())

### p-values
pvalue(it)
t.test(X ~ group, data = mydf, var.equal = TRUE)$p.value
MTP(t(as.matrix(mydf[["X"]])), Y = mydf[["group"]],
    test = "t.twosamp.equalvar")@adjp
mt.maxT(t(as.matrix(mydf[["X"]])), classlabel = as.numeric(mydf[["group"]]) -1)

set.seed(290875)
p <- numeric(1000)
for (i in 1:length(p)) {
    mydf <- dgp(n = 200, mu = 0)
    p[i] <- MTP(t(as.matrix(mydf[["X"]])), Y = mydf[["group"]],
                test = "t.twosamp.equalvar")@adjp
}

mean(p < 0.05)
