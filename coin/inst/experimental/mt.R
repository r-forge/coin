
library("coin")
### get rid of namespace
load(file.path(.find.package("coin"), "R", "all.rda"))
source("bootstrap.R")

library("multtest")
set.seed(290875)

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
mydf <- dgp(n = 40, mu = 0.5)
mydf$X <- round(mydf$X, 3) * 1000
it <- independence_test(X ~ group, data = mydf, distribution = exact())

### generate scaled and shifted bootstrap distribution
a <- mboot(it, B = 49999)
Z <- drop(a)

### pvalues
pval(a, statistic(it, "linear"))
pval(a, statistic(it, "linear"), type = "single-step")
pvalue(it)
t.test(X ~ group, data = mydf, var.equal = TRUE)$p.value
MTP(t(as.matrix(mydf[["X"]])), Y = mydf[["group"]],
    test = "t.twosamp.equalvar", B = 9999)@adjp

t.test(X ~ group, data = mydf, var.equal = TRUE, alt = "less")$p.value
pval(a, statistic(it, "linear"), alt = "less")
t.test(X ~ group, data = mydf, var.equal = TRUE, alt = "greater")$p.value
pval(a, statistic(it, "linear"), alt = "greater")

mt.maxT(t(as.matrix(mydf[["X"]])), classlabel = as.numeric(mydf[["group"]]) -1)


### quantile of 
q <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)

### exact conditional distribution
qT <- sapply(q, function(i) qperm(it, i) * sqrt(variance(it)) + expectation(it))

### bootstrap distribution
qZ <- quantile(Z, q)

### quantile bootstrap distribution very close to those of exact distr.
cbind(qZ, qT)


mtest <- function(n = 100, mu = c(0, 2), B = 1000) {

    mydf <- dgp(n = n, mu = mu)
    xnames <- names(mydf)[names(mydf) != "group"]
    ip <- new("IndependenceProblem", x = mydf[xnames], y = mydf["group"])
    it <- independence_test(ip)
    Tobs <- drop(statistic(it, "linear"))
    Z <- mboot(it, B = B)
    p <- pval(Z, Tobs, type = "step-down")
    names(p) <- xnames
    p
}

Rprof("mtest")
a <- mtest(B = 10000)
Rprof(NULL)

### small simulation: X1 and X10 under H0, the others under `increasing'
### alternative
p <- c()
for (i in 1:1000) {
    print(i)
    ##p <- rbind(p, mtest(mu = c(seq(from = 0, to = 0.8, by = 0.1), 0)))
    p <- rbind(p, mtest(mu = rep(0, 5)))
}
save(p, file = "p.rda")
colMeans(p <= 0.05)
mean(apply(p, 1, min) <= 0.05)

### comparison with multtest
B <- 1000
mydf <- dgp(n = 100, mu = c(seq(from = 0, to = 0.8, by = 0.1), 0))
xnames <- names(mydf)[names(mydf) != "group"]
fm <- as.formula(paste(paste(xnames, collapse = "+"), "~ group"))
it <- independence_test(fm, data = mydf)
Tobs <- drop(statistic(it, "linear"))
Z <- mboot(it, B = B)
pval(Z, Tobs)

MTP(t(as.matrix(mydf[,xnames])), Y = mydf[["group"]], 
    test = "t.twosamp.equalvar")@adjp
