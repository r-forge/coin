
library("coin")
library("multtest")
set.seed(290875)

mboot <- function(it, B = 1000) {

    ### get transformation g(X) and influence function h(Y)
    xtrans <- coin:::get_xtrans(it)
    ytrans <- coin:::get_ytrans(it)
    ### weights
    weights <- coin:::get_weights(it)

    ### generate B bootstrap samples (via weights from multinomial
    ### distribution)
    n <- sum(weights)
    bs <- rmultinom(B, n, weights / n)
    storage.mode(bs) <- "double"

    ### null values (\lambda_0 and \tau_0)
    lambda0 <- expectation(it)
    tau0 <- variance(it)

    ### for each bootstrap sample, calculate _multivariate_ linear statistic
    bootT <- vector(mode = "list", length = ncol(bs))
    for (i in 1:ncol(bs))
        ### multivariate linear statistic
        bootT[[i]] <- .Call("R_LinearStatistic", xtrans, ytrans, bs[,i], PACKAGE = "coin")

    ### compute Z matrix
    T <- matrix(unlist(bootT), ncol = length(bootT))
    ### bootstrap shifted test statistics
    ET <- rowMeans(T)
    VT <- apply(T, 1, var)
    fact <- sqrt(pmin(1, tau0 / VT))
    ### typo on page 42: Z = nu_0(T - ET) + lambda_0 is right
    Z <- fact * (T - ET) + lambda0
    # rownames(RET) <- rownames(statistic(it, "linear"))
    ### this is a (M x B) matrix of bootstrap T statistics
    return(Z)
}

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

cmax <- function(x) {
    do.call(pmax, as.data.frame(t(x)))
}


pval <- function(Z, Tobs, alternative = c("two.sided", "less", "greater")) {

    alternative <- match.arg(alternative)
    Tobs <- drop(Tobs)
    EZ <- rowMeans(Z)
    Z <- Z - EZ
    Tobs <- Tobs - EZ
    switch(alternative, "two.sided" = {
        maxabsZ <- cmax(abs(Z))
        sapply(abs(Tobs), function(x) mean(maxabsZ > x))
    }, "less" = {
        rowMeans(cmax(Z) > max(Tobs))
        maxZ <- cmax(Z)
        sapply(abs(Tobs), function(x) mean(maxZ > x))
    }, "greater" = {
        maxZ <- cmax(Z)
        sapply(abs(Tobs), function(x) mean(maxZ < x))
    })
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
pvalue(it)
t.test(X ~ group, data = mydf, var.equal = TRUE)$p.value
MTP(t(as.matrix(mydf[["X"]])), Y = mydf[["group"]],
    test = "t.twosamp.equalvar", B = 9999)@adjp

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
    fm <- as.formula(paste(paste(xnames, collapse = "+"), "~ group"))
    it <- independence_test(fm, data = mydf)
    Tobs <- drop(statistic(it, "linear"))
    Z <- mboot(it, B = B)
    p <- pval(Z, Tobs)
    names(p) <- xnames
    p
}

### small simulation: X1 and X10 under H0, the others under `increasing'
### alternative
p <- c()
for (i in 1:1000) {
    print(i)
    p <- rbind(p, mtest(mu = c(seq(from = 0, to = 0.8, by = 0.1), 0)))
}
colMeans(p <= 0.05)
mean(pmin(p[,1], p[,10]) <= 0.05)

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
