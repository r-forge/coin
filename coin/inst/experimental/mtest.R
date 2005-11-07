
library("coin")

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

    ### for each bootstrap sample, calculate _multivariate_ linear statistic
    ### and standardize
    cEV <- coin:::expectvaronly
    RET <- apply(bs, 2, function(b) {
        ### multivariate linear statistic
        ls <- .Call("R_LinearStatistic", xtrans, ytrans, b, PACKAGE = "coin")
        ### conditional expectation and variance
        ev <- cEV(xtrans, ytrans, b)
        ### standardization
        ss <- (ls - ev$E) / sqrt(ev$V)
        return(ss)
    })
    # rownames(RET) <- rownames(statistic(it, "linear"))
    ### this is a (M x B) matrix of bootstrap T statistics
    return(RET)
}

### example with y1 under H_0 and y2 under H_1
mydf <- data.frame(y1 = rnorm(50), 
                   y2 = rnorm(50) + rep(c(1, 0), c(25, 25)), 
                   x = gl(2, 25))

### independence test with bivariate response in two groups
it <- independence_test(y1 + y2 ~ x, data = mydf)

### observed standardized statistics
Tobs <- statistic(it, "standardized")

### B bootstrap samples of Tobs
Tboot <- mboot(it)

### bootstrap shifted test statistics
Z <- Tboot - rowMeans(Tboot)

### bootstrap unadjusted p-value
rowMeans(Z > as.vector(Tobs))
