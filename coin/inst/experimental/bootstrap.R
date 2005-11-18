

### compute expectation and variance for observations with x != 0
expectvarcontrast <- function(x, y, w) {

    if (max(abs(colSums(x))) > sqrt(.Machine$double.eps))
        stop(sQuote("x"), " does not define contrasts")

    ev <- apply(x, 2, function(xx) {
        coin:::expectvaronly(matrix(xx, ncol = 1), y, w * as.numeric(abs(xx) > 0))
    })

    list(E = sapply(ev, function(e) e$E), V = sapply(ev, function(e) e$V))
}


### bootstrap scaled and shifted test statistics
mboot <- function(it, B = 1000, xcontrast = TRUE) {

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
    if (!xcontrast) {
        lambda0 <- expectation(it)
        tau0 <- variance(it)
    } else {
        ### in case x defines contrasts (i.e. colSums(x) == 0) 
        ### compute expectation and variance for x[x != 0] only
        ### Otherwise, subset pivotality fails (I guess).
        ev <- expectvarcontrast(xtrans, ytrans, weights)
        lambda0 <- ev$E
        tau0 <- ev$V
    }

    ### for each bootstrap sample, calculate _multivariate_ linear statistic
    T <- matrix(0, nrow = length(statistic(it, "linear")), ncol = ncol(bs))
    for (i in 1:ncol(bs))
        ### multivariate linear statistic
        T[,i] <- .Call("R_LinearStatistic", xtrans, ytrans, bs[,i], 
                       PACKAGE = "coin")

    ### compute Z matrix: bootstrap scaled and shifted test statistics
    ET <- rowMeans(T)
    VT <- apply(T, 1, var)
    fact <- sqrt(pmin(1, tau0 / VT))
    ### typo on page 42: Z = nu_0(T - ET) + lambda_0 is right
    Z <- fact * (T - ET) + lambda0
    # rownames(RET) <- rownames(statistic(it, "linear"))
    ### this is a (M x B) matrix of bootstrap T statistics
    return(Z)
}

cmax <- function(x) {
    do.call(pmax, as.data.frame(t(x)))
}

### single-step and step-down maxT
pval <- function(Z, Tobs, type = c("single-step", "step-down"), 
                 alternative = c("two.sided", "less", "greater")) {

    alternative <- match.arg(alternative)
    type <- match.arg(type)
    Tobs <- drop(Tobs)

    ### standardized Z and Tobs with empirical mean and variance
    EZ <- rowMeans(Z)
    sdZ <- apply(Z, 1, sd)
    Z <- (Z - EZ)/sdZ
    Tobs <- (Tobs - EZ)/sdZ

    if (alternative == "two.sided") {
        Z <- abs(Z)
        Tobs <- abs(Tobs)
    }
    if (alternative == "less") {
        Z <- -Z
        Tobs <- -Tobs
    }

    if (type == "single-step")
        return(sapply(Tobs, function(x) mean(cmax(Z) > x)))
    if (type == "step-down")
        return(drop(coin:::sdmaxT(t(Z), matrix(Tobs))))
}

mt.bootstrap <- function(object, B = 1000, 
                      type = c("single-step", "step-down"), ...) {

    Tobs <- statistic(object, "linear")
    Z <- mboot(object, B = B, ...)
    pval(Z, Tobs, alternative = object@statistic@alternative)
}

