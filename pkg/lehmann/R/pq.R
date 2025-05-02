
statpvalPerm <- function(res, xt, B = 0, alternative) {

    rs <- rowSums(xt)
    cs <- colSums(xt)

    if (B) {
        if (!is.finite(B)) {
            B <- 10000
            if (!isTRUE(all.equal(unique(rs), 1))) {
                warning("cannot compute exact distribution")
            } else {
                i <- cbind(1, 1:nrow(xt))
                lr <- lm.wfit(x = i, y = res, w = rep(1, length(res)))
                if (sum(lr$residuals) > .Machine$double.eps) {
                    warning("cannot compute exact distribution")
                } else {
                    cf <- lr$coefficients
                    W <- sum(1:nrow(xt) * xt[,2])
                    n.x <- unname(cs)[2]
                    n.y <- unname(cs)[1]
                    STATISTIC <- c("U" = W - n.x * (n.x + 1) / 2)
                    PVAL <- switch(alternative, two.sided = {
                        p <- if (STATISTIC > (n.x * n.y/2)) 
                            pwilcox(STATISTIC - 1, n.x, n.y, lower.tail = FALSE) 
                        else 
                            pwilcox(STATISTIC, n.x, n.y)
                        min(2 * p, 1)
                    }, greater = {
                        pwilcox(STATISTIC - 1, n.x, n.y, lower.tail = FALSE)
                    }, less = {
                        pwilcox(STATISTIC, n.x, n.y)
                    })
                    return(list(STATISTIC = STATISTIC, PVAL = PVAL, TYPE = "exact (Wilcoxon)"))
                }
            }
        }

        rt <- r2dtable(B, r = rs, c = cs)
        U <- sapply(rt, function(x) sum(x[,2] * res))
        STATISTIC <- c("Score Z" = sum(res * xt[,2]))
        PVAL <- switch(alternative, two.sided = {
                    mean(abs(U) >= abs(STATISTIC))
                }, greater = mean(U >= STATISTIC),
                   less = mean(U <= STATISTIC))
        TYPE <- "approximate"
    } else {
        EV <- SW(res, xt)
        STATISTIC <- c("Score Z" = sum(xt[,2] * res))
        if (alternative == "less") {
            PVAL <- pnorm(STATISTIC, mean = EV$Expectation, sd = sqrt(EV$Covariance))
        } else if (alternative == "greater") {
            PVAL <- pnorm(STATISTIC, mean = EV$Expectation, sd = sqrt(EV$Covariance), lower.tail = FALSE)
        } else {
            PVAL <- 2 * pnorm(-abs(STATISTIC), mean = EV$Expectation, sd = sqrt(EV$Covariance))
        }
        TYPE <- "asymptotic"
    }
    return(list(STATISTIC = STATISTIC, PVAL = PVAL, TYPE = TYPE))
}

qPerm <- function(p, res, xt, B = 0) {
    rs <- rowSums(xt)
    cs <- colSums(xt)
    if (B) {
        if (!is.finite(B)) {
            B <- 10000
            if (!isTRUE(all.equal(unique(rs), 1))) {
                warning("cannot compute exact distribution")
            } else {
                i <- cbind(1, 1:nrow(xt))
                lr <- lm.wfit(x = i, y = res, w = rep(1, length(res)))
                if (sum(lr$residuals) > .Machine$double.eps) {
                    warning("cannot compute exact distribution")
                } else {
                    cf <- lr$coefficients
                    n.x <- unname(cs)[2]
                    n.y <- unname(cs)[1]
                    qU <- qwilcox(p, m = n.x, n = n.y)
                    return(-((qU + n.x * (n.x + 1) / 2) * cf[2] + cf[1] * cs[1]))
                }
            }
        }
        rt <- r2dtable(B, r = rowSums(xt), c = cs)
        U <- sapply(rt, function(x) sum(x[,2] * res))
        return(quantile(U, probs = p))
    } else {
        EV <- SW(res, xt)
        return(qnorm(p, mean = EV$Expectation, sd = sqrt(EV$Covariance)))
    }
}

SW <- function(res, xt) {

    Y <- matrix(res, ncol = 1, nrow = 2 * length(res))
    weights <- c(xt)
    x <- gl(2, nrow(xt))
    X <- matrix(x == levels(x)[2], ncol = 1)

    w. <- sum(weights)
    wX <- weights * X
    wY <- weights * Y
    ExpX <- colSums(wX)
    ExpY <- colSums(wY) / w.
    CovX <- crossprod(X, wX)
    Yc <- t(t(Y) - ExpY)
    CovY <- crossprod(Yc, weights * Yc) / w.
    T <- crossprod(X, wY)
    Exp <- kronecker(ExpY, ExpX)
    Cov <- w. / (w. - 1) * kronecker(CovY, CovX) -
           1 / (w. - 1) * kronecker(CovY, tcrossprod(ExpX))
    list(Expectation = as.vector(Exp),
         Covariance = Cov)
}
