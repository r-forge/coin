
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
                    W <- sum(1:nrow(xt) * xt[,1])
                    n.x <- unname(cs)[1]
                    n.y <- unname(cs)[2]
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
        U <- sapply(rt, function(x) sum(x[,1] * res))
        STATISTIC <- c("Score Z" = sum(res * xt[,1]))
        PVAL <- switch(alternative, two.sided = {
                    mean(abs(U) >= abs(STATISTIC))
                }, greater = mean(U >= STATISTIC),
                   less = mean(U <= STATISTIC))
        TYPE <- "approximate"
    } else {
        EV <- SW(res, xt)
        STATISTIC <- c("Score Z" = sum(res * xt[,1]))
        if (alternative == "less") {
            PVAL <- pnorm(STATISTIC, mean = EV$E, sd = sqrt(EV$V))
        } else if (alternative == "greater") {
            PVAL <- pnorm(STATISTIC, mean = EV$E, sd = sqrt(EV$V), lower.tail = FALSE)
        } else {
            PVAL <- 2 * pnorm(-abs(STATISTIC), mean = EV$E, sd = sqrt(EV$V))
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
                    qU <- qwilcox(p, m = cs[1], n = cs[2])
                    n.x <- cs[1]
                    n.y <- cs[2]
                    return(-((qU + n.x * (n.x + 1) / 2) * cf[2] + cf[1] * cs[1]))
                }
            }
        }
        rt <- r2dtable(B, r = rowSums(xt), c = cs)
        U <- sapply(rt, function(x) sum(x[,1] * res))
        return(quantile(U, probs = p))
    } else {
        EV <- SW(res, xt)
        return(qnorm(p, mean = EV$E, sd = sqrt(EV$V)))
    }
}

SW <- function(res, xt) {

    x <- gl(2, nrow(xt))
    r <- rep(res, 2)
    w <- c(xt)

    x <- (x == levels(x)[2])
    wx <- x * w
    sw <- sum(w)
    Eh <- sum(r * w) / sw
    E <- sum(wx) * Eh
    Vh <- sum(w * (r - Eh)^2) / sw
    V <- sum(w) / (sum(w) - 1) * Vh * sum(wx^2) - 1 / (sum(w) - 1) * Vh * sum(wx)^2
    list(E = E, V = V)
}
