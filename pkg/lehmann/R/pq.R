
statpvalPerm <- function(r, xt, B = 0, alternative) {
    x <- gl(2, length(r))
    r <- rep(r, 2)
    w <- c(xt)
    distribution <- "asymptotic"
    if (B) {
        if (is.finite(B))
            distribution <- approximate(nresample = B)
        else
            distribution <- "exact"
    }
    d <- data.frame(r = r, x = x, w = w)    
    it <- independence_test(r ~ x, data = d, weights = ~ w, 
                            teststat = "scalar",
                            alternative = alternative,
                            distribution = distribution)
    c(Z = statistic(it, "standardized"), pval = pvalue(it))
}

qPerm <- function(p, r, xt, B = 0) {
    x <- gl(2, nrow(xt))
    r <- rep(r, 2)
    w <- c(xt)
    distribution <- "asymptotic"
    if (B) {
        if (is.finite(B))
            distribution <- approximate(nresample = B)
        else
            distribution <- "exact"
    }
        
    d <- data.frame(r = r, x = x, w = w)    
    it <- independence_test(r ~ x, data = d, weights = ~ w,
                            teststat = "scalar",
                            distribution = distribution)
    qperm(it, p = p)
}

