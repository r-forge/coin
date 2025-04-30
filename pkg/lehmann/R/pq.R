
statpvalPerm <- function(r, x, w, B = 0, alternative) {
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

qPerm <- function(p, r, x, w, B = 0) {
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

