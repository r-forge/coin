
psmirnov <- function(q, n.x, n.y, obs = 1:(n.x + n.y), 
                     lower.tail = TRUE, log.p = FALSE) {

    ### see stats/src/ks.c line 103ff
    q <- (0.5 + floor(q * n.x * n.y - 1e-7)) / (n.x * n.y);

    z <- sort(obs)
    cmat <- umat <- matrix(0, nrow = n.x + 1, ncol = n.y + 1)

    ix <- 0:n.x
    jy <- 0:n.y

    for (i in ix) {
        for (j in jy) {
            if (abs(i / n.x - j / n.y) < q) 
                cmat[i + 1, j + 1] <- 1
            k <- i + j
            if (k > 0 && k < length(z) && z[k] == z[k + 1])
                cmat[i + 1, j + 1] <- 1
        }
    }

    for (i in ix) {
        for (j in jy) {
            if (i * j == 0) 
                fct <- 1    
            else
                fct <- umat[i + 1, j] + umat[i, j + 1]
            umat[i + 1, j + 1] <- cmat[i + 1, j + 1] * fct
        }
    }            

    term <- lgamma(n.x + n.y + 1) - lgamma(n.x + 1) - lgamma(n.y + 1)
    ret <- umat[n.x + 1, n.y + 1]
    if (lower.tail && log.p)
        return(log(ret) - term)
    if (lower.tail && !log.p)
        return(ret / exp(term))
    if (!lower.tail && !log.p)
        return(1 - ret / exp(term))
    return(log1p(-ret / exp(term)))
}
