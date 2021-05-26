
psmirnov <- function(q, n.x, n.y, obs = NULL, two.sided = TRUE,
                     lower.tail = TRUE, log.p = FALSE) {

    ###
    ### Distribution function Prob(D < q) for the Smirnov test statistic
    ###
    ### D = sup_c | ECDF_x(c) - ECDF_y(c) | 	(two.sided)
    ###
    ### D = sup_c ( ECDF_x(c) - ECDF_y(c) ) 	(!two.sided)
    ###

    if (length(q) > 1)
        return(sapply(q, psmirnov, n.x = n.x, n.y = n.y, obs = obs, 
                                   lower.tail = lower.tail, log.p = log.p))

    n.x <- as.integer(n.x)
    n.y <- as.integer(n.y)
    N <- n.x + n.y
    if (length(n.x) > 1 || length(n.y) > 1 || n.x < 2 || n.x < 2)
        stop(sQuote("n.x"), "and/or", sQuote("n.y"), "misspecified")
    umat <- matrix(0, nrow = n.x + 1L, ncol = n.y + 1L)
    ix <- 0:n.x
    jy <- 0:n.y

    ### see stats/src/ks.c line 103ff
    q <- (0.5 + floor(as.double(q) * n.x * n.y - 1e-7)) / (n.x * n.y);

    abs2 <- if (two.sided) abs else function(x) x

    if (!is.null(obs)) {
        if (length(obs) != N)
            stop(sQuote("length(obs)"), "is not equal to", sQuote("n.x + n.y"))
        z <- sort(obs)
        eqz <- diff(z) < sqrt(.Machine$double.eps)
        for (i in ix) {
            for (j in jy) {
                cmat <- 0L
                if (abs2(i / n.x - j / n.y) < q) 
                    cmat <- 1L
                k <- i + j
                if (k > 0 && k < N && eqz[k])
                    cmat <- 1L
                fct <- 1L
                if (i > 0 && j > 0) 
                    fct <- umat[i + 1, j] + umat[i, j + 1]
                umat[i + 1, j + 1] <- cmat * fct
            }
        }            
    } else {
        for (i in ix) {
            for (j in jy) {
                cmat <- 0L
                if (abs2(i / n.x - j / n.y) < q) 
                    cmat <- 1L
                fct <- 1    
                if (i > 0 && j > 0) 
                    fct <- umat[i + 1, j] + umat[i, j + 1]
                umat[i + 1, j + 1] <- cmat * fct
            }
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

### Schröer & Trenkler (1994)
1 - psmirnov(3 / 7, 5, 7)
1 - psmirnov(3 / 7, 5, 7, two.sided = FALSE)

1 - psmirnov(3 / 7, 5, 7, obs = 1:12)
1 - psmirnov(3 / 7, 5, 7, obs = 1:12, two.sided = FALSE)

#system.time(psmirnov(.5, 100, 200))

### Schröer & Trenkler
1 - psmirnov(3 / 7, 5, 7, obs = c(1, 2, 2, 3, 3, 1, 2, 3, 3, 4, 5, 6))

x <- c(1, 2, 2, 3, 3)
y <- c(1, 2, 3, 3, 4, 5, 6)

s <- replicate(10000, {z <- sample(c(x, y)); ks.test(z[1:5], z[-(1:5)])$stat})

mean(s >= 3/7)

### Wilcox (1997)
x <- c(0, 32, 9, 0, 2, 0, 41, 0, 0, 0, 6, 18, 3, 3, 0, 11, 11, 2, 0, 11)
y <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 8, 0, 3, 0, 0, 32, 12, 2, 0, 0, 0)

D <- ks.test(x, y, exact = FALSE)$stat

1 - psmirnov(D + .01, n.x = length(x), n.y = length(y), obs = c(x, y))
