
.rwilcox <- function(n, sizes, z = NULL) 
{

    if (n < 0) 
        stop("invalid arguments")
    if (n == 0L) 
        return(numeric(0))
    n <- floor(n)
    if (length(sizes) != 2L) 
        stop("argument 'sizes' must be a vector of length 2")
    n.x <- sizes[1L]
    n.y <- sizes[2L]
    if (n.x < 1) 
        stop("not enough 'x' data")
    if (n.y < 1) 
        stop("not enough 'y' data")
    n.x <- floor(n.x)
    n.y <- floor(n.y)
    if (is.null(z))
        return(rwilcox(n, n.x, n.y))

    ### sample from conditional distribution
    uz <- sort(unique(z))
    tb <- r2dtable(n, table(z), sizes)
    ret <- colSums(uz * matrix(do.call("rbind", tb)[,1L], nrow = length(uz)))
    ret <- ret - n.x * (n.x + 1) / 2
    ret
}
    
.dwilcox <- function(x, sizes, z = NULL, log = FALSE) 
{

    if (length(sizes) != 2L) 
        stop("argument 'sizes' must be a vector of length 2")
    n.x <- sizes[1L]
    n.y <- sizes[2L]
    if (n.x < 1) 
        stop("not enough 'x' data")
    if (n.y < 1) 
        stop("not enough 'y' data")
    n.x <- floor(n.x)
    n.y <- floor(n.y)
    N <- n.x + n.y
    if (is.null(z))
        return(dwilcox(x, n.x, n.y, log))


    ### x is U statistic
    if (is.numeric(x)) 
        q <- as.double(x)
    else stop("argument 'x' must be numeric")
    ret <- rep.int(0, length(q))
    ret[is.na(q)] <- NA
    IND <- which(!is.na(ret))
    if (!length(IND)) 
        return(ret)

    if (length(z) != N)
        stop(sQuote("length(z)"), "is not equal to", sQuote("sum(sizes)"))
    ### ranks x.5 are possible, algo expects integer scores
    fct <- ifelse(isTRUE(all.equal(max(z - floor(z)), 0)), 1, 2)
    sc <- as.integer(sort(floor(z * fct)))

    dens <- .Call("C_StreitbergRoehmel",
                  score = sc,
                  m = as.integer(n.x))
    ### W := sum(ranks[1:n.x])
    W <- seq_along(dens)
    ### q is U
    q <- q + n.x * (n.x + 1) / 2
    q <- q[IND] * fct
    W <- W[match(q, W)]
    ret[IND] <- ifelse(is.na(W), 0, dens[W])
    if (log) return(log(ret))
    return(ret)
}

.pwilcox <- function(q, sizes, z = NULL, exact = TRUE,
                     conditional = exact, simulate = FALSE, B = 2000, 
                     lower.tail = TRUE, log.p = FALSE) 
{

    ### q is U
    if (is.numeric(q)) 
        q <- as.double(q)
    else stop("argument 'q' must be numeric")
    ret <- rep.int(0, length(q))
    ret[is.na(q)] <- NA
    IND <- which(!is.na(ret))
    if (!length(IND)) 
        return(ret)
    if (length(sizes) != 2L) 
        stop("argument 'sizes' must be a vector of length 2")
    n.x <- sizes[1L]
    n.y <- sizes[2L]
    if (n.x < 1) 
        stop("not enough 'x' data")
    if (n.y < 1) 
        stop("not enough 'y' data")
    n.x <- floor(n.x)
    n.y <- floor(n.y)
    N <- n.x + n.y
    exact <- exact && !simulate

    if (!exact) {
        if (simulate) {
            Usim <- .rwilcox(B, sizes = sizes, z = z)
            ### need P(U <= q)
            ret[IND] <- sapply(q[IND], function(x) sum(Usim < x + sqrt(.Machine$double.eps))) 
            ret[IND] <- ret[IND] / (B + 1L)
            if (log.p & lower.tail)
                return(log(ret))
            if (!log.p & lower.tail)
                return(ret)
            if (log.p & !lower.tail)
                return(log1p(-ret))
            if (!log.p & !lower.tail)
                return(1 - ret)
        }

        ### mean and variance of U statistic
        if (is.null(z)) {
            ### w/o ties, conditional == unconditional
            E <- n.x * n.y / 2
            V <- n.x * n.x / 12 * (n.x + n.y + 1)
        } else {
            if (conditional) {
                ### mean and variance of conditional permutation
                ### distribution
                ### Strasser, H. and Weber, C.  (1999).  On the asymptotic 
                ### theory of permutation statistics.  _Mathematical Methods 
                ### of Statistics_ *8*(2), 220-250.
                E <- n.x * mean(z) - n.x * (n.x + 1) / 2
                Vh <- mean((z - mean(z))^2)
                V <- N / (N - 1) * Vh * n.x - 1 / (N - 1) * Vh * n.x^2
            } else {
                ### mean and variance of unconditional distribution
                NTIES <- table(z)
                E <- n.x * n.y / 2
                V <- n.x * n.y / 12 * 
                      ((N + 1) - sum(NTIES^3 - NTIES) / 
                      (N * (N - 1)))
            }
        }
        return(pnorm(q, mean = E, sd = sqrt(V), 
                     lower.tail = lower.tail, log.p = log.p))
    } 

    ### exact conditional distribution
    ### support of U
    U <- 0:(2 * n.x * n.y) / 2L
    ### density
    d <- .dwilcox(U, sizes = sizes, z = z)
    ### cdf
    ret[IND] <- sapply(q[IND], function(x) sum(d[U < x + sqrt(.Machine$double.eps)]))
    if (log.p & lower.tail)
        return(log(ret))
    if (!log.p & lower.tail)
        return(ret)
    if (log.p & !lower.tail)
        return(log1p(-ret))
    if (!log.p & !lower.tail)
        return(1 - ret)
}

.qwilcox <- function(p, sizes, z = NULL, exact = TRUE,
                     conditional = exact, simulate = FALSE, B = 2000,
                     lower.tail = TRUE, log.p = FALSE) 
{

    if (length(sizes) != 2L) 
        stop("argument 'sizes' must be a vector of length 2")
    n.x <- sizes[1L]
    n.y <- sizes[2L]
    if (n.x < 1) 
        stop("not enough 'x' data")
    if (n.y < 1) 
        stop("not enough 'y' data")
    n.x <- floor(n.x)
    n.y <- floor(n.y)
    N <- n.x + n.y

    if (is.null(z))
        return(qwilcox(p, n.x, n.y, lower.tail, log.p))

    if (is.numeric(p)) 
        p <- as.double(p)
    else stop("argument 'p' must be numeric")
    ret <- rep.int(0, length(p))
    if (log.p) {
        ret[is.na(p) | p > 0] <- NA
    } else {
        ret[is.na(p) | p < 0 | p > 1] <- NA
    }
    IND <- which(!is.na(ret))
    if (!length(IND)) 
        return(ret)

    if (!exact) {
        if (simulate) {
            Usim <- .rwilcox(B, sizes = sizes, z = z)
            U <- 0:(2 * n.x * n.y) / 2L
            prb <- ecdf(Usim)(U)
            if (log.p) prb <- log(prb)
            if (!lower.tail) {
                if (log.p) p <- log1p(-exp(p))
                p <- 1 - p
            }
            ret[IND] <- sapply(p[IND], function(u) U[prb > u][1L])
            return(ret)
        }

        ### see pwilcox
        if (is.null(z)) {
            E <- n.x * n.y / 2
            V <- n.x * n.x / 12 * (n.x + n.y + 1)
        } else {
            if (conditional) {
                E <- n.x * mean(z) - n.x * (n.x + 1) / 2
                Vh <- mean((z - mean(z))^2)
                V <- N / (N - 1) * Vh * n.x - 1 / (N - 1) * Vh * n.x^2
            } else {
                NTIES <- table(z)
                E <- n.x * n.y / 2
                V <- n.x * n.y / 12 * 
                      ((N + 1) - sum(NTIES^3 - NTIES) / 
                      (N * (N - 1)))
            }
        }
        return(floor(qnorm(p, mean = E, sd = sqrt(V), 
                           lower.tail = lower.tail, log.p = log.p)))
    } 

    U <- 0:(2 * n.x * n.y) / 2L
    prb <- .pwilcox(U, sizes = sizes, z = z, log.p = log.p)
    if (!lower.tail) {
        if (log.p) p <- log1p(-exp(p))
        p <- 1 - p
    }

    ret[IND] <- sapply(p[IND], function(x) U[prb > x][1L])
    return(ret)
}
