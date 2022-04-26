
qqplot <- function(x, y, plot.it = TRUE,
                   xlab = deparse1(substitute(x)),
                   ylab = deparse1(substitute(y)), ..., 
                   conf.level = NULL, exact = NULL, 
                   simulate = FALSE, B = 2000)
{
    sx <- sort(x)
    sy <- sort(y)
    lenx <- length(sx)
    leny <- length(sy)
    if( leny < lenx )
        sx <- approx(1L:lenx, sx, n = leny)$y
    if( leny > lenx )
        sy <- approx(1L:leny, sy, n = lenx)$y

    if (is.null(conf.level)) {
        if(plot.it)
            plot(sx, sy, xlab = xlab, ylab = ylab, ...)
        return(invisible(list(x = sx, y = sy)))
    }

    if (conf.level < 0 || conf.level > 1)
            stop("'conf.level' is not a probability")

    N <- lenx + leny
    if (is.null(exact)) exact <- lenx * leny < 10000
    z <- c(x, y)
    TIES <- any(duplicated(z))
    if (!TIES) z <- NULL
    ca <- qsmirnov(conf.level, sizes = c(lenx, leny), z = z,
                   exact = exact && !simulate, 
                   simulate = simulate, B = B)

    ### Switzer (1976, Biometrika) 10.1093/biomet/63.1.13
    i <- as.double(seq_along(x))
    Re <- ceiling(i * N / lenx) - i
    Rl <- ceiling(i * N / lenx - ca * as.double(leny)) - i
    Rl[Rl < 1] <- NA
    Rl[Rl > leny] <- NA
    Rr <- floor(i * N / lenx - leny / lenx + ca * as.double(leny)) - i + 1
    Rr[Rr < 1] <- NA
    Rr[Rr > leny] <- NA
    lwr <- sy[Rl]
    upr <- sy[Rr]

    if (plot.it) {
        plot(sx, sy, xlab = xlab, ylab = ylab, type = "n", ...)
        points(sx, sy, ...)
        lines(sx, lwr, type = "S", col = "grey", ...)
        lines(sx, upr, type = "S", col = "grey", ...)
    }	

    return(invisible(list(x = sx, y = sy, lwr = lwr, upr = upr)))
}

x <- rnorm(50)
y <- rnorm(50, sd = .95)
ex <- TRUE
ks <- ks.test(x, y, exact = ex)
pval <- ks$p.value
qqplot(x, y, conf.level = 1 - pval, exact = ex)
abline(a = 0, b = 1)

qsmirnov(1 - pval, sizes = c(length(x), length(y)), exact = ex)
ks$stat
