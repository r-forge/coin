
confband <- function(object, level = .95, ...)
    UseMethod("confband")

confband.ks.test <- function(object, level = .95, ...) {

    x <- object$data$x
    y <- object$data$y
    n.x <- length(x)
    n.y <- length(y)
    N <- n.x + n.y
    TIES <- length(unique(obs <- sort(c(x, y)))) != N

    ### <FIXME>: handle simulate.p.value = TRUE and also B </FIXME>
    if (TIES) {
        ca <- qsmirnov(level, m = n.x, n = n.y, 
                       z = obs, exact = object$exact)
    } else {
        ca <- qsmirnov(level, m = n.x, n = n.y, 
                       exact = object$exact)
    }

    if (object$alternative != "two-sided") 
        warning("computing two-sided confidance band from one-sided test")

    i <- as.double(seq_along(x))
    Re <- ceiling(i * N / n.x) - i
    Rl <- ceiling(i * N / n.x - ca * as.double(n.y)) - i
    Rl[Rl < 1] <- NA
    Rl[Rl > n.y] <- NA
    Rr <- floor(i * N / n.x - n.y / n.x + ca * as.double(n.y)) - i + 1
    Rr[Rr < 1] <- NA
    Rr[Rr > n.y] <- NA

    sx <- sort(x)
    sy <- sort(y)
    ret <- cbind(x = sx, Estimate = sy[Re], 
                 lwr = sy[Rl], upr = sy[Rr])
    if (!is.null(object$response) && !is.null(object$groups)) {
        attr(ret, "xlab") <- paste0(object$response, 
                                    " (", object$groups[1L], ")")
        attr(ret, "ylab") <- paste0(object$response, 
                                    " (", object$groups[2L], ")")
    }
    class(ret) <- c("confband.ks.test", class(ret))
    ret
}

plot.confband.ks.test <- function(x, y = NULL, 
    xlim = range(x[, "x"]), ylim = range(x[, c("lwr", "upr")], na.rm = TRUE), 
    ...) {

    xlab <- list(...)$xlab
    ylab <- list(...)$ylab
    if (is.null(xlab)) xlab <- attr(x, "xlab")
    if (is.null(ylab)) ylab <- attr(x, "ylab")
    if (is.null(xlab)) xlab <- "x"
    if (is.null(ylab)) ylab <- "y"

    plot(x[, "x"], x[, "Estimate"], xlim = xlim, ylim = ylim,
         type = "S", xlab = xlab, ylab = ylab, ...)
    abline(a = 0, b = 1, col = "lightgrey")
    lines(x[, "x"], x[, "lwr"], type = "S", lty = 2)
    lines(x[, "x"], x[, "upr"], type = "S", lty = 2)
    invisible(x)
}

