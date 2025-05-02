
ppplot <- function(x, ...)
    UseMethod("ppplot")

ppplot.formula <- function(formula, data, subset, na.action = na.pass, ...)
{
    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    if (length(attr(terms(formula[-2L]), "term.labels")) != 1L)
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    y <- mf[[response]]
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
    ## Call the default method.
    DATA <- split(mf[[response]], g)
    RVAL <- ppplot(x = DATA[[1L]], y = DATA[[2L]], xlab = levels(g)[1], ylab = levels(g)[2], ...)
    return(invisible(RVAL))
}

ppplot.default <- function(x, y, plot.it = TRUE,
            xlab = deparse1(substitute(x)),
            ylab = deparse1(substitute(y)), ...,
            conf.level = NULL, conf.args = list(col = NA, border = NULL)) {

    force(xlab)
    force(ylab)
    if (xlab == ylab) {
        xlab <- paste0("x = ", xlab)
        ylab <- paste0("y = ", ylab)
    }

    ex <- ecdf(x)
    sy <- sort(unique(y))
    py <- ecdf(y)(sy)
    px <- ex(sy)
    ret <- list(x = px, y = py)
    if (!plot.it)
        return(ret)

    plot(px, py, xlim = c(0, 1), ylim = c(0, 1), 
         xlab = xlab, ylab = ylab, type = "n", ...)

    if (!is.null(conf.level)) {
        prb <- seq_len(1000) / 1001
        res <- c(x, y)
        grp <- gl(2, 1, labels = c(xlab, ylab))
        grp <- grp[rep(1:2, c(length(x), length(y)))]
        args <- conf.args
        args$y <- res
        args$x <- grp
        args$conf.level <- conf.level
        args$border <- args$col <- NULL
        tt <- do.call("trafo.test", args)

        lwr <- tt$link$p(tt$link$q(prb) - tt$conf.int[1])
        upr <- tt$link$p(tt$link$q(prb) - tt$conf.int[2])
        x <- c(prb, rev(prb))
        y <- c(lwr, rev(upr))
        xn <- c(x[1L], rep(x[-1L], each = 2))
        yn <- c(rep(y[-length(y)], each = 2), y[length(y)])
        polygon(x = xn, y = yn, col = conf.args$col, border = conf.args$border)
        lines(prb, tt$link$p(tt$link$q(prb) - tt$estimate))
    }
    points(px, py, ...)
    return(invisible(ret)) 
}
