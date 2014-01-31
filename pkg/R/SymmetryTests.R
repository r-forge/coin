### Friedman Test
friedman_test <- function(object, ...) UseMethod("friedman_test")

friedman_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    ft("friedman_test", "SymmetryProblem", formula, data, subset,
       frame = parent.frame(), ...)
}

friedman_test.SymmetryProblem <- function(object,
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"),
                 " does not represent a K-sample problem",
                 " (maybe the grouping variable is not a factor?)")
        return(TRUE)
    }

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(distribution = distribution,
                       ytrafo = function(data)
                           trafo(data, numeric_trafo = rank_trafo,
                                 block = object@block),
                       check = check)
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <- if (is.ordered(object@x[[1]])) "scalar"
                     else "quad"

    RET <- do.call("symmetry_test", c(list(object = object), args))

    if (is_singly_ordered(RET@statistic))
        RET@method <- "Page Test"
    else
        RET@method <- "Friedman Test"

    return(RET)
}


### Wilcoxon signed-rank test
wilcoxsign_test <- function(object, ...) UseMethod("wilcoxsign_test")

wilcoxsign_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    if (is.null(d$bl))
        d <- list(y = data.frame(c(d$y[[1]], d$x[[1]])),
                  x = data.frame(gl(2, length(d$x[[1]]))),
                  block = factor(rep.int(1:length(d$x[[1]]), 2)))
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("wilcoxsign_test", c(list(object = sp), list(...)))
    return(RET)
}

wilcoxsign_test.SymmetryProblem <- function(object,
    zero.method = c("Pratt", "Wilcoxon"), ...) {

    zero.method <- match.arg(zero.method)

    y <- object@y[[1]]
    x <- object@x[[1]]
    block <- object@block

    if (!is.numeric(y))
        stop(sQuote("y"), " is not a numeric variable")
    if (is.factor(x)) {
        if (nlevels(x) != 2)
            stop(sQuote("x"), " is not a factor with two levels")
        diffs <- tapply(1:length(y), block, function(b)
            y[b][x[b] == levels(x)[1]] - y[b][x[b] == levels(x)[2]]
        )
    } else {
        stop(sQuote("x"), " is not a factor")
    }

    abs_diffs <- abs(diffs)
    if (all(abs_diffs < .Machine$double.eps))
        stop("all pairwise differences equal zero")

    pos_abs_diffs <- abs_diffs > 0
    if (zero.method == "Pratt") {
        rank_abs_diffs <- rank(abs_diffs)
        pos <- (rank_abs_diffs * (diffs > 0))[pos_abs_diffs]
        neg <- (rank_abs_diffs * (diffs < 0))[pos_abs_diffs]
    } else {
        diffs <- diffs[pos_abs_diffs]
        abs_diffs <- abs_diffs[pos_abs_diffs]
        rank_abs_diffs <- rank(abs_diffs)
        pos <- rank_abs_diffs * (diffs > 0)
        neg <- rank_abs_diffs * (diffs < 0)
    }
    n <- length(pos)

    y <- as.vector(rbind(pos, neg))
    x <- factor(rep.int(0:1, n), labels = c("pos", "neg"))
    block <- gl(n, 2)

    ip <- new("IndependenceProblem", x = data.frame(x = x),
              y = data.frame(y = y), block = block)

    args <- setup_args(teststat = "scalar", paired = TRUE)

    RET <- do.call("independence_test", c(list(object = ip), args))

    if (zero.method == "Pratt")
        RET@method <- "Wilcoxon-Pratt Signed-Rank Test"
    else
        RET@method <- "Wilcoxon Signed-Rank Test"
    RET@nullvalue <- 0

    return(RET)
}


### sign test
sign_test <- function(object, ...) UseMethod("sign_test")

sign_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    if (is.null(d$bl))
        d <- list(y = data.frame(c(d$y[[1]], d$x[[1]])),
                  x = data.frame(gl(2, length(d$x[[1]]))),
                  block = factor(rep.int(1:length(d$x[[1]]), 2)))
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("sign_test", c(list(object = sp), list(...)))
    return(RET)
}

sign_test.SymmetryProblem <- function(object, ...) {

    y <- object@y[[1]]
    x <- object@x[[1]]
    block <- object@block

    if (!is.numeric(y))
        stop(sQuote("y"), " is not a numeric variable")
    if (is.factor(x)) {
        if (nlevels(x) != 2)
            stop(sQuote("x"), " is not a factor with two levels")
        diffs <- tapply(1:length(y), block, function(b)
            y[b][x[b] == levels(x)[1]] - y[b][x[b] == levels(x)[2]]
        )
    } else {
        stop(sQuote("x"), " is not a factor")
    }

    abs_diffs <- abs(diffs)
    if (all(abs_diffs < .Machine$double.eps))
        stop("all pairwise differences equal zero")

    diffs <- diffs[abs_diffs > 0]
    n <- length(diffs)

    y <- as.vector(rbind(as.numeric(diffs > 0), as.numeric(diffs < 0)))
    x <- factor(rep.int(0:1, n), labels = c("pos", "neg"))
    block <- gl(n, 2)

    ip <- new("IndependenceProblem", x = data.frame(x = x),
              y = data.frame(y = y), block = block)

    args <- setup_args(teststat = "scalar", paired = TRUE)

    RET <- do.call("independence_test", c(list(object = ip), args))

    RET@method <- "Sign Test"
    RET@nullvalue <- 0

    return(RET)
}
