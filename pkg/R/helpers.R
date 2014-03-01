asymptotic <- function(maxpts = 25000, abseps = 0.001, releps = 0) {
    function(object)
        AsymptNullDistribution(object, maxpts = maxpts,
                               abseps = abseps, releps = releps)
}

approximate <- function(B = 10000) {
    function(object)
        ApproxNullDistribution(object, B = B)
}

exact <- function(algorithm = c("auto", "shift", "split-up"), fact = NULL) {
    algorithm <- match.arg(algorithm)
    function(object)
        ExactNullDistribution(object, algorithm = algorithm, fact = fact)
}

LinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_LinearStatistic", x, y, weights, PACKAGE = "coin")
}

ExpectCovarInfluence <- function(y, weights) {
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_ExpectCovarInfluence", y, weights, PACKAGE = "coin")
}

expectvaronly <- function(x, y, weights) {
    indx <- rep.int(seq_len(nrow(x)), weights)
    x <- x[indx, , drop = FALSE]
    y <- y[indx, , drop = FALSE]
    n <- nrow(x)
    Ey <- colMeans(y)
    Vy <- rowMeans((t(y) - Ey)^2)

    rSx <- colSums(x)
    rSx2 <- colSums(x^2)
    ## in case rSx _and_ Ey are _both_ vectors
    E <- .Call("R_kronecker", Ey, rSx, PACKAGE = "coin")
    V <- n / (n - 1) * .Call("R_kronecker", Vy, rSx2, PACKAGE = "coin")
    V <- V - 1 / (n - 1) * .Call("R_kronecker", Vy, rSx^2, PACKAGE = "coin")
    list(E = drop(E), V = matrix(V, nrow = 1L))
}

ExpectCovarLinearStatistic <- function(x, y, weights, varonly = FALSE) {
    if (varonly) {
        ev <- expectvaronly(x, y, weights)
        new("ExpectCovar", expectation = ev$E, covariance = ev$V)
    } else {
        storage.mode(x) <- "double"
        storage.mode(y) <- "double"
        storage.mode(weights) <- "double"
        expcovinf <- ExpectCovarInfluence(y, weights)
        .Call("R_ExpectCovarLinearStatistic", x, y, weights, expcovinf,
               PACKAGE = "coin")
    }
}

### copied from package MASS
MPinv <- function (X, tol = eps())
{
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
        stop("X must be a numeric or complex matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    RET <- if (all(Positive))
               Xsvd$v %*% (1 / Xsvd$d * t(Xsvd$u))
           else if (!any(Positive))
               array(0, dim(X)[2L:1L])
           else
               Xsvd$v[, Positive, drop = FALSE] %*%
                 ((1 / Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    list(MPinv = RET, rank = sum(Positive))
}

copyslots <- function(source, target) {
    slots <- names(getSlots(class(source)))
    slots <- slots[(slots %in% names(getSlots(class(target))))]
    if (length(slots) == 0)
        stop("no common slots to copy to")
    for (s in slots)
        eval(parse(text = paste("target@", s, " <- source@", s)))
    return(target)
}

ft <- function(test, class, formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    object <- formula2data(formula, data, subset, weights = weights, ...)
    object <- new(class, x = object$x, y = object$y, block = object$bl,
                  weights = object$w)
    args <- list(...)
    args$frame <- NULL

    ## warn users of weighted rank tests
    if (test %in% ranktests && !is.null(object@weights) &&
        !is_unity(object@weights))
        warning("Rank transformation doesn't take weights into account")

    do.call(test, c(list(object = object), args))
}

ranktests <-
    c("wilcox_test", "kruskal_test", "normal_test", "median_test",
      "savage_test", "taha_test", "klotz_test", "mood_test", "ansari_test",
      "fligner_test", "conover_test", "surv_test", "quade_test",
      "friedman_test", "wilcoxsign_test", "spearman_test", "fisyat_test",
      "quadrant_test", "koziol_test")

formula2data <- function(formula, data, subset, weights = NULL, ...) {

    other <- list()
    if (!is.null(weights)) other = list(weights = weights)

    ## in case `data' is an ExpressionSet object
    if (extends(class(data), "ExpressionSet")) {
        dat <- ModelEnvFormula(formula = formula,
                               data = Biobase::pData(Biobase::phenoData(data)),
                               subset = subset, other = other,
                               designMatrix = FALSE, responseMatrix = FALSE,
                               na.action = na.omit,
                               ...)

        ## x are _all_ expression levels, always
        x <- as.data.frame(t(Biobase::exprs(data)))
    } else {
        dat <- ModelEnvFormula(formula = formula,
                               data = data,
                               subset = subset, other = other,
                               na.action = na.omit,
                               designMatrix = FALSE, responseMatrix = FALSE,
                               ...)

        ## rhs of formula
        if (has(dat, "input"))
            x <- dat@get("input")
        else
            stop("missing right hand side of formula")
    }

    ## ~ x + y is allowed
    if (has(dat, "response"))
        y <- dat@get("response")
    else {
        if (ncol(x) == 2L) {
            y <- x[2L]
            x <- x[1L]
        } else
            stop("missing left hand side of formula")
    }

    ## y ~ x | block or ~ y + x | block
    if (has(dat, "blocks")) {
        block <- dat@get("blocks")
        attr(block[[1L]], "blockname") <- colnames(block)
    } else
        block <- NULL

    RET <- list(x = x, y = y, block = block, bl = block[[1L]], w = NULL)
    if (!is.null(weights))
        RET$w <- dat@get("weights")[[1L]]

    return(RET)
}

setscores <- function(x, scores) {

    if (is.null(scores)) return(x)

    varnames <- names(scores)
    if (!is.list(scores) || is.null(varnames))
       stop(sQuote("scores"), " is not a named list")

    missing <- varnames[!varnames %in% c(colnames(x@x), colnames(x@y))]
    if (length(missing) > 0L)
        stop("Variable(s) ", paste(missing, sep = ", "),
             " not found in ", sQuote("x"))

    for (var in varnames) {
        if (!is.null(x@x[[var]])) {
            if (!is.factor(x@x[[var]]))
                stop(sQuote(var), " is not a factor")
            if (nlevels(x@x[[var]]) != length(scores[[var]]))
                stop("scores for variable ", sQuote(var), " don't match")
            if (!is_monotone(scores[[var]]))
                stop("scores for variable ", sQuote(var), " aren't monotone")
            x@x[[var]] <- ordered(x@x[[var]], levels = levels(x@x[[var]]))
            attr(x@x[[var]], "scores") <- scores[[var]]
        }
        if (!is.null(x@y[[var]])) {
            if (!is.factor(x@y[[var]]))
                stop(sQuote(var), " is not a factor")
            if (nlevels(x@y[[var]]) != length(scores[[var]]))
                stop("scores for variable ", sQuote(var), " don't match")
            if (!is_monotone(scores[[var]]))
                stop("scores for variable ", sQuote(var), " aren't monotone")
            x@y[[var]] <- ordered(x@y[[var]], levels = levels(x@y[[var]]))
            attr(x@y[[var]], "scores") <- scores[[var]]
        }
    }
    return(x)
}

### user-supplied trafo functions may return a vector or matrices
### with NROW being equal for the x and y variables
check_trafo <- function(tx, ty) {

    if (!(is.numeric(tx) || is.logical(tx)))
        stop(sQuote("xtrafo"), " does not return a numeric or logical vector")
    if (!(is.numeric(ty) || is.logical(ty)))
        stop(sQuote("ytrafo"), " does not return a numeric or logical vector")
    if (NROW(tx) != NROW(ty))
        stop("Dimensions of returns of ", sQuote("xtrafo"), " and ",
             sQuote("ytrafo"), " don't match")
    if (!is.matrix(tx)) tx <- matrix(tx, ncol = 1L)
    if (!is.matrix(ty)) ty <- matrix(ty, ncol = 1L)
    storage.mode(tx) <- "double"
    storage.mode(ty) <- "double"
    list(xtrafo = tx, ytrafo = ty)
}

table2df <- function(x) {
    if (!is.table(x))
        stop(sQuote("x"), " is not of class ", sQuote("table"))
    x <- as.data.frame(x)
    freq <- x[["Freq"]]
    x <- x[rep.int(seq_len(nrow(x)), freq), , drop = FALSE]
    rownames(x) <- seq_len(nrow(x))
    return(x[, colnames(x) != "Freq"])
}

table2df_sym <- function(x) {
    x <- table2df(x)
    lx <- levels(x[[1L]])
    if (!all(vapply(x, function(x) all(levels(x) == lx), NA)))
        stop("table ", sQuote("x"), " does not represent a symmetry problem")
    data.frame(conditions = factor(rep.int(colnames(x),
                                           rep.int(nrow(x), ncol(x)))),
               response = factor(unlist(x, recursive = FALSE,
                                        use.names = FALSE),
                                 labels = lx))
}

table2IndependenceProblem <- function(object) {

    df <- as.data.frame(object)
    if (ncol(df) == 3L)
        ip <- new("IndependenceProblem", x = df[1L], y = df[2L],
                  block = NULL, weights = df[["Freq"]])
    if (ncol(df) == 4L) {
        attr(df[[3L]], "blockname") <- colnames(df)[3L]
        ip <- new("IndependenceProblem", x = df[1L], y = df[2L],
                  block = df[[3]], weights = df[["Freq"]])
    }
    ip
}

is_2sample <- function(object) {
    groups <- nlevels(droplevels(object@x)[[1L]]) == 2L &&
                  ncol(object@xtrans) == 1L
    return(is_Ksample(object) && groups)
}

is_Ksample <- function(object) {
    groups <- (ncol(object@x) == 1L && is.factor(object@x[[1L]]))
###    values <- (ncol(object@y) == 1L && ncol(object@ytrans) == 1L)
    values <- ncol(object@ytrans) == 1L
    return(groups && values)
}

is_numeric_y <- function(object)
    ncol(object@y) == 1L && is.numeric(object@y[[1L]])

is_censored_y <- function(object)
    ncol(object@y) == 1L && is.Surv(object@y[[1L]])

is_ordered_x <- function(object)
###    all(vapply(object@x, function(x) is.numeric(x) || is.ordered(x), NA))
    ncol(object@x) == 1L && is.ordered(object@x[[1L]])

is_corr <- function(object)
    (is.numeric(object@x[[1L]]) && is.numeric(object@y[[1L]])) &&
        (ncol(object@xtrans) == 1L && ncol(object@ytrans) == 1L)

## is_contingency <- function(object) {
##     x <- object@x
##     groups <- (ncol(x) == 1L && is.factor(x[[1L]]))
##     y <- object@y
##     values <- (ncol(y) == 1L && is.factor(y[[1L]]))
##     ### trans <- all(rowSums(object@xtrans) %in% c(0, 1)) &&
##     ###          all(rowSums(object@ytrans) %in% c(0, 1))
##     ### hm, two ordinal variables are a contingency problem as well (???)
##     return((groups && values)) ### && trans)
## }

is_contingency <- function(object)
    (ncol(object@x) == 1L && is.factor(object@x[[1L]])) &&
        (ncol(object@y) == 1L && is.factor(object@y[[1L]]))

is_ordered <- function(object)
    (is_Ksample(object) || is_contingency(object)) &&
        (is.ordered(object@x[[1L]]) || is.ordered(object@y[[1L]]))

is_singly_ordered <- function(object) {
    x <- object@x[[1L]]
    y <- object@y[[1L]]
    (is_Ksample(object) || is_contingency(object)) &&
        ((is.ordered(x) && (is.numeric(y) || (!is.ordered(y) && nlevels(y) > 2L))) ||
         (is.ordered(y) && (is.numeric(x) || (!is.ordered(x) && nlevels(x) > 2L))))
}

is_doubly_ordered <- function(object) {
    x <- object@x[[1L]]
    y <- object@y[[1L]]
    (is_Ksample(object) || is_contingency(object)) &&
        ((is.ordered(x) && is.ordered(y)) ||
         ((is.ordered(x) && nlevels(y) == 2L) ||
          (is.ordered(y) && nlevels(x) == 2L)))
}

is_completeblock <- function(object)
    all(table(object@x[[1L]], object@block) == 1L)

is_scalar <- function(object)
    ncol(object@xtrans) == 1L && ncol(object@ytrans) == 1L

is_integer <- function(x, fact = NULL) {
    if (is.null(fact))
        fact <- c(1, 2, 10, 100, 1000, 10000, 100000)
    f <- vapply(fact, function(f) max(abs(round(x * f) - (x * f))) < eps(), NA)
    if (RET <- any(f))
        attr(RET, "fact") <- min(fact[f])
    RET
}

is_monotone <- function (x)
    all(x == cummax(x)) || all(x == cummin(x))

isequal <- function(a, b) {
    attributes(a) <- NULL
    attributes(b) <- NULL
    if (!isTRUE(all.equal(a, b))) {
        print(a, digits = 10)
        print(b, digits = 10)
        return(FALSE)
    } else {
        return(TRUE)
    }
}

check_distribution_arg <- function(distribution,
    values = c("asymptotic", "approximate", "exact")) {
    if (is.character(distribution)) {
        distribution <- match.arg(distribution[1], values)
        distribution <- eval(parse(text = paste0(distribution, "()")))
    }
    distribution
}

setup_args <- function(...) {
    cl <- match.call(independence_test.IndependenceProblem,
                     call = sys.call(sys.parent()), expand.dots = FALSE)
    ## find default arguments and values
    args <- formals(independence_test.IndependenceProblem)
    args$object <- args$... <- NULL
    nm <- names(args)
    ## replace default values with user-specified values
    for (i in nm[nm %in% names(cl)])
        args[[i]] <- cl[[i]]
    ## override default and user-specified values
    for (i in nm[nm %in% names(list(...))])
        args[[i]] <- list(...)[[i]]
    lapply(args, eval.parent)
}

statnames <- function(object) {
    nc <- ncol(object@ytrans)
    nr <- ncol(object@xtrans)
    dn <- list(colnames(object@xtrans),
               colnames(object@ytrans))
    if (is.null(dn[[1L]])) {
        if (nr == 1L) {
            dn[[1L]] <- ""
        } else {
            dn[[1L]] <- paste0("X", seq_len(nr))
        }
    }
    if (is.null(dn[[2L]])) {
        if (nc == 1L) {
            dn[[2L]] <- ""
        } else {
            dn[[2L]] <- paste0("Y", seq_len(nc))
        }
    }
    list(dimnames = dn,
         names = paste(rep.int((dn[[1L]]), nc),
                       rep.int((dn[[2L]]), rep.int(nr, nc)),
                       sep = ifelse(dn[[1L]] == "" || dn[[2L]] == "", "", ":")))
}

eps <- function() sqrt(.Machine$double.eps)

EQ <- function(x, y)
    abs(x - y) < eps()

GE <- function(x, y)
    x > y | abs(x - y) < eps()

LE <- function(x, y)
    x < y | abs(x - y) < eps()

### don't use! never!
get_weights <- function(object) object@statistic@weights
get_xtrans <- function(object) object@statistic@xtrans
get_ytrans <- function(object) object@statistic@ytrans

is_unity <- function(x)
    max(abs(x - 1.0)) < eps()

setColnames <- function (object, nm) {
    colnames(object) <- nm
    object
}

if(getRversion() < "2.15") paste0 <- function(...) paste(..., sep = "")
