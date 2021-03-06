

asymptotic <- function(maxpts = 25000, abseps = 0.001, releps = 0) {
    RET <- function(object)
        AsymptNullDistribution(object, maxpts = maxpts, 
                               abseps = abseps, releps = releps)
    RET
}

approximate <- function(B = 1000) {
    RET <- function(object)
        ApproxNullDistribution(object, B = B)
    RET
}

exact <- function(algorithm = c("shift", "split-up"), fact = NULL) {
    algorithm <- match.arg(algorithm)
    RET <- function(object)
        ExactNullDistribution(object, algorithm = algorithm, fact = fact)
    RET
}

LinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_LinearStatistic", x, y, weights, 
          PACKAGE = "coin")
}

ExpectCovarInfluence <- function(y, weights) {
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_ExpectCovarInfluence", y, weights, 
          PACKAGE = "coin")
}

expectvaronly <- function(x, y, weights) {
    indx <- rep(1:nrow(x), weights)
    x <- x[indx,,drop = FALSE]
    y <- y[indx,,drop = FALSE]
    n <- nrow(x)
    Ey <- colMeans(y)
    Vy <- rowMeans((t(y) - Ey)^2)

    rSx <- colSums(x)   
    rSx2 <- colSums(x^2)
    ### in case rSx _and_ Ey are _both_ vectors
    E <- .Call("R_kronecker", Ey, rSx, package = "coin") 
          ### as.vector(kronecker(Ey, rSx))
    V <- n / (n - 1) * .Call("R_kronecker", Vy, rSx2, package = "coin") 
                        ### kronecker(Vy, rSx2)
    V <- V - 1 / (n - 1) * .Call("R_kronecker", Vy, rSx^2, package = "coin") 
                        ### kronecker(Vy, rSx^2)
    list(E = drop(E), V = matrix(V, nrow = 1))
}

ExpectCovarLinearStatistic <- function(x, y, weights, varonly = FALSE) {
    if (varonly) {
        ev <- expectvaronly(x, y, weights)
        RET <- new("ExpectCovar")
        RET@expectation <- ev$E
        RET@covariance <- ev$V
        return(RET)
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
MPinv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X))) 
        stop("X must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if (all(Positive)) 
        RET <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        RET <- array(0, dim(X)[2:1])
    else RET <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
    return(list(MPinv = RET, rank = sum(Positive)))
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

formula2data <- function(formula, data, subset, weights = NULL, ...) {

    other <- list()
    if (!is.null(weights)) other = list(weights = weights)

    ### in case `data' is an ExpressionSet object 
    if (extends(class(data), "ExpressionSet")) {
        dat <- ModelEnvFormula(formula = formula, 
                               data = Biobase::pData(Biobase::phenoData(data)),
                               subset = subset, other = other,
                               designMatrix = FALSE, responseMatrix = FALSE,
                               na.action = na.omit, 
                               ...)

        ### x are _all_ expression levels, always
        x <- as.data.frame(t(Biobase::exprs(data)))

    } else {

        dat <- ModelEnvFormula(formula = formula, 
                               data = data,
                               subset = subset, other = other, 
                               na.action = na.omit, 
                               designMatrix = FALSE, responseMatrix = FALSE,
                               ...)

        ### rhs of formula
        if (has(dat, "input"))
            x <- dat@get("input")
        else 
            stop("missing right hand side of formula")
    }

    ### ~ x + y is allowed
    if (has(dat, "response"))
        y <- dat@get("response")
    else {
        if (ncol(x) == 2) {
            y <- x[2]
            x <- x[1]
        } else 
        stop("missing left hand side of formula")
    }

    ### y ~ x | block or ~ y + x | block
    if (has(dat, "blocks")) {
        block <- dat@get("blocks")
        attr(block[[1]], "blockname") <- colnames(block)
    } else 
        block <- NULL

    RET <- list(x = x, y = y, block = block, bl = block[[1]], w = NULL)
    if (!is.null(weights)) RET$w <- dat@get("weights")[[1]]

    return(RET)
}

setscores <- function(x, scores) {

    if (is.null(scores)) return(x)

    if (!is.list(scores) || is.null(names(scores)))
       stop(sQuote("scores"), " is not a named list")

    varnames <- names(scores)

    missing <- varnames[!varnames %in% c(colnames(x@x), colnames(x@y))]
    if (length(missing) > 0)
        stop("Variable(s)", paste(missing, sep = ", "), 
             " not found in ", sQuote("x"))


    for (var in varnames) {
        if (!is.null(x@x[[var]])) {

            if (!is.factor(x@x[[var]]))
                stop(var, " is not a factor")

            if (nlevels(x@x[[var]]) != length(scores[[var]]))
                stop("scores for variable ", var, " don't match")
            
            x@x[[var]] <- ordered(x@x[[var]], levels = levels(x@x[[var]]))
            attr(x@x[[var]], "scores") <- scores[[var]]
        }
        if (!is.null(x@y[[var]])) {

            if (!is.factor(x@y[[var]]))
                stop(var, " is not a factor")

            if (nlevels(x@y[[var]]) != length(scores[[var]]))
                stop("scores for variable ", var, " don't match")
            
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
    if (!is.matrix(tx)) tx <- matrix(tx, ncol = 1)
    if (!is.matrix(ty)) ty <- matrix(ty, ncol = 1)
    storage.mode(tx) <- "double"
    storage.mode(ty) <- "double"
    list(xtrafo = tx, ytrafo = ty)
}

table2df <- function(x) {
    if (!is.table(x))
        stop(sQuote("x"), " is not of class ", sQuote("table"))
    x <- as.data.frame(x)
    freq <- x[["Freq"]]
    x <- x[rep(1:nrow(x), freq), ,drop = FALSE]
    rownames(x) <- 1:nrow(x)
    return(x[,colnames(x) != "Freq"])
}

table2df_sym <- function(x) {
    x <- table2df(x)
    lx <- levels(x[[1]])
    if (!all(sapply(x, function(x) all(levels(x) == lx))))
        stop("table ", sQuote("x"), " does not represent a symmetry problem")
    n <- nrow(x)
    p <- ncol(x)
    y <- data.frame(groups = factor(rep(colnames(x), rep(n, p))),
                    response = factor(unlist(x), labels = lx))
    rownames(y) <- 1:(n*p)
    y
}

table2IndependenceProblem <- function(object) {

    df <- as.data.frame(object)
    if (ncol(df) == 3)
        ip <- new("IndependenceProblem", x = df[1], y = df[2],
                  block = NULL, weights = df[["Freq"]])
    if (ncol(df) == 4) {
        attr(df[[3]], "blockname") <- colnames(df)[3]
        ip <- new("IndependenceProblem", x = df[1], y = df[2],
                  block = df[[3]], weights = df[["Freq"]])
    }
    ip
}

is_2sample <- function(object) {
    groups <- nlevels(object@x[[1]]) == 2 && ncol(object@xtrans) == 1
    return(is_Ksample(object) && groups)
}

is_Ksample <- function(object) {
    groups <- (ncol(object@x) == 1 && is.factor(object@x[[1]]))
#    values <- (ncol(object@y) == 1 && ncol(object@ytrans) == 1)
    values <- ncol(object@ytrans) == 1
    return(groups && values)
}

is_numeric_y <- function(object) {
    is.numeric(object@y[[1]])
}

is_censored_y <- function(object) {
    ncol(object@y) == 1 && class(object@y[[1]]) == "Surv"
}

is_corr <- function(object) {
    (is.numeric(object@x[[1]]) && is.numeric(object@y[[1]])) &&
     (ncol(object@xtrans) == 1 && ncol(object@ytrans) == 1)
}

is_contingency <- function(object) {
    x <- object@x
    groups <- (ncol(x) == 1 && is.factor(x[[1]]))
    y <- object@y
    values <- (ncol(y) == 1 && is.factor(y[[1]]))
    ### trans <- all(rowSums(object@xtrans) %in% c(0,1)) && 
    ###          all(rowSums(object@ytrans) %in% c(0, 1))
    ### hm, two ordinal variables are a contingency problem as well (???)
    return((groups && values)) ### && trans)       
}

is_ordered <- function(object) {
    x <- object@x
    y <- object@y
    (is_Ksample(object) || is_contingency(object)) && 
    (is.ordered(x[[1]]) || is.ordered(y[[1]]))
}

is_completeblock <- function(object) {
    all(table(object@x[[1]], object@block) == 1)
}

is_scalar <- function(object) {
    ncol(object@xtrans) == 1 && ncol(object@ytrans) == 1
}

is_ordered_x <- function(object) {
    all(sapply(object@x, function(x) is.numeric(x) || is.ordered(x)))
}

is_integer <- function(x, fact = c(1, 2, 10, 100, 1000))
    sapply(fact, function(f) max(abs(round(x * f) - (x * f))) < eps())

is_censored <- function(object) {
    ncol(object@y) == 1 && class(object@y[[1]]) == "Surv"
}

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
        distribution <- eval(parse(text = 
                                   paste(distribution, "()", sep = "")))
    }  
    distribution
}

statnames <- function(object) {
    nc <- ncol(object@ytrans)
    nr <- ncol(object@xtrans)
    dn <- list(colnames(object@xtrans),
               colnames(object@ytrans))
    if (is.null(dn[[1]])) {
        if (nr == 1) {
            dn[[1]] <- ""
        } else {
            dn[[1]] <- paste("X", 1:nr, sep = "")
        }
    }
    if (is.null(dn[[2]])) {
        if (nc == 1) {
            dn[[2]] <- ""
        } else {
            dn[[2]] <- paste("Y", 1:nc, sep = "")
        }
    }
    list(dimnames = dn, 
         names = paste(rep((dn[[1]]), nc), 
                       rep((dn[[2]]), rep(nr, nc)), 
                       sep = ifelse(dn[[1]] == "" || dn[[2]] == "", "", ":")))
}

eps <- function() sqrt(.Machine$double.eps)

GE <- function(x, y)
    x > y | abs(x - y) < eps()

LE <- function(x, y)
    x < y | abs(x - y) < eps()

### don't use! never!
get_weights <- function(object) object@statistic@weights
get_xtrans <- function(object) object@statistic@xtrans
get_ytrans <- function(object) object@statistic@ytrans

chkone <- function(w)
    !(max(abs(w - 1.0)) < eps())

ranktests <- function()
    c("wilcox_test", "normal_test", "median_test",
      "ansari_test", "surv_test", "kruskal_test",
      "fligner_test", "spearman_test", "friedman_test",
      "wilcoxsign_test")
