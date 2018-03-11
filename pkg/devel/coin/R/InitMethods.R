### new("ExpectCovar", ...)
setMethod("initialize",
    signature = "ExpectCovar",
    definition = function(.Object, pq = 1) {
        pq <- as.integer(pq)
        .Object@expectation <- rep(0, pq)
        .Object@covariance <- matrix(0, nrow = pq, ncol = pq)
        .Object@dimension  <- as.integer(pq)

        .Object
    }
)

### new("ExpectCovarInfluence", ...)
setMethod("initialize",
    signature = "ExpectCovarInfluence",
    definition = function(.Object, q) {
        .Object@expectation <- rep(0, q)
        .Object@covariance <- matrix(0, nrow = q, ncol = q)
        .Object@dimension  <- as.integer(q)
        .Object@sumweights <- log(1) # was 'as.double(0.0)' but there seem to be
                                     # protection issues (in 'party')
        .Object
    }
)

### new("CovarianceMatrix", ...)
setMethod("initialize",
    signature = "CovarianceMatrix",
    definition = function(.Object, covariance, ...) {
        callNextMethod(.Object, covariance = covariance, ...)
    }
)

### new("Variance", ...)
setMethod("initialize",
    signature = "Variance",
    definition = function(.Object, variance, ...) {
        callNextMethod(.Object, variance = variance, ...)
    }
)

### new("IndependenceProblem", ...)
### initialized data
setMethod("initialize",
    signature = "IndependenceProblem",
    definition = function(.Object, x, y, block = NULL, weights = NULL, ...) {

        if (NROW(x) == 0L && NROW(y) == 0L)
            stop(sQuote("x"), " and ", sQuote("y"),
                 " do not contain data")
        if (length(x) == 0L) {
            dn <- dimnames(x)
            x <- data.frame(x = rep.int(1L, nrow(x)))
            dimnames(x) <- dn
        }
        if (anyNA(x))
            stop(sQuote("x"), " contains missing values")
        if (anyNA(y))
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && !is.factor(block))
            stop(sQuote("block"), " is not a factor")
        if (!is.null(block) && anyNA(block))
            stop(sQuote("block"), " contains missing values")
        if (!is.null(weights) && anyNA(weights))
            stop(sQuote("weights"), " contains missing values")

        .Object@x <- droplevels(x)
        .Object@y <- droplevels(y)
        .Object@block <- if (is.null(block))
                             factor(rep.int(0L, nrow(x)))
                         else {
                             blockname <- attr(block, "blockname", exact = TRUE)
                             block <- droplevels(block)
                             if (!is.null(blockname))
                                 attr(block, "blockname") <- blockname
                             if (any(table(block) < 2L))
                                 stop(sQuote("block"), " contains levels with",
                                      " less than two observations")
                             block
                         }
        .Object@weights <- if (is.null(weights))
                               rep.int(1.0, nrow(x))
                           else
                               as.double(weights)

        if (!validObject(.Object))
            stop("not a valid object of class ", sQuote("IndependenceProblem"))

        .Object
    }
)

### new("IndependenceTestProblem", ...)
### set up test problem, i.e., transformations of the data
setMethod("initialize",
    signature = "IndependenceTestProblem",
    definition = function(.Object, object, xtrafo = trafo, ytrafo = trafo, ...) {

        if (!extends(class(object), "IndependenceProblem"))
            stop("Argument ", sQuote("object"), " is not of class ",
                 sQuote("IndependenceProblem"))

        tr <- check_trafo(xtrafo(object@x), ytrafo(object@y))

        .Object <- copyslots(object, .Object)
        .Object@xtrans <- tr$xtrafo
        .Object@ytrans <- tr$ytrafo
        .Object@xtrafo <- xtrafo
        .Object@ytrafo <- ytrafo

        .Object
    }
)

### new("IndependenceLinearStatistic", ...)
### compute test statistics and their expectation / covariance matrix
setMethod("initialize",
    signature = "IndependenceLinearStatistic",
    definition = function(.Object, object, varonly = FALSE, ...) {

        if (!extends(class(object), "IndependenceTestProblem"))
            stop("Argument ", sQuote("object"), " is not of class ",
                  sQuote("IndependenceTestProblem"))

        lsev <- LinStatExpCov(X = object@xtrans, Y = object@ytrans,
                              weights = as.integer(object@weights),
                              block = object@block,
                              varonly = varonly)
        nm <- statnames(object)$names # pretty names

        .Object <- copyslots(object, .Object)
        .Object@linearstatistic <- drop(lsev$LinearStatistic)
        .Object@expectation <- setNames(lsev$Expectation, nm)
        .Object@covariance <-
            if (varonly) {
                new("Variance", setNames(drop(lsev$Variance), nm))
            } else {
                cov <- matrix(0, nrow = length(nm), ncol = length(nm),
                              dimnames = list(nm, nm))
                cov[lower.tri(cov, diag = TRUE)] <- lsev$Covariance
                cov <- cov + t(cov)
                diag(cov) <- diag(cov) / 2
                new("CovarianceMatrix", cov)
            }

        if (any(variance(.Object) < eps()))
            warning("The conditional covariance matrix has ",
                    "zero diagonal elements")

        .Object
    }
)

### new("ScalarIndependenceTestStatistic", ...)
### the basis of well known univariate tests
setMethod("initialize",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(.Object, object,
        alternative = c("two.sided", "less", "greater"), paired = FALSE, ...) {

        if (!extends(class(object), "IndependenceLinearStatistic"))
            stop("Argument ", sQuote("object"), " is not of class ",
                  sQuote("IndependenceLinearStatistic"))

        ss <- (object@linearstatistic - expectation(object)) /
                  sqrt(variance(object))

        .Object <- copyslots(object, .Object)
        .Object@teststatistic <- .Object@standardizedlinearstatistic <- drop(ss)
        .Object@alternative <- match.arg(alternative)
        .Object@paired <- paired

        .Object
    }
)

### new("MaxTypeIndependenceTestStatistic", ...)
setMethod("initialize",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(.Object, object,
        alternative = c("two.sided", "less", "greater"), ...) {

        if (!extends(class(object), "IndependenceLinearStatistic"))
            stop("Argument ", sQuote("object"), " is not of class ",
                  sQuote("IndependenceLinearStatistic"))

        ss <- (object@linearstatistic - expectation(object)) /
                  sqrt(variance(object))

        .Object <- copyslots(object, .Object)
        .Object@teststatistic <-
            switch(alternative,
                "less"      = drop(min(ss)),
                "greater"   = drop(max(ss)),
                "two.sided" = drop(max(abs(ss))))
        .Object@standardizedlinearstatistic <- ss
        .Object@alternative <- match.arg(alternative)

        .Object
    }
)

### new("QuadTypeIndependenceTestStatistic", ...)
setMethod("initialize",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(.Object, object, paired = FALSE, ...) {

        if (!extends(class(object), "IndependenceLinearStatistic"))
            stop("Argument ", sQuote("object"), " is not of class ",
                  sQuote("IndependenceLinearStatistic"))

        cs <- object@linearstatistic - expectation(object)
        mp <- MPinv(covariance(object), ...)

        .Object <- copyslots(object, .Object)
        .Object@teststatistic <- drop(cs %*% mp$MPinv %*% cs)
        .Object@standardizedlinearstatistic <- cs / sqrt(variance(object))
        .Object@covarianceplus <- mp$MPinv
        .Object@df <- mp$rank
        .Object@paired <- paired

        .Object
    }
)

### new("SymmetryProblem", ...)
### initialized data
setMethod("initialize",
    signature = "SymmetryProblem",
    definition = function(.Object, x, y, block = NULL, weights = NULL, ...) {

        if (anyNA(x))
            stop(sQuote("x"), " contains missing values")
        if (!is.factor(x[[1L]]) || length(unique(table(x[[1L]]))) != 1L)
            stop(sQuote("x"), " is not a balanced factor")
        if (anyNA(y))
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && anyNA(y))
            stop(sQuote("block"), " contains missing values")

        .Object@x <- x
        .Object@y <- y
        .Object@block <- if (is.null(block))
                             factor(rep.int(seq_len(nrow(x) / nlevels(x[[1L]])),
                                            nlevels(x[[1L]])))
                         else
                             block
        .Object@weights <- if (is.null(weights))
                               rep.int(1.0, nrow(x))
                           else
                               as.double(weights)

        if (!validObject(.Object))
            stop("not a valid object of class ", sQuote("SymmetryProblem"))

        .Object
    }
)
