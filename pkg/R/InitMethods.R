
### new("CovarianceMatrix", ...)
setMethod(f = "initialize",
    signature = "CovarianceMatrix",
    definition = function(.Object, x) {
        .Object@covariance <- x
        .Object
    }
)

### new("Variance", ...)
setMethod(f = "initialize",
    signature = "Variance",
    definition = function(.Object, x) {
        .Object@variance <- x
        .Object
    }
)

### new("IndependenceProblem", ...)
### initialized data
setMethod(f = "initialize",
    signature = "IndependenceProblem",
    definition = function(.Object, x, y, block = NULL, weights = NULL) {

        if (NROW(x) == 0 && NROW(y) == 0)
            stop(sQuote("x"), " and ", sQuote("y"),
                 " do not contain data")
        if (length(x) == 0) {
            dn <- dimnames(x)
            x <- data.frame(x = rep.int(1, nrow(x)))
            dimnames(x) <- dn
        }
        if (any(is.na(x)))
            stop(sQuote("x"), " contains missing values")
        if (any(is.na(y)))
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && !is.factor(block))
            stop(sQuote("block"), " is not a factor")
        if (!is.null(block) && any(is.na(block)))
            stop(sQuote("block"), " contains missing values")
        if (!is.null(weights) && any(is.na(weights)))
            stop(sQuote("weights"), " contains missing values")
        .Object@x <- x
        .Object@y <- y
        if (is.null(block)) {
            .Object@block <- factor(rep.int(0, nrow(x)))
        } else {
            if (any(table(block) < 2))
                stop(sQuote("block"),
                     " contains levels with less than two observations")
            .Object@block <- block
        }
        if (is.null(weights)) {
            .Object@weights <- rep.int(1.0, nrow(x))
        } else {
            .Object@weights <- as.double(weights)
        }
        if (!validObject(.Object))
            stop("not a valid object of class ",
                 sQuote("IndependenceProblem"))
        .Object
    }
)

### new("IndependenceTestProblem", ...)
### set up test problem, i.e., transformations of the data
setMethod(f = "initialize",
    signature = "IndependenceTestProblem",
    definition = function(.Object, ip, xtrafo = trafo, ytrafo = trafo, ...) {

        if (!extends(class(ip), "IndependenceProblem"))
            stop("Argument ", sQuote("ip"), " is not of class ",
                  sQuote("IndependenceProblem"))

        .Object <- copyslots(ip, .Object)

        tr <- check_trafo(xtrafo(ip@x), ytrafo(ip@y))
        .Object@xtrans <- tr$xtrafo
        .Object@ytrans <- tr$ytrafo
        .Object@xtrafo <- xtrafo
        .Object@ytrafo <- ytrafo

        .Object
    }
)

### new("IndependenceLinearStatistic", ...)
### compute test statistics and their expectation / covariance matrix
setMethod(f = "initialize",
    signature = "IndependenceLinearStatistic",
    definition = function(.Object, itp, varonly = FALSE) {

        if (!extends(class(itp), "IndependenceTestProblem"))
            stop("Argument ", sQuote("itp"), " is not of class ",
                  sQuote("IndependenceTestProblem"))

        .Object <- copyslots(itp, .Object)

        .Object@linearstatistic <-
            drop(LinearStatistic(itp@xtrans, itp@ytrans, itp@weights))

        ### <REMINDER>
        ### for teststat = "max" and distribution = "approx"
        ### we don't need to covariance matrix but the variances only
        ### </REMINDER>

        ### possibly stratified by block
        if (nlevels(itp@block) == 1) {
            expcov <-
                ExpectCovarLinearStatistic(itp@xtrans, itp@ytrans, itp@weights,
                                           varonly = varonly)
            exp <- expcov@expectation
            cov <- expcov@covariance
        } else {
            exp <- 0
            cov <- 0
            for (lev in levels(itp@block)) {
                indx <- (itp@block == lev)
                ec <- ExpectCovarLinearStatistic(itp@xtrans[indx,,drop = FALSE],
                                                 itp@ytrans[indx,,drop = FALSE],
                                                 itp@weights[indx],
                                                 varonly = varonly)
                exp <- exp + ec@expectation
                cov <- cov + ec@covariance
            }
        }

        .Object@expectation <- drop(exp)
        if (varonly) {
            .Object@covariance <- new("Variance", drop(cov))
        } else {
            .Object@covariance <- new("CovarianceMatrix", cov)
        }

        ### pretty names
        nm <- statnames(itp)$names
        names(.Object@expectation) <- nm

        if (extends(class(.Object@covariance), "CovarianceMatrix")) {
                dimnames(.Object@covariance@covariance) <- list(nm, nm)
        }
        if (extends(class(.Object@covariance), "Variance"))
                names(.Object@covariance@variance) <- nm

        if (any(variance(.Object) < eps()))
            warning("The conditional covariance matrix has ",
                    "zero diagonal elements")
        .Object
    }
)


### new("IndependenceTestStatistic", ...)
### compute test statistics and their expectation / covariance matrix
setMethod(f = "initialize",
    signature = "IndependenceTestStatistic",
    definition = function(.Object, itp, varonly = FALSE) {

        copyslots(new("IndependenceLinearStatistic", itp, varonly = varonly),
                  .Object)
    }
)

### new("ScalarIndependenceTestStatistic", ...)
### the basis of well known univariate tests
setMethod(f = "initialize",
    signature = "ScalarIndependenceTestStatistic",
    definition = function(.Object, its,
        alternative = c("two.sided", "less", "greater"), paired = FALSE) {

        if (!extends(class(its), "IndependenceTestStatistic"))
            stop("Argument ", sQuote("its"), " is not of class ",
                  sQuote("IndependenceTestStatistic"))

        .Object <- copyslots(its, .Object)
        .Object@alternative <- match.arg(alternative)
        .Object@paired <- paired

        standstat <- (its@linearstatistic - expectation(its)) /
                     sqrt(variance(its))
        .Object@teststatistic <- drop(standstat)
        .Object@standardizedlinearstatistic <- drop(standstat)

        .Object
    }
)

### new("MaxTypeIndependenceTestStatistic", ...)
setMethod(f = "initialize",
    signature = "MaxTypeIndependenceTestStatistic",
    definition = function(.Object, its,
        alternative = c("two.sided", "less", "greater")) {

        if (!extends(class(its), "IndependenceTestStatistic"))
            stop("Argument ", sQuote("its"), " is not of class ",
                  sQuote("IndependenceTestStatistic"))

        .Object <- copyslots(its, .Object)

        .Object@alternative <- match.arg(alternative)
        standstat <- (its@linearstatistic - expectation(its)) /
                      sqrt(variance(its))
        .Object@teststatistic <- switch(alternative,
            "less" = drop(min(standstat)),
            "greater" = drop(max(standstat)),
            "two.sided" = drop(max(abs(standstat)))
         )
        .Object@standardizedlinearstatistic <- standstat

        .Object
    }
)

### new("QuadTypeIndependenceTestStatistic", ...)
setMethod(f = "initialize",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(.Object, its, ...) {

        if (!extends(class(its), "IndependenceTestStatistic"))
            stop("Argument ", sQuote("its"), " is not of class ",
                  sQuote("IndependenceTestStatistic"))

        .Object <- copyslots(its, .Object)

        mp <- MPinv(covariance(its), ...)
        .Object@covarianceplus <- mp$MPinv
        .Object@df <- mp$rank

        stand <- (its@linearstatistic - expectation(its))
        .Object@teststatistic <-
            drop(stand %*% .Object@covarianceplus %*% stand)
        .Object@standardizedlinearstatistic <-
            (its@linearstatistic - expectation(its)) / sqrt(variance(its))

        .Object
    }
)

### new("SymmetryProblem", ...)
### initialized data
setMethod(f = "initialize",
    signature = "SymmetryProblem",
    definition = function(.Object, x, y, block = NULL, weights = NULL) {

        if (any(is.na(x)))
            stop(sQuote("x"), " contains missing values")
        if (!is.factor(x[[1]]) || length(unique(table(x[[1]]))) != 1)
            stop(sQuote("x"), " is not a balanced factor")
        if (any(is.na(y)))
            stop(sQuote("y"), " contains missing values")
        if (!is.null(block) && any(is.na(y)))
            stop(sQuote("block"), " contains missing values")

        .Object@x <- x
        .Object@y <- y

        if (is.null(block))
            .Object@block <- factor(rep.int(seq_len(nrow(x) / nlevels(x[[1]])),
                                            nlevels(x[[1]])))
        else
            .Object@block <- block

        if (is.null(weights))
            .Object@weights <- rep.int(1.0, nrow(x))
        else
            .Object@weights <- as.double(weights)

        .Object
    }
)
