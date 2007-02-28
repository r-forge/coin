
ft <- function(test, formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    d <- formula2data(formula, data, subset, weights = weights, ...)
    ip <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl,
              weights = d$w)
    args <- list(...)
    args$frame <- NULL
    RET <- do.call(test, c(list(object = ip), args))
    return(RET)
}

### a generic test procedure for classical (and not so classical) tests
independence_test <- function(object, ...) UseMethod("independence_test")

independence_test.formula <- function(formula, data = list(), subset = NULL, 
    weights = NULL, ...) {

    ft("independence_test", formula, data, subset, weights, 
       frame = parent.frame(), ...)
}

independence_test.table <- function(object, 
    distribution = c("asymptotic", "approximate"), ...) {

    distribution <- check_distribution_arg(distribution, 
                                           c("asymptotic", "approximate"))
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (class(distribution) == "asymptotic") {
        df <- as.data.frame(object)
        if (ncol(df) == 3)
            ip <- new("IndependenceProblem", x = df[1], y = df[2], block = NULL, 
                      weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = df[[3]], weights = df[["Freq"]])
        }
    } else {
        df <- table2df(object)
        if (ncol(df) == 3) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = df[[3]])
        } else {
            ip <- new("IndependenceProblem", x = df[1], y = df[2], 
                      block = NULL) 
        }
    }
    ### </FIXME>
    RET <- do.call("independence_test", 
                   c(list(object = ip, distribution = distribution), 
                   list(...)))
    return(RET)
}

independence_test.IndependenceProblem <- function(object,
    teststat = c("max", "quad", "scalar"),
    distribution = c("asymptotic", "approximate", "exact"), 
    alternative = c("two.sided", "less", "greater"), 
    xtrafo = trafo, ytrafo = trafo, scores = NULL, check = NULL, ...) {

    addargs <- list(...)
    if (length(addargs) > 0) 
        warning("additional arguments ", 
                paste(names(addargs), collapse = ", "),
                " will be ignored")

    ### just for backward compatibility
    teststat <- match.arg(teststat, choices = c("maxtype", "quadtype", "scalar"), 
                          several.ok = TRUE)
    if (teststat[1] == "maxtype") teststat <- "max"
    if (teststat[1] == "quadtype") teststat <- "quad"
    teststat <- match.arg(teststat)
    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    ### expand weights if conditional MC is requested
    if (class(distribution) == "approximate") {
        w <- object@weights
        if (chkone(w)) {
            indx <- rep(1:length(w), w)
            object <- new("IndependenceProblem", 
                          x = object@x[indx,,drop = FALSE],
                          y = object@y[indx,,drop = FALSE], 
                          block = object@block[indx])
        }
    }

    ### convert factors to ordered and attach scores if requested
    object <- setscores(object, scores)

    ### transform data if requested and setup a test problem
    itp <- new("IndependenceTestProblem", object, xtrafo = xtrafo, 
               ytrafo = ytrafo, ...)

    if (!is.null(check)) {
        if (is.function(check)) {
            if (!check(itp))
                stop(sQuote("check"), " failed")
        } else {
            stop(sQuote("check"), " is not a function")
        }
    }

    ### check type of test statistic and alternative
    scalar <- is_scalar(itp)

    if (!scalar) {
        if (teststat == "scalar") {
            warning("Length linear statistic > 1, using ",
                    sQuote("max"), "-type test statistic")
            teststat <- "max"
        }
    } else {
        if (teststat == "max") teststat <- "scalar"
    }
    if (alternative != "two.sided" && teststat == "quad")
        warning(sQuote("alternative"), " is ignored for ", 
                teststat, " type test statistics")

    ### compute linear statistic, conditional expectation and
    ### conditional covariance
    its <- new("IndependenceTestStatistic", itp, 
        varonly = class(distribution) == "approximate" && teststat == "max")

    ### compute test statistic and corresponding null distribution
    RET <- switch(teststat,
        "scalar" = {
            ts <- new("ScalarIndependenceTestStatistic", its, 
                      alternative = alternative)

            nd <- switch(class(distribution),
                "asymptotic" = do.call("AsymptNullDistribution", 
                                       c(list(object = ts), distribution)),
                "exact"  = do.call("ExactNullDistribution", 
                                   c(list(object = ts), distribution)),
                "approximate" = do.call("ApproxNullDistribution", 
                                        c(list(object = ts), distribution))
            )
            new("ScalarIndependenceTest", statistic = ts, distribution = nd)
        },
        "max" = {
            ts <- new("MaxTypeIndependenceTestStatistic", its, 
                      alternative = alternative)
            nd <- switch(class(distribution),
                "asymptotic" = do.call("AsymptNullDistribution", 
                                       c(list(object = ts), distribution)),
                "exact"  = do.call("ExactNullDistribution", 
                                   c(list(object = ts), distribution)),
                "approximate" = do.call("ApproxNullDistribution", 
                                        c(list(object = ts), distribution))
                      )
            new("MaxTypeIndependenceTest", statistic = ts, distribution = nd)
        },
        "quad" = {
            ts <- new("QuadTypeIndependenceTestStatistic", its)
            nd <- switch(class(distribution),
                "asymptotic" = do.call("AsymptNullDistribution", 
                                       c(list(object = ts), distribution)),
                "exact"  = do.call("ExactNullDistribution", 
                                   c(list(object = ts), distribution)),
                "approximate" = do.call("ApproxNullDistribution", 
                                        c(list(object = ts), distribution))
            )
            new("QuadTypeIndependenceTest", statistic = ts, 
                distribution = nd)
        })

    RET@method <- "General Independence Test"

    ### return object inheriting from class `IndependenceTest'
    return(RET)
}

