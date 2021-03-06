
ft <- function(test, formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    d <- formula2data(formula, data, subset, weights = weights, ...)
    ip <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl,
              weights = d$w)
    args <- list(...)
    args$frame <- NULL

    ### warn users of weighted rank tests
    w <- d$w
    if (is.null(w)) w <- 1L
    if (test %in% ranktests() & (max(abs(w - 1)) > .Machine$double.eps))
        warning("Rank transformation doesn't take weights into account")

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
    ip <- table2IndependenceProblem(object)
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
    its <- new("IndependenceTestStatistic", itp, varonly = TRUE) 
#        varonly = class(distribution) == "approximate" && teststat == "max")

    ### compute test statistic and corresponding null distribution
    RET <- switch(teststat,
        "scalar" = {
            ts <- new("ScalarIndependenceTestStatistic", its, 
                      alternative = alternative)
            nd <- distribution(ts)
            new("ScalarIndependenceTest", statistic = ts, distribution = nd)
        },
        "max" = {
            ts <- new("MaxTypeIndependenceTestStatistic", its, 
                      alternative = alternative)
            nd <- distribution(ts)
            new("MaxTypeIndependenceTest", statistic = ts, distribution = nd)
        },
        "quad" = {
            ts <- new("QuadTypeIndependenceTestStatistic", its)
            nd <- distribution(ts)
            new("QuadTypeIndependenceTest", statistic = ts, 
                distribution = nd)
        })

    ### return object inheriting from class `IndependenceTest'
    return(RET)
}

