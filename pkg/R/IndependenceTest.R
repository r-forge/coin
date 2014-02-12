### a generic test procedure for classical (and not so classical) tests
independence_test <- function(object, ...) UseMethod("independence_test")

independence_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("independence_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

independence_test.table <- function(object, ...) {

    object <- table2IndependenceProblem(object)
    object <- do.call("independence_test", c(list(object = object), list(...)))
    return(object)
}

independence_test.IndependenceProblem <- function(object,
    teststat = c("max", "quad", "scalar"),
    distribution = c("asymptotic", "approximate", "exact"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, scores = NULL, check = NULL, paired = FALSE,
    ...) {

    addargs <- list(...)
    if (length(addargs) > 0)
        warning("additional arguments ",
                paste0(names(addargs), collapse = ", "),
                " will be ignored")

    teststat <- match.arg(teststat)
    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    ## convert factors to ordered and attach scores if requested
    if (!is.null(scores))
        object <- setscores(object, scores)

    ## transform data if requested and setup a test problem
    object <- new("IndependenceTestProblem", object, xtrafo = xtrafo,
                  ytrafo = ytrafo, ...)

    if (!is.null(check)) {
        if (is.function(check)) {
            if (!check(object))
                stop(sQuote("check"), " failed")
        } else {
            stop(sQuote("check"), " is not a function")
        }
    }

    ## check type of test statistic and alternative
    scalar <- is_scalar(object)

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

    ## compute linear statistic, conditional expectation and
    ## conditional covariance
    object <- new("IndependenceTestStatistic", object, varonly = TRUE)
###         varonly = class(distribution) == "approximate" && teststat == "max")

    ## compute test statistic and corresponding null distribution
    object <- switch(teststat,
        "scalar" = {
            object <- new("ScalarIndependenceTestStatistic", object,
                          alternative = alternative, paired = isTRUE(paired))
            new("ScalarIndependenceTest", statistic = object,
                distribution = distribution(object))
        },
        "max" = {
            object <- new("MaxTypeIndependenceTestStatistic", object,
                          alternative = alternative)
            new("MaxTypeIndependenceTest", statistic = object,
                distribution = distribution(object))
        },
        "quad" = {
            object <- new("QuadTypeIndependenceTestStatistic", object)
            new("QuadTypeIndependenceTest", statistic = object,
                distribution = distribution(object))
        })

    object@call <- match.call()
    ## return object inheriting from class `IndependenceTest'
    return(object)
}
