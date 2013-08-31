### Pearson's chi-squared test
chisq_test <- function(object, ...) UseMethod("chisq_test")

chisq_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("chisq_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

chisq_test.table <- function(object, ...) {

    ip <- table2IndependenceProblem(object)
    RET <- do.call("chisq_test", c(list(object = ip), list(...)))
    return(RET)
}

chisq_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!is_contingency(object))
            stop(sQuote("object"),
                 " does not represent a contingency problem")
        if (nlevels(object@block) != 1)
            stop(sQuote("object"), " contains blocks: use ",
                 sQuote("cmh_test"), " instead")
        return(TRUE)
    }
    n <- sum(object@weights)

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args()
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <-
        if ((is.ordered(object@x[[1]]) && is.ordered(object@y[[1]])) ||
                ((is.ordered(object@x[[1]]) && nlevels(object@y[[1]]) == 2) ||
                 (is.ordered(object@y[[1]]) && nlevels(object@x[[1]]) == 2)))
            "scalar"
        else "quad"
    args$alternative <- match.arg(args$alternative,
                                  c("two.sided", "less", "greater"))

    ## transform data if requested and setup a test problem
    itp <- new("IndependenceTestProblem", object)

    if (!check(itp))
        stop(sQuote("check"), " failed")

    its <- new("IndependenceTestStatistic", itp, varonly = FALSE)

    ## use the classical chisq statistic based on Pearson
    ## residuals (O - E)^2 / E
    ## see Th. 3.1 and its proof in Strasser & Weber (1999).
    RET <-
        if (args$teststat == "scalar") {
            ts <- new("ScalarIndependenceTestStatistic", its, args$alternative)
            ts@teststatistic <- ts@teststatistic * sqrt(n / (n - 1))
            ts@standardizedlinearstatistic <-
                ts@standardizedlinearstatistic * sqrt(n / (n - 1))
            ts@covariance <- new("CovarianceMatrix", covariance(ts) * (n - 1) / n)
            new("ScalarIndependenceTest", statistic = ts,
                distribution = distribution(ts))
        } else {
            if (args$alternative != "two.sided")
                warning(sQuote("alternative"),
                        " is ignored for quad type test statistics")
            ts <- new("QuadTypeIndependenceTestStatistic", its)
            ts@teststatistic <- ts@teststatistic * n / (n - 1)
            ts@standardizedlinearstatistic <-
                ts@standardizedlinearstatistic * sqrt(n / (n - 1))
            ts@covariance <- new("CovarianceMatrix", covariance(ts) * (n - 1) / n)
            ts@covarianceplus <- MPinv(covariance(ts))$MPinv
            new("QuadTypeIndependenceTest", statistic = ts,
                distribution = distribution(ts))
        }

    if (is_doubly_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else if (is_singly_ordered(RET@statistic))
        RET@method <- "Generalized Pearson's Chi-Squared Test"
    else
        RET@method <- "Pearson's Chi-Squared Test"

    RET@call <- match.call()

    return(RET)
}


### generalized Cochran-Mantel-Haenzel test
cmh_test <- function(object, ...) UseMethod("cmh_test")

cmh_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("cmh_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

cmh_test.table <- function(object, ...) {

    ip <- table2IndependenceProblem(object)
    RET <- do.call("cmh_test", c(list(object = ip), list(...)))
    return(RET)
}

cmh_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!is_contingency(object))
            stop(sQuote("object"),
                 " does not represent a contingency problem")
        return(TRUE)
    }

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(distribution = distribution,
                       check = check)
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <-
        if ((is.ordered(object@x[[1]]) && is.ordered(object@y[[1]])) ||
                ((is.ordered(object@x[[1]]) && nlevels(object@y[[1]]) == 2) ||
                 (is.ordered(object@y[[1]]) && nlevels(object@x[[1]]) == 2)))
            "scalar"
        else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (is_doubly_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else
        RET@method <- "Generalized Cochran-Mantel-Haenszel Test"

    return(RET)
}


### linear-by-linear association test
lbl_test <- function(object, ...) UseMethod("lbl_test")

lbl_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("lbl_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

lbl_test.table <- function(object, ...) {

    ip <- table2IndependenceProblem(object)
    RET <- do.call("lbl_test", c(list(object = ip), list(...)))
    return(RET)
}

lbl_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!is_doubly_ordered(object))
            stop(sQuote("object"),
                 " does not represent a problem with ordered data")
        return(TRUE)
    }

    ## convert factors to ordered
    object@x[] <- lapply(object@x,
        function(x) if (is.factor(x) && nlevels(x) > 2) {
                        return(ordered(x))
                    } else {
                        return(x)
                    })
    object@y[] <- lapply(object@y,
        function(x) if (is.factor(x) && nlevels(x) > 2) {
                        return(ordered(x))
                    } else {
                        return(x)
                    })

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(teststat = "scalar",
                       distribution = distribution,
                       check = check)

    RET <- do.call("independence_test", c(list(object = object), args))

    RET@method <- "Linear-by-Linear Association Test"

    return(RET)
}
