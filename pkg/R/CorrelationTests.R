### Spearman test
spearman_test <- function(object, ...) UseMethod("spearman_test")

spearman_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("spearman_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

spearman_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate"), ...) {

    check <- function(object) {
        if (!is_corr(object))
            stop(sQuote("object"),
                 " does not represent a univariate correlation problem")
        return(TRUE)
    }

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(teststat = "scalar",
                       distribution = distribution,
                       xtrafo = function(data)
                           trafo(data, numeric_trafo = rank_trafo),
                       ytrafo = function(data)
                           trafo(data, numeric_trafo = rank_trafo),
                       check = check)

    RET <- do.call("independence_test", c(list(object = object), args))

    RET@method <- "Spearman Correlation Test"
    RET@parameter <- "rho"
    RET@nullvalue <- 0

    return(RET)
}


### Fisher-Yates correlation test (based on van der Waerden scores)
fisyat_test <- function(object, ...) UseMethod("fisyat_test")

fisyat_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("fisyat_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

fisyat_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate"),
    ties.method = c("mid-ranks", "average-scores"), ...) {

    check <- function(object) {
        if (!is_corr(object))
            stop(sQuote("object"),
                 " does not represent a univariate correlation problem")
        return(TRUE)
    }

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(teststat = "scalar",
                       distribution = distribution,
                       xtrafo = function(data)
                           trafo(data, numeric_trafo = function(x)
                               normal_trafo(x, ties.method = ties.method)),
                       ytrafo = function(data)
                           trafo(data, numeric_trafo = function(y)
                               normal_trafo(y, ties.method = ties.method)),
                       check = check)

    RET <- do.call("independence_test", c(list(object = object), args))

    RET@method <- "Fisher-Yates Correlation Test"
    RET@parameter <- "rho"
    RET@nullvalue <- 0

    return(RET)
}


## quadrant test (Hajek, Sidak & Sen, pp. 124--125)
quadrant_test <- function(object, ...) UseMethod("quadrant_test")

quadrant_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("quadrant_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

quadrant_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate"),
    mid.score = c("0", "0.5", "1"), ...) {
    ## <FIXME> in principle is "exact" also possible, unless mid.score == "0.5",
    ## since the data is effectively reduced to a 2x2 table.  But...
    ## </FIXME>

    check <- function(object) {
        if (!is_corr(object))
            stop(sQuote("object"),
                 " does not represent a univariate correlation problem")
        return(TRUE)
    }

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(teststat = "scalar",
                       distribution = distribution,
                       xtrafo = function(data)
                           trafo(data, numeric_trafo = function(x)
                               median_trafo(x, mid.score = mid.score)),
                       ytrafo = function(data)
                           trafo(data, numeric_trafo = function(y)
                               median_trafo(y, mid.score = mid.score)),
                       check = check)

    RET <- do.call("independence_test", c(list(object = object), args))

    RET@method <- "Quadrant Test"
    RET@parameter <- "rho"
    RET@nullvalue <- 0

    return(RET)
}


## Koziol-Nemec test
koziol_test <- function(object, ...) UseMethod("koziol_test")

koziol_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("koziol_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

koziol_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate"),
    ties.method = c("mid-ranks", "average-scores"), ...) {

    check <- function(object) {
        if (!is_corr(object))
            stop(sQuote("object"),
                 " does not represent a univariate correlation problem")
        return(TRUE)
    }

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(teststat = "scalar",
                       distribution = distribution,
                       xtrafo = function(data)
                           trafo(data, numeric_trafo = function(x)
                               koziol_trafo(x, ties.method = ties.method)),
                       ytrafo = function(data)
                           trafo(data, numeric_trafo = function(y)
                               koziol_trafo(y, ties.method = ties.method)),
                       check = check)

    RET <- do.call("independence_test", c(list(object = object), args))

    RET@method <- "Koziol-Nemec Test"
    RET@parameter <- "rho"
    RET@nullvalue <- 0

    return(RET)
}
