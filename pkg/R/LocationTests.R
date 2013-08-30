### permutation test without transformations
oneway_test <- function(object, ...) UseMethod("oneway_test")

oneway_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("oneway_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

oneway_test.IndependenceProblem <- function(object, ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"),
                 " does not represent a K-sample problem",
                 " (maybe the grouping variable is not a factor?)")
        return(TRUE)
    }

    twosamp <- nlevels(object@x[[1]]) == 2

    args <- setup_args(check = check)

    RET <- do.call("independence_test", c(list(object = object), args))

    if (is_singly_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        RET@method <- "2-Sample Permutation Test"
        RET@nullvalue <- 0
    } else
        RET@method <- "K-Sample Permutation Test"

    return(RET)
}


### OK, OK, here is the most prominent one ...
wilcox_test <- function(object, ...) UseMethod("wilcox_test")

wilcox_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("wilcox_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

wilcox_test.IndependenceProblem <- function(object,
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_2sample(object) && is_numeric_y(object)))
            stop(sQuote("object"),
                 " does not represent a two-sample problem",
                 " (maybe the grouping variable is not a factor?)")
        return(TRUE)
    }

    args <- setup_args(teststat = "scalar",
                       ytrafo = function(data)
                           trafo(data, numeric_trafo = rank_trafo),
                       check = check)

    RET <- do.call("independence_test", c(list(object = object), args))

    RET@method <- "Wilcoxon-Mann-Whitney Test"
    RET@nullvalue <- 0

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_location(RET@statistic, RET@distribution,
                             level = level)
        RET@conf.level <- conf.level
    }

    return(RET)
}


### Kruskal-Wallis test
kruskal_test <- function(object, ...) UseMethod("kruskal_test")

kruskal_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("kruskal_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

kruskal_test.IndependenceProblem <- function(object,
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
                           trafo(data, numeric_trafo = rank_trafo),
                       check = check)
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <- if (is.ordered(object@x[[1]])) "scalar"
                     else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (is_singly_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else
        RET@method <- "Kruskal-Wallis Test"

    return(RET)
}


### normal quantiles (van der Waerden) test
normal_test <- function(object, ...) UseMethod("normal_test")

normal_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("normal_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

normal_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"),
                 " does not represent a K-sample problem",
                 " (maybe the grouping variable is not a factor?)")
        return(TRUE)
    }

    twosamp <- nlevels(object@x[[1]]) == 2

    args <- setup_args(ytrafo = function(data)
                           trafo(data, numeric_trafo = function(y)
                               normal_trafo(y, ties.method = ties.method)),
                       check = check)
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <- if (is.ordered(object@x[[1]]) || twosamp) "scalar"
                     else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (is_singly_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        RET@method <- "2-Sample Normal Quantile (van der Waerden) Test"
        RET@nullvalue <- 0
        if (conf.int) {
            RET <- new("ScalarIndependenceTestConfint", RET)
            RET@confint <- function(level)
                confint_location(RET@statistic, RET@distribution,
                                 level = level)
            RET@conf.level <- conf.level
        }
    } else
        RET@method <- "K-Sample Normal Quantile (van der Waerden) Test"

    return(RET)
}


### median test
median_test <- function(object, ...) UseMethod("median_test")

median_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("median_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

median_test.IndependenceProblem <- function(object,
    mid.score = c("0", "0.5", "1"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"),
                 " does not represent a two-sample problem",
                 " (maybe the grouping variable is not a factor?)")
        return(TRUE)
    }

    twosamp <- nlevels(object@x[[1]]) == 2

    args <- setup_args(ytrafo = function(data)
                           trafo(data, numeric_trafo = function(y)
                               median_trafo(y, mid.score = mid.score)),
                       check = check)
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <- if (is.ordered(object@x[[1]]) || twosamp) "scalar"
                     else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (is_singly_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        RET@method <- "2-Sample Median Test"
        RET@nullvalue <- 0
        if (conf.int) {
            RET <- new("ScalarIndependenceTestConfint", RET)
            RET@confint <- function(level)
                confint_location(RET@statistic, RET@distribution,
                                 level = level)
            RET@conf.level <- conf.level
        }
    } else
        RET@method <- "K-Sample Median Test"

    return(RET)
}


### Savage test
savage_test <- function(object, ...) UseMethod("savage_test")

savage_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("savage_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

savage_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"),
                 " does not represent a K-sample problem",
                 " (maybe the grouping variable is not a factor?)")
        return(TRUE)
    }

    twosamp <- nlevels(object@x[[1]]) == 2

    args <- setup_args(ytrafo = function(data)
                           trafo(data, numeric_trafo = function(y)
                               savage_trafo(y, ties.method = ties.method)),
                       check = check)
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <- if (is.ordered(object@x[[1]]) || twosamp) "scalar"
                     else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (is_singly_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        RET@method <- "2-Sample Savage Test"
        RET@nullvalue <- 0
        if (conf.int) {
            RET <- new("ScalarIndependenceTestConfint", RET)
            RET@confint <- function(level)
                confint_location(RET@statistic, RET@distribution,
                                 level = level)
            RET@conf.level <- conf.level
        }
    } else
        RET@method <- "K-Sample Savage Test"

    return(RET)
}
