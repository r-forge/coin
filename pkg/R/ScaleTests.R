### Ansari-Bradley test
ansari_test <- function(object, ...) UseMethod("ansari_test")

ansari_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("ansari_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

ansari_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"),
                 " does not represent a K-sample problem",
                 " (maybe the grouping variable is not a factor?)")
        if (is_singly_ordered(object))
            stop(colnames(object@x), " is an ordered factor")
        return(TRUE)
    }

    twosamp <- nlevels(object@x[[1]]) == 2

    args <- setup_args(ytrafo = function(data)
                           trafo(data, numeric_trafo = function(y)
                               ansari_trafo(y, ties.method = ties.method)),
                       check = check)
    args$teststat <- if (twosamp) "scalar" else "quad"
    args$alternative <- match.arg(args$alternative,
                                  c("two.sided", "less", "greater"))
    if (args$alternative == "less")
        args$alternative <- "greater"
    else if (args$alternative == "greater")
        args$alternative <- "less"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        RET@method <- "2-Sample Ansari-Bradley Test"
        RET@parameter <- "ratio of scales"
        RET@nullvalue <- 1
        if (conf.int) {
            RET <- new("ScalarIndependenceTestConfint", RET)
            RET@confint <- function(level)
                confint_scale(RET@statistic, RET@distribution,
                                 level = level)
            RET@conf.level <- conf.level
        }
    } else
        RET@method <- "K-Sample Ansari-Bradley Test"

    return(RET)
}


### Fligner-Killeen test
fligner_test <- function(object, ...) UseMethod("fligner_test")

fligner_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("fligner_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

fligner_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(object) {
        if (!(is_Ksample(object) && is_numeric_y(object)))
            stop(sQuote("object"),
                 " does not represent a K-sample problem",
                 " (maybe the grouping variable is not a factor?)")
        if (is_singly_ordered(object))
            stop(colnames(object@x), " is an ordered factor")
        return(TRUE)
    }

    twosamp <- nlevels(object@x[[1]]) == 2

    args <- setup_args(ytrafo = function(data)
                           trafo(data, numeric_trafo = function(y)
                               fligner_trafo(y, ties.method = ties.method)),
                       check = check)
    args$teststat <- if (twosamp) "scalar" else "quad"

    ## eliminate location differences (see 'stats/R/fligner.test')
    object@y[[1]] <- object@y[[1]] -
        tapply(object@y[[1]], object@x[[1]], median)[object@x[[1]]]

    RET <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        RET@method <- "2-Sample Fligner-Killeen Test"
        RET@parameter <- "ratio of scales"
        RET@nullvalue <- 1
        if (conf.int) {
            RET <- new("ScalarIndependenceTestConfint", RET)
            RET@confint <- function(level)
                confint_scale(RET@statistic, RET@distribution,
                                 level = level)
            RET@conf.level <- conf.level
        }
    } else
        RET@method <- "K-Sample Fligner-Killeen Test"

    return(RET)
}
