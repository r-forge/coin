### Taha test
taha_test <- function(object, ...) UseMethod("taha_test")

taha_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("taha_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

taha_test.IndependenceProblem <- function(object,
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
                               rank_trafo(y)^2),
                       check = check)
    args$teststat <- if (twosamp) "scalar" else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        RET@method <- "2-Sample Taha Test"
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
        RET@method <- "K-Sample Taha Test"

    return(RET)
}


### Klotz Test
klotz_test <- function(object, ...) UseMethod("klotz_test")

klotz_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("klotz_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

klotz_test.IndependenceProblem <- function(object,
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
                               klotz_trafo(y, ties.method = ties.method)),
                       check = check)
    args$teststat <- if (twosamp) "scalar" else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        RET@method <- "2-Sample Klotz Test"
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
        RET@method <- "K-Sample Klotz Test"

    return(RET)
}


### Mood Test
mood_test <- function(object, ...) UseMethod("mood_test")

mood_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("mood_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

mood_test.IndependenceProblem <- function(object,
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
                               mood_trafo(y, ties.method = ties.method)),
                       check = check)
    args$teststat <- if (twosamp) "scalar" else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        RET@method <- "2-Sample Mood Test"
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
        RET@method <- "K-Sample Mood Test"

    return(RET)
}


### Ansari-Bradley test
ansari_test <- function(object, ...) UseMethod("ansari_test")

ansari_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("ansari_test", "IndependenceProblem", formula, data, subset, weights,
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
    alternative <- match.arg(args$alternative,
                             c("two.sided", "less", "greater"))
    if (alternative == "less")
        args$alternative <- "greater"
    else if (alternative == "greater")
        args$alternative <- "less"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        RET@method <- "2-Sample Ansari-Bradley Test"
        RET@parameter <- "ratio of scales"
        RET@nullvalue <- 1
        RET@statistic@alternative <- alternative
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

    ft("fligner_test", "IndependenceProblem", formula, data, subset, weights,
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


### Conover-Iman (1978) test
conover_test <- function(object, ...) UseMethod("conover_test")

conover_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("conover_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

conover_test.IndependenceProblem <- function(object,
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
                               rank_trafo(abs(y))^2),
                       check = check)
    args$teststat <- if (twosamp) "scalar" else "quad"

    ## eliminate location differences
    object@y[[1]] <- object@y[[1]] -
        tapply(object@y[[1]], object@x[[1]], mean)[object@x[[1]]]

    RET <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        RET@method <- "2-Sample Conover-Iman Test"
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
        RET@method <- "K-Sample Conover-Iman Test"

    return(RET)
}
