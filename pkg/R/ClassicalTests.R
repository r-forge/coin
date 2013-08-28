
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


### weighted logrank test
surv_test <- function(object, ...) UseMethod("surv_test")

surv_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("surv_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

surv_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "Hothorn-Lausen", "average-scores"),
    type = c("logrank", "Gehan-Breslow", "Tarone-Ware", "Prentice",
             "Prentice-Marek", "Andersen-Borgan-Gill-Keiding",
             "Fleming-Harrington", "Self"),
    rho = NULL, gamma = NULL, ...) {

    type <- match.arg(type)[1]

    check <- function(object) {
        if (!(is_Ksample(object) && is_censored_y(object)))
            stop(sQuote("object"),
                 " does not represent a K-sample problem with censored data",
                 " (maybe the grouping variable is not a factor?)")
        return(TRUE)
    }

    twosamp <- nlevels(object@x[[1]]) == 2

    args <- setup_args(ytrafo = function(data)
                           trafo(data, surv_trafo = function(y)
                               logrank_trafo(y, ties.method = ties.method,
                                             type = type, rho = rho,
                                             gamma = gamma)),
                       check = check)
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <- if (is.ordered(object@x[[1]]) || twosamp) "scalar"
                     else "quad"

    RET <- do.call("independence_test", c(list(object = object), args))

    if (is_singly_ordered(RET@statistic))
        RET@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        RET@method <- paste("2-Sample",
                            if (type == "logrank") "Logrank" else type, "Test")
        ## theta = lambda_2 / lambda_1
        RET@parameter <- "theta"
        RET@nullvalue <- 1 # theta = 1 => Equal hazards
    } else
        RET@method <- paste("K-Sample",
                            if (type == "logrank") "Logrank" else type, "Test")

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


### contrast test
contrast_test <- function(object, ...) UseMethod("contrast_test")

contrast_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("contrast_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

contrast_test.IndependenceProblem <- function(object,
    cmatrix, distribution = c("asymptotic", "approximate"), ...) {

    if (!(ncol(object@x) == 1 && is.factor(object@x[[1]])))
        stop(sQuote("object@x"), " is not univariate or a factor")

    if  (!is.matrix(cmatrix) || nrow(cmatrix) != nlevels(object@x[[1]]))
        stop(sQuote("cmatrix"), " is not a matrix with ",
             nlevels(object@x), " rows")

    if (is.null(colnames(cmatrix)))
        colnames(cmatrix) <- paste0("C", 1:ncol(cmatrix))

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(teststat = "max",
                       distribution = distribution,
                       xtrafo = function(data)
                           trafo(data) %*% cmatrix)

    RET <- do.call("independence_test", c(list(object = object), args))

    RET@method <- "General Contrast Test"

    return(RET)
}


### maxstat test
maxstat_test <- function(object, ...) UseMethod("maxstat_test")

maxstat_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("maxstat_test", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

maxstat_test.table <- function(object, ...) {

    ip <- table2IndependenceProblem(object)
    RET <- do.call("maxstat_test", c(list(object = ip), list(...)))
    return(RET)
}

maxstat_test.IndependenceProblem <- function(object,
    teststat = c("max", "quad"),
    distribution = c("asymptotic", "approximate"),
    minprob = 0.1, maxprob = 1 - minprob, ...) {

    teststat <- match.arg(teststat)
    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    xtrafo <- function(data)
        trafo(data,
              numeric_trafo = function(x)
                  maxstat_trafo(x, minprob = minprob, maxprob = maxprob),
              factor_trafo = function(x)
                  fmaxstat_trafo(x, minprob = minprob, maxprob = maxprob),
              ordered_trafo = function(x)
                  ofmaxstat_trafo(x, minprob = minprob, maxprob = maxprob))

    args <- setup_args(teststat = teststat,
                       distribution = distribution,
                       xtrafo = xtrafo)

    ## convert factors to ordered and attach scores if requested
    object <- setscores(object, args$scores)
    args$scores <- NULL

    RET <- do.call("independence_test", c(list(object = object), args))

    RET@method <- "Generalized Maximally Selected Statistics"

    ## estimate cutpoint
    wm <- which.max(apply(abs(statistic(RET, "standardized")), 1, max))
    whichvar <- attr(RET@statistic@xtrans, "assign")[wm]
    maxcontr <- RET@statistic@xtrans[, wm]
    if (is.factor(RET@statistic@x[[whichvar]])) {
        cp <- levels(RET@statistic@x[[whichvar]][maxcontr > 0][, drop = TRUE])
        cp0 <- levels(RET@statistic@x[[whichvar]][maxcontr == 0][, drop = TRUE])
        lab <- paste0("{", paste0(cp, collapse = ", "), "} vs. {",
                      paste0(cp0, collapse = ", "), "}")
    } else {
        cp <- max(RET@statistic@x[[whichvar]][maxcontr > 0])
        lab <- paste0("<= ", format(cp, digits = getOption("digits")))
    }
    if (ncol(object@x) > 1)
        estimate <- list(covariable = colnames(RET@statistic@x)[whichvar],
                         cutpoint = cp, label = lab)
    else
        estimate <- list(cutpoint = cp, label = lab)
    class(estimate) <- c("list", "cutpoint")
    RET@estimates <- list(estimate = estimate)

    return(RET)
}

### and now symmetry tests

### Friedman Test
friedman_test <- function(object, ...) UseMethod("friedman_test")

friedman_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("friedman_test", c(list(object = sp), list(...)))
    return(RET)
}

friedman_test.SymmetryProblem <- function(object,
    distribution = c("asymptotic", "approximate"), ...) {

    if (!is_completeblock(object))
        stop("Not an unreplicated complete block design")

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(distribution = distribution,
                       ytrafo = function(data)
                           trafo(data, numeric_trafo = rank_trafo,
                                 block = object@block))
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <- if (is.ordered(object@x[[1]])) "scalar"
                     else "quad"

    RET <- do.call("symmetry_test", c(list(object = object), args))

    if (is_singly_ordered(RET@statistic))
        RET@method <- "Page Test"
    else
        RET@method <- "Friedman Test"

    return(RET)
}


### marginal homogeneity test
mh_test <- function(object, ...) UseMethod("mh_test")

mh_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("mh_test", c(list(object = sp), list(...)))
    return(RET)
}

mh_test.table <- function(object, ...) {
    df <- table2df_sym(object)
    sp <- new("SymmetryProblem", x = df["conditions"], y = df["response"])
    RET <- do.call("mh_test", c(list(object = sp), list(...)))
    return(RET)
}

mh_test.SymmetryProblem <- function(object,
    distribution = c("asymptotic", "approximate"), ...) {

    if (!is_completeblock(object))
        stop("Not an unreplicated complete block design")
    if (ncol(object@y) != 1 || !is.factor(object@y[[1]]))
        stop("Response variable is not a factor")

    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))

    args <- setup_args(distribution = distribution)
    object <- setscores(object, args$scores)
    args$scores <- NULL
    args$teststat <-
        if ((is.ordered(object@x[[1]]) && is.ordered(object@y[[1]])) ||
                ((is.ordered(object@x[[1]]) && nlevels(object@y[[1]]) == 2) ||
                 (is.ordered(object@y[[1]]) && nlevels(object@x[[1]]) == 2)))
            "scalar"
        else "quad"

    RET <- do.call("symmetry_test", c(list(object = object), args))

    if (is_ordered(RET@statistic))
        RET@method <- "Marginal Homogeneity Test for Ordered Data"
    else
        RET@method <- "Marginal Homogeneity Test"

    return(RET)
}


### Wilcoxon signed-rank test
wilcoxsign_test <- function(object, ...) UseMethod("wilcoxsign_test")

wilcoxsign_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    if (is.null(d$bl))
        d <- list(y = data.frame(c(d$y[[1]], d$x[[1]])),
                  x = data.frame(gl(2, length(d$x[[1]]))),
                  block = factor(rep(1:length(d$x[[1]]), 2)))
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("wilcoxsign_test", c(list(object = sp), list(...)))
    return(RET)
}

wilcoxsign_test.SymmetryProblem <- function(object,
    zero.method = c("Pratt", "Wilcoxon"), ...) {

    zero.method <- match.arg(zero.method)

    y <- object@y[[1]]
    x <- object@x[[1]]
    block <- object@block

    if (!is.numeric(y))
        stop(sQuote("y"), " is not a numeric variable")
    if (is.factor(x)) {
        if (nlevels(x) != 2)
            stop(sQuote("x"), " is not a factor with two levels")
        diffs <- tapply(1:length(y), block, function(b)
            y[b][x[b] == levels(x)[1]] - y[b][x[b] == levels(x)[2]]
        )
    } else {
        stop(sQuote("x"), " is not a factor")
    }

    abs_diffs <- abs(diffs)
    if (all(abs_diffs < .Machine$double.eps))
        stop("all pairwise differences equal zero")

    pos_abs_diffs <- abs_diffs > 0
    if (zero.method == "Pratt") {
        rank_abs_diffs <- rank(abs_diffs)
        pos <- (rank_abs_diffs * (diffs > 0))[pos_abs_diffs]
        neg <- (rank_abs_diffs * (diffs < 0))[pos_abs_diffs]
    } else {
        diffs <- diffs[pos_abs_diffs]
        abs_diffs <- abs_diffs[pos_abs_diffs]
        rank_abs_diffs <- rank(abs_diffs)
        pos <- rank_abs_diffs * (diffs > 0)
        neg <- rank_abs_diffs * (diffs < 0)
    }
    n <- length(pos)

    y <- as.vector(rbind(pos, neg))
    x <- factor(rep(0:1, n), labels = c("pos", "neg"))
    block <- gl(n, 2)

    ip <- new("IndependenceProblem", x = data.frame(x = x),
              y = data.frame(y = y), block = block)

    args <- setup_args(teststat = "scalar", paired = TRUE)

    RET <- do.call("independence_test", c(list(object = ip), args))

    if (zero.method == "Pratt")
        RET@method <- "Wilcoxon-Pratt Signed-Rank Test"
    else
        RET@method <- "Wilcoxon Signed-Rank Test"
    RET@nullvalue <- 0

    return(RET)
}


### Sign test
sign_test <- function(object, ...) UseMethod("sign_test")

sign_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    if (is.null(d$bl))
        d <- list(y = data.frame(c(d$y[[1]], d$x[[1]])),
                  x = data.frame(gl(2, length(d$x[[1]]))),
                  block = factor(rep(1:length(d$x[[1]]), 2)))
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("sign_test", c(list(object = sp), list(...)))
    return(RET)
}

sign_test.SymmetryProblem <- function(object, ...) {

    y <- object@y[[1]]
    x <- object@x[[1]]
    block <- object@block

    if (!is.numeric(y))
        stop(sQuote("y"), " is not a numeric variable")
    if (is.factor(x)) {
        if (nlevels(x) != 2)
            stop(sQuote("x"), " is not a factor with two levels")
        diffs <- tapply(1:length(y), block, function(b)
            y[b][x[b] == levels(x)[1]] - y[b][x[b] == levels(x)[2]]
        )
    } else {
        stop(sQuote("x"), " is not a factor")
    }

    abs_diffs <- abs(diffs)
    if (all(abs_diffs < .Machine$double.eps))
        stop("all pairwise differences equal zero")

    diffs <- diffs[abs_diffs > 0]
    n <- length(diffs)

    y <- as.vector(rbind(as.numeric(diffs > 0), as.numeric(diffs < 0)))
    x <- factor(rep(0:1, n), labels = c("pos", "neg"))
    block <- gl(n, 2)

    ip <- new("IndependenceProblem", x = data.frame(x = x),
              y = data.frame(y = y), block = block)

    args <- setup_args(teststat = "scalar", paired = TRUE)

    RET <- do.call("independence_test", c(list(object = ip), args))

    RET@method <- "Sign Test"
    RET@nullvalue <- 0

    return(RET)
}
