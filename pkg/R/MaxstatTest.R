### generalized maximally selected statistics
maxstat_test <- function(object, ...) UseMethod("maxstat_test")

maxstat_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("maxstat_test", "IndependenceProblem", formula, data, subset, weights,
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
