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
