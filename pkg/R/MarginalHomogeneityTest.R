### marginal homogeneity test
mh_test <- function(object, ...) UseMethod("mh_test")

mh_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    ft("mh_test", "SymmetryProblem", formula, data, subset,
       frame = parent.frame(), ...)
}

mh_test.table <- function(object, ...) {

    object <- table2df_sym(object)
    object <- new("SymmetryProblem", x = object["conditions"],
                  y = object["response"])
    object <- do.call("mh_test", c(list(object = object), list(...)))
    return(object)
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
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    args$teststat <-
        if ((is.ordered(object@x[[1]]) && is.ordered(object@y[[1]])) ||
                ((is.ordered(object@x[[1]]) && nlevels(object@y[[1]]) == 2) ||
                 (is.ordered(object@y[[1]]) && nlevels(object@x[[1]]) == 2)))
            "scalar"
        else "quad"

    object <- do.call("symmetry_test", c(list(object = object), args))

    if (is_ordered(object@statistic))
        object@method <- "Marginal Homogeneity Test for Ordered Data"
    else
        object@method <- "Marginal Homogeneity Test"

    return(object)
}
