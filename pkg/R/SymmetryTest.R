### a generic test procedure for classical (and not so classical) tests
symmetry_test <- function(object, ...) UseMethod("symmetry_test")

symmetry_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    ft("symmetry_test", "SymmetryProblem", formula, data, subset,
       frame = parent.frame(), ...)
}

symmetry_test.table <- function(object, ...) {

    object <- table2df_sym(object)
    object <- new("SymmetryProblem", x = object["conditions"],
                  y = object["response"])
    object <- do.call("symmetry_test", c(list(object = object), list(...)))
    return(object)
}

symmetry_test.SymmetryProblem <- function(object,
    teststat = c("maximum", "quadratic", "scalar"),
    distribution = c("asymptotic", "approximate", "exact"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, scores = NULL,
    check = NULL, ...) {

    args <- setup_args(distribution = check_distribution_arg(distribution,
                           match.arg(distribution)))

    class(object) <- "IndependenceProblem"
    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "General Symmetry Test"

    return(object)
}
