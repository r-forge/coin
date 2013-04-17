### a generic test procedure for classical (and not so classical) tests
symmetry_test <- function(object, ...) UseMethod("symmetry_test")

symmetry_test.formula <- function(formula, data = list(), subset = NULL,
    ...) {

    d <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    sp <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("symmetry_test", c(list(object = sp), list(...)))
    return(RET)
}

symmetry_test.table <- function(object, ...) {
    df <- table2df_sym(object)
    sp <- new("SymmetryProblem", x = df["conditions"], y = df["response"])
    RET <- do.call("symmetry_test", c(list(object = sp), list(...)))
    return(RET)
}

symmetry_test.SymmetryProblem <- function(object,
    teststat = c("max", "quad", "scalar"),
    distribution = c("asymptotic", "approximate"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, scores = NULL,
    check = NULL, ...) {
    distribution <- check_distribution_arg(distribution,
        values = c("asymptotic", "approximate"))
    class(object) <- "IndependenceProblem"
    RET <- independence_test(object, teststat, distribution, alternative,
                             xtrafo, ytrafo, scores, check, ...)
    RET@method <- "General Symmetry Test"
    return(RET)
}
