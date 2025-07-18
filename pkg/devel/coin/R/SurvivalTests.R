### weighted logrank test
logrank_test <- function(object, ...) UseMethod("logrank_test")

logrank_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("logrank_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

logrank_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "Hothorn-Lausen", "average-scores"),
    type = c("logrank", "Gehan-Breslow", "Tarone-Ware", "Peto-Peto", "Prentice",
             "Prentice-Marek", "Andersen-Borgan-Gill-Keiding",
             "Fleming-Harrington", "Gaugler-Kim-Liao", "Self"),
    rho = NULL, gamma = NULL, ...) {

    type <- match.arg(type)

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, surv_trafo = function(y)
                logrank_trafo(y, ties.method = ties.method,
                              type = type, rho = rho,
                              gamma = gamma)),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_censored_y(object))
                stop(sQuote(colnames(object@y)), " is not a censored variable")
            if (!is_unity(object@weights))
                warning("rank transformation doesn't take weights into account")
            TRUE
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear tests
    args$teststat <- if (is_ordered_x(object) || twosamp) "scalar"
                     else "quadratic"

    object <- do.call(independence_test, c(object = object, args))

    if (is_ordered_x(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        object@method <- paste("Two-Sample",
                            if (type == "logrank") "Logrank" else type, "Test")
        object@parameter <- "theta" # theta = lambda_2 / lambda_1
        object@nullvalue <- 1       # theta = 1 => Equal hazards
    } else
        object@method <- paste("K-Sample",
                            if (type == "logrank") "Logrank" else type, "Test")

    object
}
