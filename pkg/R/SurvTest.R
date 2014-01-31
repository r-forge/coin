### weighted logrank test
surv_test <- function(object, ...) UseMethod("surv_test")

surv_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("surv_test", "IndependenceProblem", formula, data, subset, weights,
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
