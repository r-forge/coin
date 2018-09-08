setMethod("show",
    signature = "IndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )

        RET <- list(
            statistic = setNames(statistic(object), nm = "c"),
            p.value = pvalue(object),
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "MaxTypeIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )

        RET <- list(
            statistic = setNames(statistic(object), nm = "maxT"),
            p.value = pvalue(object),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        estimates <- object@estimates
        if (length(estimates) > 0)
            RET <- c(RET, estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "QuadTypeIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )

        RET <- list(
            statistic = setNames(statistic(object), nm = "chi-squared"),
            p.value = pvalue(object),
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        parameters <- object@distribution@parameters
        if (length(parameters) == 1 && names(parameters) == "df")
            RET$parameter <- setNames(parameters[[1]], nm = "df")
        estimates <- object@estimates
        if (length(estimates) > 0)
            RET <- c(RET, estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "ScalarIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )

        RET <- list(
            statistic = setNames(statistic(object), nm = "Z"),
            p.value = pvalue(object),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        nullvalue <- object@nullvalue
        if (length(nullvalue) > 0)
            RET$null.value <- setNames(nullvalue, nm = object@parameter)
        estimates <- object@estimates
        if (length(estimates) > 0)
            RET <- c(RET, estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "ScalarIndependenceTestConfint",
    definition = function(object) {
        distname <- switch(class(object@distribution),
                        "AsymptNullDistribution" = "Asymptotic",
                        "ApproxNullDistribution" = "Approximative",
                        "ExactNullDistribution"  = "Exact"
                    )
        ci <- confint(object, level = object@conf.level)

        RET <- list(
            statistic = setNames(statistic(object), nm = "Z"),
            p.value = pvalue(object),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method),
            conf.int = ci$conf.int,
            estimate = ci$estimate
        )
        nullvalue <- object@nullvalue
        if (length(nullvalue) > 0)
            RET$null.value <- setNames(nullvalue, nm = object@parameter)
        estimates <- object@estimates
        if (length(estimates) > 0)
            RET <- c(RET, estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

print.ci <- function(x, ...) {
    if (hasName(x, "conf.int")) {
        cat(format(100 * attr(x$conf.int, "conf.level")),
            "percent confidence interval:\n",
            format(c(x$conf.int[1], x$conf.int[2]), ...), "\n")
    }
    if (hasName(x, "estimate")) {
        cat("sample estimates:\n")
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
}

print.MCp <- function(x, ...) {
    p <- x
    attributes(p) <- NULL
    print(p, ...)
    ci <- list(conf.int = attr(x, "conf.int"))
    class(ci) <- "ci"
    print(ci, ...)
}

print.cutpoint <- function(x, ...) {
    cat(paste0("  ", dQuote("best"), " cutpoint: ", x$label, "\n"))
    if (hasName(x, "covariable"))
        cat(paste0("       covariable: ", x$covariable, "\n"))
}
