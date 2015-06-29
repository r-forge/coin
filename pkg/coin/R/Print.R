setMethod("show",
    signature = "IndependenceTest",
    definition = function(object) {
        stat <- object@statistic@teststatistic
        names(stat) <- "c"

        dataname <- varnames(object@statistic)

        dist <- object@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(statistic = stat,
                    p.value = object@distribution@pvalue(stat),
                    data.name = dataname,
                    method = paste(distname, object@method))
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "MaxTypeIndependenceTest",
    definition = function(object) {
        stat <- object@statistic@teststatistic
        names(stat) <- "maxT"
        dist <- object@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        dataname <- varnames(object@statistic)

        RET <- list(statistic = stat,
                    p.value = object@distribution@pvalue(stat),
                    alternative = object@statistic@alternative,
                    data.name = dataname,
                    method = paste(distname, object@method))
        if (length(object@estimates) > 0)
            RET <- c(RET, object@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "QuadTypeIndependenceTest",
    definition = function(object) {
        stat <- object@statistic@teststatistic
        names(stat) <- "chi-squared"
        dist <- object@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        dataname <- varnames(object@statistic)

        RET <- list(statistic = stat,
                    p.value = dist@pvalue(stat),
                    data.name = dataname,
                    method = paste(distname, object@method))
        if (length(dist@parameters) == 1 &&
                   names(dist@parameters) == "df") {
            RET$parameter <- dist@parameters[[1]]
            names(RET$parameter) <- "df"
        }
        if (length(object@estimates) > 0)
            RET <- c(RET, object@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "ScalarIndependenceTest",
    definition = function(object) {
        stat <- object@statistic@teststatistic
        names(stat) <- "Z"

        dataname <- varnames(object@statistic)

        dist <- object@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(statistic = stat,
                    p.value = object@distribution@pvalue(stat),
                    alternative = object@statistic@alternative,
                    data.name = dataname,
                    method = paste(distname, object@method))
        if (length(object@nullvalue) > 0) {
            RET$null.value = object@nullvalue
            names(RET$null.value) <- object@parameter
        }
        if (length(object@estimates) > 0)
            RET <- c(RET, object@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "ScalarIndependenceTestConfint",
    definition = function(object) {
        stat <- object@statistic@teststatistic
        names(stat) <- "Z"
        dist <- object@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        dataname <- varnames(object@statistic)

        ci <- confint(object, level = object@conf.level)

        RET <- list(statistic = stat,
                    p.value = object@distribution@pvalue(stat),
                    alternative = object@statistic@alternative,
                    data.name = dataname,
                    method = paste(distname, object@method),
                    conf.int = ci$conf.int,
                    estimate = ci$estimate)
        if (length(object@nullvalue) > 0) {
            RET$null.value = object@nullvalue
            names(RET$null.value) <- object@parameter
        }
        if (length(object@estimates) > 0)
            RET <- c(RET, object@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

print.ci <- function(x, ...) {
    if(!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")),
            "percent confidence interval:\n",
            format(c(x$conf.int[1], x$conf.int[2])), "\n")
    }
    if(!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
}

print.MCp <- function(x, ...) {
    p <- x
    attributes(p) <- NULL
    print(p)
    ci <- list(conf.int = attr(x, "conf.int"))
    class(ci) <- "ci"
    print(ci)
}

print.cutpoint <- function(x, ...) {
    cat(paste0("  ", dQuote("best"), " cutpoint: ", x$label, "\n"))
    if (!is.null(x$covariable))
        cat(paste0("       covariable: ", x$covariable, "\n"))
}
