
varnames <- function(object) {

    x <- object@x
    y <- object@y

    yordered <- sapply(y, is.ordered)
    ynames <- paste(colnames(y), ifelse(yordered, " (ordered)", ""), sep = "", 
                    collapse = ", ")

    if (length(x) == 1) {
        if (is.ordered(x[[1]])) {
            xnames <- paste("groups", paste(levels(x[[1]]), collapse = " < "))
        } else {
            if (is.factor(x[[1]])) {
                xnames <- paste("groups", 
                                paste(levels(x[[1]]), collapse = ", "))
            } else {
                xnames <- colnames(x)
            }
        }
    } else {
        xordered <- sapply(x, is.ordered)
        xnames <- paste(colnames(x), ifelse(xordered, "(ordered)", ""), 
                        sep = "", collapse = ", ")
    }

    if (nlevels(object@block) > 1)
        xnames <- paste(xnames,  
            paste("stratified by", attr(object@block, "blockname")))

    return(paste(ynames, "by", xnames, collapse = ""))
}

setMethod(f = "show", signature = "QuadTypeIndependenceTest", 
    definition = function(object) {

        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "T"
        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotical",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact") 

        dataname <- varnames(x@statistic)

        RET <- list(statistic = stat,
                    p.value = dist@pvalue(stat),
                    data.name = dataname,
                    method = paste(distname, x@method))
        if (length(dist@parameters) == 1 && 
                   names(dist@parameters) == "df") {
            RET$parameter <- dist@parameters[[1]]
            names(RET$parameter) <- "df"
        }
        if (length(x@statistic@estimates) > 0)
            RET <- c(RET, x@statistic@estimates)
        class(RET) <- "htest"
        print(RET)
    }
)

setMethod(f = "show", signature = "MaxTypeIndependenceTest",
    definition = function(object) {

        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "T"
        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotical",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        dataname <- varnames(x@statistic)

        RET <- list(statistic = stat,
                    p.value = x@distribution@pvalue(stat),
                    data.name = dataname,
                    method = paste(distname, x@method))
        if (length(x@statistic@estimates) > 0)
            RET <- c(RET, x@statistic@estimates)
        class(RET) <- "htest"
        print(RET)
    }
)

setMethod(f = "show", signature = "ScalarIndependenceTest",
    definition = function(object) {
          
        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "T"

        dataname <- varnames(x@statistic)

        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotical",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(statistic = stat,
                    p.value = x@distribution@pvalue(stat),
                    null.value = c(mu = x@nullvalue),
                    alternative = x@statistic@alternative,
                    data.name = dataname,
                    method = paste(distname, x@method))
        if (length(x@statistic@estimates) > 0)
            RET <- c(RET, x@statistic@estimates)
        class(RET) <- "htest"
        print(RET)
    }
)

setMethod(f = "show", signature = "ScalarIndependenceTestConfint",
    definition = function(object) {

        x <- object
        stat <- x@statistic@teststatistic
        names(stat) <- "T"
        dist <- x@distribution
        cld <- class(dist)
        attributes(cld) <- NULL
        distname <- switch(cld,
            "AsymptNullDistribution" = "Asymptotical",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        dataname <- varnames(x@statistic)

        ci <- confint(object, level = object@conf.level)

        RET <- list(statistic = stat,
                    p.value = x@distribution@pvalue(stat),
                    null.value = c(mu = x@nullvalue),
                    alternative = x@statistic@alternative,
                    method = paste(distname, x@method),
                    data.name = dataname,
                    conf.int = ci$conf.int,
                    estimate = ci$estimate)
        if (length(x@statistic@estimates) > 0)
            RET <- c(RET, x@statistic@estimates)
        class(RET) <- "htest"
        print(RET)
    }
)
