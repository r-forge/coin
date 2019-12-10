### methods for extracting p-values
setGeneric("pvalue",
    function(object, ...) {
        standardGeneric("pvalue")
    }
)

### <DEPRECATED>
### The "PValue" class was deprecated in 1.3-0 and at the same time this
### method was added as a temporary solution.  To be removed in 2.0-0.
setMethod("pvalue",
    signature = "PValue",
    definition = function(object, q, ...) {
        RET <- object@pvalue(q)
        class(RET) <- "pvalue"
        RET
    }
)
### </DEPRECATED>

setMethod("pvalue",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        RET <- object@pvalue(q)
        class(RET) <- c("pvalue", "numeric")
        RET
    }
)

setMethod("pvalue",
    signature = "ApproxNullDistribution",
    definition = function(object, q, ...) {
        RET <- callNextMethod(object, q, ...)
        attr(RET, "nresample") <- object@nresample
        RET
    }
)

setMethod("pvalue",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@distribution, object@statistic@teststatistic, ...)
    }
)

setMethod("pvalue",
    signature = "MaxTypeIndependenceTest",
    definition = function(object,
        method = c("global", "single-step", "step-down", "unadjusted"),
###     combinations = c("free", "restricted"), # placeholder
        distribution = c("joint", "marginal"),
        type = c("Bonferroni", "Sidak"), ...) {
            method <- match.arg(method,
                          choices = c("global", "single-step", "step-down",
                                      "unadjusted", "discrete"),
                          several.ok = TRUE)[1]
            if (method == "discrete")
                stop(sQuote(paste("method =", dQuote(method))),
                        " is defunct; see ", sQuote("?pvalue"))
            distribution <- match.arg(distribution)
            type <- match.arg(type)

            C <- attr(object@statistic@xtrans, "contrast")
            if (!is.null(C) && method != "global")
                warning("p-values may be incorrect due to violation\n",
                        "  of the subset pivotality condition")
            ## NOTE: Two ^^ spaces needed for correct rendering

            if (method == "global")
                callNextMethod(object, ...)
            else if (method == "single-step") {
                if (distribution == "joint")
                    joint(object, stepdown = FALSE, ...)
                else {
                    if (type == "Bonferroni")
                        marginal(object, stepdown = FALSE,
                                 bonferroni = TRUE, ...)
                    else
                        marginal(object, stepdown = FALSE,
                                 bonferroni = FALSE, ...)
                }
            } else if (method == "step-down") {
                if (distribution == "joint")
                    joint(object, stepdown = TRUE, ...)
                else {
                    if (type == "Bonferroni")
                        marginal(object, stepdown = TRUE,
                                 bonferroni = TRUE, ...)
                    else
                        marginal(object, stepdown = TRUE,
                                 bonferroni = FALSE, ...)
                }
            }
            else
                unadjusted(object, ...)
    }
)


### methods for extracting mid-p-values
setGeneric("midpvalue",
    function(object, ...) {
        standardGeneric("midpvalue")
    }
)

setMethod("midpvalue",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        RET <- object@midpvalue(q)
        class(RET) <- c("pvalue", "numeric")
        RET
    }
)

setMethod("midpvalue",
    signature = "ApproxNullDistribution",
    definition = function(object, q, ...) {
        RET <- callNextMethod(object, q, ...)
        attr(RET, "nresample") <- object@nresample
        RET
    }
)

setMethod("midpvalue",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@distribution, object@statistic@teststatistic, ...)
    }
)


### methods for extracting p-value intervals
setGeneric("pvalue_interval",
    function(object, ...) {
        standardGeneric("pvalue_interval")
    }
)

setMethod("pvalue_interval",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        object@pvalueinterval(q)
    }
)

setMethod("pvalue_interval",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@distribution, object@statistic@teststatistic, ...)
    }
)


### methods for extracting test size
setGeneric("size",
    function(object, ...) {
        standardGeneric("size")
    }
)

setMethod("size",
    signature = "NullDistribution",
    definition = function(object,
        alpha, type = c("p-value", "mid-p-value"), ...) {
            type <- match.arg(type)
            object@size(alpha, type)
    }
)

setMethod("size",
    signature = "IndependenceTest",
    definition = function(object,
        alpha, type = c("p-value", "mid-p-value"), ...) {
            callGeneric(object@distribution, alpha, type, ...)
    }
)


### methods for extracting the density function
setGeneric("dperm",
    function(object, x, ...) {
        standardGeneric("dperm")
    }
)

setMethod("dperm",
    signature = "NullDistribution",
    definition = function(object, x, ...) {
        object@d(x)
    }
)

setMethod("dperm",
    signature = "IndependenceTest",
    definition = function(object, x, ...) {
        callGeneric(object@distribution, x, ...)
    }
)


### methods for extracting the distribution function
setGeneric("pperm",
    function(object, q, ...) {
        standardGeneric("pperm")
    }
)

setMethod("pperm",
    signature = "NullDistribution",
    definition = function(object, q, ...) {
        object@p(q)
    }
)

setMethod("pperm",
    signature = "IndependenceTest",
    definition = function(object, q, ...) {
        callGeneric(object@distribution, q, ...)
    }
)


### methods for extracting the quantile function
setGeneric("qperm",
    function(object, p, ...) {
        standardGeneric("qperm")
    }
)

setMethod("qperm",
    signature = "NullDistribution",
    definition = function(object, p, ...) {
        object@q(p)
    }
)

setMethod("qperm",
    signature = "IndependenceTest",
    definition = function(object, p, ...) {
        callGeneric(object@distribution, p, ...)
    }
)


### methods for extracting random deviates
setGeneric("rperm",
    function(object, n, ...) {
        standardGeneric("rperm")
    }
)

setMethod("rperm",
    signature = "NullDistribution",
    definition = function(object, n, ...) {
        object@q(runif(n))
    }
)

setMethod("rperm",
    signature = "IndependenceTest",
    definition = function(object, n, ...) {
        callGeneric(object@distribution, n, ...)
    }
)


### methods for extracting the support
setGeneric("support",
    function(object, ...) {
        standardGeneric("support")
    }
)

setMethod("support",
    signature = "NullDistribution",
    definition = function(object, ...) {
        object@support(...)
    }
)

setMethod("support",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@distribution, ...)
    }
)


### methods for extracting test statistics etc.
setGeneric("statistic",
    function(object, ...) {
        standardGeneric("statistic")
    }
)

setMethod("statistic",
    signature = "IndependenceLinearStatistic",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"), ...) {
            type <- match.arg(type)
            nr <- ncol(object@xtrans)
            nc <- ncol(object@ytrans)
            dn <- statnames(object)$dimnames
            switch(type,
                "test" = stop(
                    sQuote(paste("type =", dQuote("test"))),
                    " not defined for objects of class ",
                    dQuote("IndependenceLinearStatistic")
                ),
                "linear" = matrix(
                    object@linearstatistic,
                    nrow = nr, ncol = nc, dimnames = dn
                ),
                "centered" = matrix(
                    object@linearstatistic - object@expectation,
                    nrow = nr, ncol = nc, dimnames = dn
                ),
                "standardized" = matrix(
                    (object@linearstatistic - object@expectation) /
                        sqrt(variance(object)),
                    nrow = nr, ncol = nc, dimnames = dn
                )
            )
    }
)

setMethod("statistic",
    signature = "IndependenceTestStatistic",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"), ...) {
            type <- match.arg(type)
            switch(type,
                "test" = {
                    object@teststatistic
                },
                "standardized" = {
                    nr <- ncol(object@xtrans)
                    nc <- ncol(object@ytrans)
                    dn <- statnames(object)$dimnames
                    matrix(
                        object@standardizedlinearstatistic,
                        nrow = nr, ncol = nc, dimnames = dn
                    )
                },
                callNextMethod(object, type, ...)
            )
    }
)

setMethod("statistic",
    signature = "IndependenceTest",
    definition = function(object,
        type = c("test", "linear", "centered", "standardized"), ...) {
            callGeneric(object@statistic, type, ...)
    }
)


### methods for extracting expectations
setGeneric("expectation",
    function(object, ...) {
        standardGeneric("expectation")
    }
)

setMethod("expectation",
    signature = "IndependenceLinearStatistic",
    definition = function(object, ...) {
        nr <- ncol(object@xtrans)
        nc <- ncol(object@ytrans)
        dn <- statnames(object)$dimnames
        matrix(
            object@expectation,
            nrow = nr, ncol = nc, dimnames = dn
        )
    }
)

setMethod("expectation",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@statistic, ...)
    }
)


### methods for extracting the covariance matrix
setGeneric("covariance",
    function(object, ...) {
        standardGeneric("covariance")
    }
)

### <DEPRECATED>
### Note: The "CovarianceMatrix", "Variance" and "VarCovar" classes were
### deprecated in 1.4-0.  To be removed in 2.0-0.
setMethod("covariance",
    signature = "CovarianceMatrix",
    definition = function(object, ...) {
        object@covariance
    }
)
### </DEPRECATED>

setMethod("covariance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, invert = FALSE, ...) {
        nm <- statnames(object)$names
        if (invert) {
            mp <- .Call(R_MPinv_sym, object@covariance, 0L, sqrt_eps)
            .Call(R_unpack_sym, mp$MPinv, nm, 0L)
        } else
            .Call(R_unpack_sym, object@covariance, nm, 0L)
    }
)

setMethod("covariance",
    signature = "QuadTypeIndependenceTestStatistic",
    definition = function(object, invert = FALSE, ...) {
        nm <- statnames(object)$names
        if (invert)
            .Call(R_unpack_sym, object@covarianceplus, nm, 0L)
        else
            callNextMethod(object, invert, ...)
    }
)

setMethod("covariance",
    signature = "IndependenceTest",
    definition = function(object, invert = FALSE, ...) {
        callGeneric(object@statistic, invert, ...)
    }
)


### methods for extracting the variances
setGeneric("variance",
    function(object, ...) {
        standardGeneric("variance")
    }
)

### <DEPRECATED>
### Note: The "CovarianceMatrix", "Variance" and "VarCovar" classes were
### deprecated in 1.4-0.  To be removed in 2.0-0.
setMethod("variance",
    signature = "Variance",
    definition = function(object, ...) {
        object@variance
    }
)

setMethod("variance",
    signature = "CovarianceMatrix",
    definition = function(object, ...) {
        diag(object@covariance)
    }
)
### </DEPRECATED>

setMethod("variance",
    signature = "IndependenceLinearStatistic",
    definition = function(object, ...) {
        nr <- ncol(object@xtrans)
        nc <- ncol(object@ytrans)
        dn <- statnames(object)$dimnames
        matrix(
            .Call(R_unpack_sym, object@covariance, NULL, 1L),
            nrow = nr, ncol = nc, dimnames = dn
        )
    }
)

setMethod("variance",
    signature = "IndependenceTest",
    definition = function(object, ...) {
        callGeneric(object@statistic, ...)
    }
)
