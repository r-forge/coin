
### a generic test procedure for classical (and not so classical) tests
independence_test <- function(x, ...) UseMethod("independence_test")

independence_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("independence_test", c(list(x = x), list(...)))
    return(RET)

}

independence_test.table <- function(x, distribution = c("asympt", "approx"), ...) {

    distribution <- match.arg(distribution)
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (distribution == "asympt") {
        df <- as.data.frame(x)
        if (ncol(df) == 3)
            x <- new("IndependenceProblem", x = df[1], y = df[2], block = NULL, 
                     weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            x <- new("IndependenceProblem", x = df[1], y = df[2], block = df[[3]], 
                     weights = df[["Freq"]])
        }
    } else {
        x <- table2df(x)
        if (ncol(x) == 3) {
            attr(x[[3]], "blockname") <- colnames(x)[3]
            x <- new("IndependenceProblem", x = x[1], y = x[2], block = x[[3]])
        } else {
            x <- new("IndependenceProblem", x = x[1], y = x[2], block = NULL) 
        }
    }
    ### </FIXME>
    RET <- do.call("independence_test", c(list(x = x, distribution = distribution), 
                   list(...)))
    return(RET)
}


independence_test.IndependenceProblem <- function(x,
    teststat = c("maxtype", "quadtype", "scalar"),
    distribution = c("asympt", "approx", "exact"), 
    alternative = c("two.sided", "less", "greater"), 
    xtrafo = trafo, ytrafo = trafo, check = NULL, ...) {

    teststat <- match.arg(teststat)
    alternative <- match.arg(alternative)
    distribution <- match.arg(distribution) 

    ### transform data if requested and setup a test problem
    itp <- new("IndependenceTestProblem", x, xtrafo = xtrafo, 
               ytrafo = ytrafo, ...)

    if (!is.null(check)) {
        if (is.function(check)) {
            if (!check(itp))
                stop(sQuote("check"), " failed")
        } else {
            stop(sQuote("check"), " is not a function")
        }
    }

    ### check type of test statistic and alternative
    scalar <- is_scalar(itp)

    if (!scalar) {
        if (teststat == "scalar") {
            warning("Length linear statistic > 1, using ",
                    sQuote("maxtype"), " test statistic")
            teststat <- "maxtype"
        }
        if (alternative != "two.sided")
            warning(sQuote("alternative"), " is ignored for ", 
                    teststat, " type test statistics")
    } else {
        if (teststat == "maxtype") teststat <- "scalar"
    }

    ### compute linear statistic, conditional expectation and
    ### conditional covariance
    its <- new("IndependenceTestStatistic", itp)

    ### compute test statistic and corresponding null distribution
    RET <- switch(teststat,
        "scalar" = {
            ts <- new("ScalarIndependenceTestStatistic", its, 
                      alternative = alternative)

            nd <- switch(distribution,
                "asympt" = AsymptNullDistribution(ts, ...),
                "exact"  = ExactNullDistribution(ts, ...),
                "approx" = ApproxNullDistribution(ts, ...)
            )
            new("ScalarIndependenceTest", statistic = ts, distribution = nd)
        },
        "maxtype" = {
            ts <- new("MaxTypeIndependenceTestStatistic", its)
            nd <- switch(distribution,
                "asympt" = AsymptNullDistribution(ts, ...),
                "exact"  = ExactNullDistribution(ts, ...),
                "approx" = ApproxNullDistribution(ts, ...)
            )
            new("MaxTypeIndependenceTest", statistic = ts, distribution = nd)
        },
        "quadtype" = {
            ts <- new("QuadTypeIndependenceTestStatistic", its)
            nd <- switch(distribution,
                "asympt" = AsymptNullDistribution(ts, ...),
                "exact"  = ExactNullDistribution(ts, ...),
                "approx" = ApproxNullDistribution(ts, ...)
            )
            new("QuadTypeIndependenceTest", statistic = ts, 
                distribution = nd)
        })

    ### return object inheriting from class `IndependenceTest'
    return(RET)
}


### OK, OK, here is the most prominent one ...
wilcox_test <- function(x, ...) UseMethod("wilcox_test")

wilcox_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("wilcox_test", c(list(x = x), list(...)))
    return(RET)

}

wilcox_test.IndependenceProblem <- function(x,  
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), 
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(x) {
        if (!(is_2sample(x) && is_numeric_y(x)))
            stop(sQuote("x"), " does not represent a two sample problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- match.arg(distribution) 

    RET <- independence_test(x, teststat = "scalar", 
        alternative = alternative, distribution = distribution, 
        ytrafo = function(data) trafo(data, numeric_trafo = rank), 
        check = check, ...)

    RET@nullvalue <- 0
    RET@method <- "Wilcoxon Mann-Whitney Rank Sum Test"

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_location(RET@statistic, RET@distribution, 
                             level = level, approx = (distribution == "asympt"))
        RET@conf.level <- conf.level
    }
    return(RET)
}


### normal quantiles (van der Waerden) test
normal_test <- function(x, ...) UseMethod("normal_test")

normal_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("normal_test", c(list(x = x), list(...)))
    return(RET)
}   

normal_test.IndependenceProblem <- function(x,  
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), 
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(x) {
        if (!(is_2sample(x) && is_numeric_y(x)))
            stop(sQuote("x"), " does not represent a two sample problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- match.arg(distribution) 

    RET <- independence_test(x, teststat = "scalar", 
        alternative = alternative, distribution = distribution, 
        ytrafo = function(data) trafo(data, numeric_trafo = normal_trafo), 
        check = check, ...)

    RET@nullvalue <- 0
    RET@method <- "Normal Quantile (van der Waerden) Test"

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_location(RET@statistic, RET@distribution,
                             level = level, approx = (distribution == "asympt"))
        RET@conf.level <- conf.level
    }
    return(RET)
}


### median test
median_test <- function(x, ...) UseMethod("median_test")

median_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("median_test", c(list(x = x), list(...)))
    return(RET)
}   

median_test.IndependenceProblem <- function(x,     
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), 
    conf.int = FALSE, conf.level = 0.95, ...) {

    check <- function(x) {
        if (!(is_2sample(x) && is_numeric_y(x)))
            stop(sQuote("x"), " does not represent a two sample problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- match.arg(distribution) 

    RET <- independence_test(x, teststat = "scalar", 
        alternative = alternative, distribution = distribution, 
        ytrafo = function(data) trafo(data, numeric_trafo = median_trafo), 
        check = check, ...)
 
    RET@nullvalue <- 0
    RET@method <- "Median Test"

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_location(RET@statistic, RET@distribution,
                             level = level, approx = (distribution == "asympt"))
        RET@conf.level <- conf.level
    }
    return(RET)
}


### Ansari-Bradley test
ansari_test <- function(x, ...) UseMethod("ansari_test")

ansari_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)  
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("ansari_test", c(list(x = x), list(...)))
    return(RET)
}   

ansari_test.IndependenceProblem <- function(x,
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),  
    conf.int = FALSE, conf.level = 0.95, ...) {     

    check <- function(x) {
        if (!(is_2sample(x) && is_numeric_y(x)))
            stop(sQuote("x"), " does not represent a two sample problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- match.arg(distribution)

    RET <- independence_test(x, teststat = "scalar",
        alternative = alternative, distribution = distribution,
        ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
        check = check, ...)
 
    RET@nullvalue <- 1
    RET@method <- "Ansari-Bradley Test"

    if (conf.int) {
        RET <- new("ScalarIndependenceTestConfint", RET)
        RET@confint <- function(level)
            confint_scale(RET@statistic, RET@distribution,
                          level = level, approx = (distribution == "asympt"))
        RET@conf.level <- conf.level
    }
    return(RET)
}


### Logrank test
logrank_test <- function(x, ...) UseMethod("logrank_test")

logrank_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)  
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("logrank_test", c(list(x = x), list(...)))
    return(RET)
}
    
logrank_test.IndependenceProblem <- function(x,  
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), ...) {

    alternative <- match.arg(alternative)
    distribution <- match.arg(distribution)

    check <- function(x) {
        if (!(is_Ksample(x) && is_censored_y(x)))
            stop(sQuote("x"), 
                 " does not represent a K sample problem with censored data")
        return(TRUE)
    }

    scalar <- FALSE
    if (is.factor(x@x[[1]])) scalar <- nlevels(x@x[[1]]) == 2

    RET <- independence_test(x, 
        teststat = ifelse(scalar, "scalar", "quadtype"), 
        distribution = distribution, check = check, ...)
 
    if (extends(class(RET@statistic), "ScalarIndependenceTest"))
        RET@nullvalue <- 0

    RET@method <- "Logrank Test"
    return(RET)
}


### Kruskal-Wallis test
kruskal_test <- function(x, ...) UseMethod("kruskal_test")

kruskal_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("kruskal_test", c(list(x = x), list(...)))
    return(RET)
}   

kruskal_test.IndependenceProblem <- function(x,  
    distribution = c("asympt", "approx"), ...) {

    check <- function(x) {
        if (!(is_Ksample(x) && is_numeric_y(x)))
            stop(sQuote("x"), " does not represent a K sample problem")
        return(TRUE)
    }
 
    distribution <- match.arg(distribution)

    RET <- independence_test(x, 
        distribution = distribution, teststat = "quadtype",
        ytrafo = function(data) trafo(data, numeric_trafo = rank), 
        check = check, ...)

    RET@method <- paste("Kruskal-Wallis Test")
    return(RET)
}


### Fligner test
fligner_test <- function(x, ...) UseMethod("fligner_test")

fligner_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("fligner_test", c(list(x = x), list(...)))
    return(RET)
}   

fligner_test.IndependenceProblem <- function(x,  
    distribution = c("asympt", "approx"), ...) {

    check <- function(x) {
        if (!(is_Ksample(x) && is_numeric_y(x)))
            stop(sQuote("x"), " does not represent a K sample problem")
        return(TRUE)
    }
 
    distribution <- match.arg(distribution)

    ### eliminate location differences (see `stats/R/fligner.test')
    x@y[[1]] <- x@y[[1]] - tapply(x@y[[1]], x@x[[1]], median)[x@x[[1]]]

    RET <- independence_test(x,  
        distribution = distribution, teststat = "quadtype",
        ytrafo = function(data) trafo(data, numeric_trafo = fligner_trafo), 
        check = check, ...)

    RET@method <- paste("Fligner-Killeen Test")
    return(RET)
}


### Spearman test
spearman_test <- function(x, ...) UseMethod("spearman_test")

spearman_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("spearman_test", c(list(x = x), list(...)))
    return(RET)
}   

spearman_test.IndependenceProblem <- function(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx"), ...) {

    check <- function(x) {
        if (!is_corr(x))
            stop(sQuote("x"), " does not represent a univariate correlation problem")
        return(TRUE)
    }

    alternative <- match.arg(alternative)
    distribution <- match.arg(distribution) 

    RET <- independence_test(x, 
        teststat = "scalar", alternative = alternative, 
        distribution = distribution, 
        xtrafo = function(data) trafo(data, numeric_trafo = rank),
        ytrafo = function(data) trafo(data, numeric_trafo = rank), 
        check = check, ...)

    RET@nullvalue <- 0
    RET@method <- paste("Spearman Correlation Test")
    return(RET)
}


### Generalised Cochran-Mantel-Haenzel Test
cmh_test <- function(x, ...) UseMethod("cmh_test")

cmh_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("cmh_test", c(list(x = x), list(...)))
    return(RET)
}   

cmh_test.table <- function(x, distribution = c("asympt", "approx"), ...) {

    distribution <- match.arg(distribution)
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (distribution == "asympt") {
        df <- as.data.frame(x)
        if (ncol(df) == 3)
            x <- new("IndependenceProblem", x = df[1], y = df[2], block = NULL, 
                     weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            x <- new("IndependenceProblem", x = df[1], y = df[2], block = df[[3]], 
                     weights = df[["Freq"]])
        }
    } else {
        x <- table2df(x)
        if (ncol(x) == 3) {
            attr(x[[3]], "blockname") <- colnames(x)[3]
            x <- new("IndependenceProblem", x = x[1], y = x[2], block = x[[3]])
        } else {
            x <- new("IndependenceProblem", x = x[1], y = x[2], block = NULL) 
        }
    }
    ### </FIXME>
    RET <- do.call("cmh_test", c(list(x = x, distribution = distribution), list(...)))
    return(RET)
}

cmh_test.IndependenceProblem <- function(x, 
    distribution = c("asympt", "approx"), ...) {

    check <- function(x) {
        if (!is_contingency(x))
            stop(sQuote("x"), " does not represent a contingency problem")
        return(TRUE)
    }
    n <- nrow(x@x)

    distribution <- match.arg(distribution)

    RET <- independence_test(x, 
        teststat = "quadtype", distribution = distribution, check = check, 
        ...)

    RET@method <- paste("Generalised Cochran-Mantel-Haenszel Test")
    return(RET)
}


### Pearsons Chi-Squared Test
chisq_test <- function(x, ...) UseMethod("chisq_test")

chisq_test.formula <- function(formula, data = list(), subset = NULL, 
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("chisq_test", c(list(x = x), list(...)))
    return(RET)
}   

chisq_test.table <- function(x, distribution = c("asympt", "approx"), ...) {

    distribution <- match.arg(distribution)
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (distribution == "asympt") {
        df <- as.data.frame(x)
        if (ncol(df) == 3)
            x <- new("IndependenceProblem", x = df[1], y = df[2], block = NULL, 
                     weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            x <- new("IndependenceProblem", x = df[1], y = df[2], block = df[[3]], 
                     weights = df[["Freq"]])
        }
    } else {
        x <- table2df(x)
        if (ncol(x) == 3) {
            attr(x[[3]], "blockname") <- colnames(x)[3]
            x <- new("IndependenceProblem", x = x[1], y = x[2], block = x[[3]])
        } else {
            x <- new("IndependenceProblem", x = x[1], y = x[2], block = NULL) 
        }
    }
    ### </FIXME>
    RET <- do.call("chisq_test", c(list(x = x, distribution = distribution), list(...)))
    return(RET)
}

chisq_test.IndependenceProblem <- function(x,  
    distribution = c("asympt", "approx"), ...) {

    check <- function(x) {
        if (!is_contingency(x))
            stop(sQuote("x"), " does not represent a contingency problem")
        if (nlevels(x@block) != 1)
            stop(sQuote("x"), " contains blocks: use ", 
                 sQuote("cmh_test"), " instead")
        return(TRUE)
    }
    n <- sum(x@weights)

    distribution <- match.arg(distribution)

    RET <- independence_test(x, 
        teststat = "quadtype", distribution = "asympt", check = check, ...)

    ### use the classical chisq statistic based on Pearson 
    ### residuals (O - E)^2 / E
    ### see Th. 3.1 and its proof in Strasser & Weber (1999).

    RET@statistic@teststatistic <- 
        RET@statistic@teststatistic * n / (n - 1)
    RET@statistic@covariance <- 
        RET@statistic@covariance * (n - 1) / n
    RET@statistic@covarianceplus <- MPinv(RET@statistic@covariance)$MPinv

    if (distribution == "approx") {
        nd <- ApproxNullDistribution(RET@statistic, ...)
        RET <- new("QuadTypeIndependenceTest", statistic = RET@statistic,
                distribution = nd)
    }

    RET@method <- paste("Pearson's Chi-Squared Test")
    return(RET)
}


### Linear-by-Linear Association Test
lbl_test <- function(x, ...) UseMethod("lbl_test")

lbl_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("lbl_test", c(list(x = x), list(...)))
    return(RET)
}   

lbl_test.table <- function(x, distribution = c("asympt", "approx"), ...) {

    distribution <- match.arg(distribution)
    ### <FIXME> approx must be able to deal with weights </FIXME>
    if (distribution == "asympt") {
        df <- as.data.frame(x)
        if (nlevels(df[[1]]) > 2) df[[1]] <- ordered(df[[1]])
        if (nlevels(df[[2]]) > 2) df[[2]] <- ordered(df[[2]])
        if (ncol(df) == 3)
            x <- new("IndependenceProblem", x = df[1], y = df[2], block = NULL, 
                     weights = df[["Freq"]])
        if (ncol(df) == 4) {
            attr(df[[3]], "blockname") <- colnames(df)[3]
            x <- new("IndependenceProblem", x = df[1], y = df[2], block = df[[3]], 
                     weights = df[["Freq"]])
        }
    } else {
        x <- table2df(x)
        if (nlevels(x[[1]]) > 2) x[[1]] <- ordered(x[[1]])
        if (nlevels(x[[2]]) > 2) x[[2]] <- ordered(x[[2]])
        if (ncol(x) == 3) {
            attr(x[[3]], "blockname") <- colnames(x)[3]
            x <- new("IndependenceProblem", x = x[1], y = x[2], block = x[[3]])
        } else {
            x <- new("IndependenceProblem", x = x[1], y = x[2], block = NULL) 
        }
    }
    ### </FIXME>
    RET <- do.call("lbl_test", c(list(x = x, distribution = distribution), list(...)))
    return(RET)
}


lbl_test.IndependenceProblem <- function(x, xscores = NULL, 
    yscores = NULL, distribution = c("asympt", "approx"), ...) {

    check <- function(x) {
        if (!is_ordered(x))
            stop(sQuote("x"), 
                 " does not represent a problem with ordered data")
        return(TRUE)
    }

    if (!is.null(xscores)) {
        if (length(xscores) != nlevels(x@x[[1]]))
            stop(sQuote("xscores"), " don't match")
        attr(x@x[[1]], "scores") <- xscores
    }

    if (!is.null(yscores)) {
        if (length(yscores) != nlevels(x@y[[1]]))
            stop(sQuote("yscores"), " don't match")
        attr(x@y[[1]], "scores") <- yscores
    }

    distribution <- match.arg(distribution)

    RET <- independence_test(x, 
        teststat = "quadtype", distribution = distribution, check = check, ...)

    RET@method <- paste("Linear-by-Linear Association Test")
    return(RET)
}


### permutation test without transformations
perm_test <- function(x, ...) UseMethod("perm_test")

perm_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("perm_test", c(list(x = x), list(...)))  
    return(RET)
}

perm_test.IndependenceProblem <- function(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), ...) {

    alternative <- match.arg(alternative)
    distribution <- match.arg(distribution) 

    check <- function(x) {
        if (!(is_Ksample(x) && is_numeric_y(x)))
            stop(sQuote("x"), " does not represent a K sample problem")
        return(TRUE)
    }

    RET <- independence_test(x, 
        alternative = alternative, distribution = distribution, 
        check = check, ...)

    if (is_scalar(RET@statistic))
        RET@nullvalue <- 0
    RET@method <- paste(ifelse(is_scalar(RET@statistic), "2-", "K-"), 
                        paste("Sample Permutation Test"), sep = "")
    return(RET)
}

perm_test.SymmetryProblem <- function(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), ...) {

    class(x) <- "IndependenceProblem"
    RET <- perm_test(x, alternative = alternative, 
                     distribution = distribution, ...)
   return(RET)
}

### Contrast test
contrast_test <- function(x, ...) UseMethod("contrast_test")

contrast_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("contrast_test", c(list(x = x), list(...)))
    return(RET)
}

contrast_test.IndependenceProblem <- function(x, 
    cmatrix, distribution = c("asympt", "approx"), ...) {

    if (!(ncol(x@x) == 1 && is.factor(x@x[[1]])))
        stop(sQuote("x@x"), " is not univariate or a factor")

    if  (!is.matrix(cmatrix) || nrow(cmatrix) != nlevels(x@x[[1]]))
        stop(sQuote("cmatrix"), " is not a matrix with ", nlevels(x), " rows")

    if (is.null(colnames(cmatrix)))
        colnames(cmatrix) <- paste("C", 1:ncol(cmatrix), sep = "")

    distribution <- match.arg(distribution)

    xtrafo <- function(data) trafo(data) %*% cmatrix

    RET <- independence_test(x, teststat = "maxtype",
        distribution = distribution, xtrafo = xtrafo, ...)
    RET@method <- paste("Contrast Test")
    
    return(RET)
}


### Maxstat test
maxstat_test <- function(x, ...) UseMethod("maxstat_test")

maxstat_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("maxstat_test", c(list(x = x), list(...)))
    return(RET)
}

maxstat_test.IndependenceProblem <- function(x, 
    distribution = c("asympt", "approx"), ...) {

    check <- function(x) {
        if (!is_ordered_x(x))
            stop("all input variables need to be of class ", sQuote("ordered"),
                 " or ", sQuote("numeric"))
        return(TRUE)
    }

    distribution <- match.arg(distribution)

    xtrafo <- function(data) trafo(data, numeric_trafo = maxstat_trafo)

    RET <- independence_test(x, teststat = "maxtype",
        distribution = distribution, xtrafo = xtrafo, check = check, ...)

    wm <- which.max(apply(abs(statistic(RET, "standardized")), 1, max))
    whichvar <- attr(RET@statistic@xtrans, "assign")[wm]
    maxcontr <- RET@statistic@xtrans[,wm]
    estimate <- max(x@x[[whichvar]][maxcontr > 0])
    names(estimate) <- colnames(x@x)[whichvar]
    RET@statistic@estimates <- list(estimate = estimate)
    RET@method <- paste("Maxstat Test")
    
    return(RET)
}


### EXPERIMENTAL ###

### a generic test procedure for classical (and not so classical) tests
symmetry_test <- function(x, ...) UseMethod("symmetry_test")

symmetry_test.formula <- function(formula, data = list(), subset = NULL,
    ...) {

    d <- formula2data(formula, data, subset, ...)
    x <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("symmetry_test", c(list(x = x), list(...)))
    return(RET)

}

symmetry_test.SymmetryProblem <- function(x,
    teststat = c("maxtype", "quadtype", "scalar"),
    distribution = c("asympt", "approx", "exact"), 
    alternative = c("two.sided", "less", "greater"), 
    xtrafo = trafo, ytrafo = trafo, check = NULL, ...) {
    class(x) <- "IndependenceProblem"
    independence_test(x, teststat, distribution, alternative, xtrafo,
                      ytrafo, check, ...)
}

symmetry_test.table <- function(x, ...) {
    x <- table2df_sym(x)
    x <- new("SymmetryProblem", x = x["groups"], y = x["response"])
    RET <- do.call("symmetry_test", c(list(x = x), list(...))) 
    return(RET)
}

### Friedman-Test
friedman_test <- function(x, ...) UseMethod("friedman_test")

friedman_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    x <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("friedman_test", c(list(x = x), list(...)))
    return(RET)
}   

friedman_test.SymmetryProblem <- function(x, 
    distribution = c("asympt", "approx"), ...) {
    
    if (!is_completeblock(x))
        stop("Not an unreplicated complete block design")

    for (lev in levels(x@block))
        x@y[[1]][x@block == lev] <- rank(x@y[[1]][x@block == lev])

    distribution <- match.arg(distribution)

    RET <- symmetry_test(x, 
        distribution = distribution, teststat = "quadtype", ...)

    RET@method <- paste("Friedman Test")
    return(RET)
}

### Bowker-Test
bowker_test <- function(x, ...) UseMethod("bowker_test")

bowker_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    x <- new("SymmetryProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("bowker_test", c(list(x = x), list(...)))
    return(RET)
}   

bowker_test.table <- function(x, yscores = NULL, ...) {
    x <- table2df_sym(x)
    if (!is.null(yscores)) {
        x$response <- ordered(x$response)
        attr(x$response, "scores") <- yscores
    }
    x <- new("SymmetryProblem", x = x["groups"], y = x["response"])
    RET <- do.call("bowker_test", c(list(x = x), list(...))) 
    return(RET)
}

bowker_test.SymmetryProblem <- function(x, 
    distribution = c("asympt", "approx"), ...) {
    
    if (!is_completeblock(x))
        stop("Not an unreplicated complete block design")
    if (ncol(x@y) != 1 || !is.factor(x@y[[1]]))
        stop("Response variable is not a factor")

    distribution <- match.arg(distribution)

    RET <- symmetry_test(x, 
        distribution = distribution, teststat = "quadtype", ...)

    RET@method <- paste("Bowker Test")
    return(RET)
}


### Wilcoxon-Signed-Rank Test
wilcoxsign_test <- function(x, ...) UseMethod("wilcoxsign_test")

wilcoxsign_test.formula <- function(formula, data = list(), 
                                    subset = NULL, ...)
{
    d <- formula2data(formula, data, subset, ...)
    x <- new("IndependenceProblem", x = d$x, y = d$y, block = d$bl)
    RET <- do.call("wilcoxsign_test", c(list(x = x), list(...)))
    return(RET)
}   

wilcoxsign_test.IndependenceProblem <- function(x, 
    distribution = c("asympt", "approx"), ...) {

    y <- x@y[[1]]
    x <- x@x[[1]]
    if (!is.numeric(x))
        stop(sQuote("x"), " is not a numeric variable")
    if (!is.numeric(y))
        stop(sQuote("y"), " is not a numeric variable")

    block <- gl(length(x), 2)
    diffs <- x - y
    pos <- rank(abs(diffs)) * (diffs > 0)[x != y]
    neg <- rank(abs(diffs)) * (diffs < 0)[x != y]
    yy <- drop(as.vector(t(cbind(pos, neg))))
    xx <- factor(rep(c("pos", "neg"), length(x)))

    distribution <- match.arg(distribution)

    ip <- new("IndependenceProblem", x = data.frame(x = xx), 
              y = data.frame(y = yy), block = block)

    RET <- independence_test(ip,
        distribution = distribution, teststat = "scalar", ...)

    RET@method <- paste("Wilcoxon-Signed-Rank Test")
    RET@nullvalue <- 0
    return(RET)
}
