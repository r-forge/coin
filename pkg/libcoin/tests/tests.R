
library("libcoin")
library("coin")

.tr1d <- function(x) {
    tr <- switch(class(x)[1], 
                 "matrix" = x,
                 "integer" = as.double(x),
                 "numeric" = as.double(x),
                 "factor" = model.matrix(~ x - 1),
                 "ordered" = as.double((1:nlevels(x))[x]))
    storage.mode(tr) <- "double"
    list(tr = tr, whichNA = which(is.na(x)))
}

.tr2d <- function(x) {
    if (class(x)[1] == "matrix") {
        stop("not yet implemented")
    } else if (class(x)[1] %in% c("integer", "numeric")) {
        ux <- sort(unique(x))
        ix <- match(x, ux)
        tr <- matrix(c(0, ux), ncol = 1)
    } else if (class(x)[1] == "factor") {
        ux <- factor(levels(x))
        ix <- unclass(x)
        tr <- cbind(0, model.matrix(~ ux - 1))
    } else if (class(x)[1] == "ordered") {
        ux <- factor(levels(x), ordered = TRUE)
        ix <- unclass(x)
        tr <- matrix(c(0, 1:nlevels(x)), ncol = 1)
    }
    storage.mode(tr) <- "double"
    ix <- as.integer(ix)
    attr(ix, "levels") <- 1:max(ix)
    list(ix = ix, tr = tr, ux = ux, whichNA = which(is.na(x)))
}

.test1d <- function(x, y, weights = integer(0), subset = integer(0), block = integer(0),
                    xtrafo = c("id", "maxstat"), ...) {

    xtrafo <- match.arg(xtrafo)
    if (xtrafo == "id") {
        lev <- LinStatExpCov(.tr1d(x)$tr, .tr1d(y)$tr, weights = weights, subset = subset, block = block)
    } else {
        xtr <- .tr2d(x)
        lev <- LinStatExpCov(xtr$ix, .tr1d(y)$tr, weights = weights, subset = subset, block = block)
    }
    tst <- Test(lev, xtrafo = xtrafo, ordered = is.ordered(x) || is.numeric(x), ...)
    if (xtrafo == "maxstat")
        tst$index <- xtr$ux[tst$index]
    tst$LinearStatistic <- lev$LinearStatistic
    tst$Expectation <- lev$Expectation
    tst$Covariance <- lev$Covariance
    tst
}

.testcoin <- function(it) {
    ret <- list(statistic = c(statistic(it)), p.value = c(pvalue(it)),
                LinearStatistic = c(statistic(it, "linear")),
                Expectation = c(expectation(it)),
                Covariance = c(covariance(it)[!upper.tri(covariance(it))]))
    names(ret$Expectation) <- NULL
    ret
}

(tlc <- .test1d(iris$Sepal.Width, iris$Species, type = "quad"))
(tc <- .testcoin(independence_test(Species ~ Sepal.Width, data = iris, teststat = "quad")))
all.equal(tlc, tc)

set.seed(29)
(tlc <- .test1d(iris$Sepal.Width, iris$Species, type = "max"))
set.seed(29)
(tc <- .testcoin(independence_test(Species ~ Sepal.Width, data = iris, teststat = "max")))
all.equal(tlc, tc)


.test1d(iris$Sepal.Width, iris$Species, type = "max", xtrafo = "maxstat")
maxstat_test(Species ~ Sepal.Width, data = iris)

