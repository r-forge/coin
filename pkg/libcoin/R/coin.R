
.testcoin <- function(it) {
    ret <- list(statistic = c(statistic(it)), p.value = c(pvalue(it)),
                LinearStatistic = c(statistic(it, "linear")),
                Expectation = c(expectation(it)),
                Covariance = c(covariance(it)[!upper.tri(covariance(it))]))
    names(ret$Expectation) <- NULL
    names(ret$satistic) <- NULL
    ret
}

### perl -pe 's/([a-z]*_test)\(/lc(\"\1\",/g;' < coin-Ex.R > tmp.R
lc <- function(FUN, ...) {

    object <- do.call(FUN, list(...))
    blk <- object@statistic@block
    if (nlevels(blk) == 1) blk <- integer(0)
    B <- 10000L
    d <- object@distribution
    if (inherits(d, "AsymptNullDistribution")) B <- 0L
    w  <- as.integer(object@statistic@weights)
    s <- which(w > 0)
    lev <- LinStatExpCov(object@statistic@xtrans, 
                         object@statistic@ytrans, 
                         subset = s,
                         weights = w, 
                         block = blk, B = B)
    teststat <- "max"
    if (inherits(object@statistic, "QuadTypeIndependenceTestStatistic")) 
        teststat <- "quad"
        
    tst <- Test(lev, xtrafo = "id", type = teststat)
    tst$LinearStatistic <- lev$LinearStatistic
    tst$Expectation <- lev$Expectation
    tst$Covariance <- lev$Covariance  
    print(all.equal(tst, .testcoin(object)))
    object
}
