
.testcoin <- function(it) {
    ret <- list(TestStatistic = c(statistic(it)), p.value = c(pvalue(it)),
                LinearStatistic = c(statistic(it, "linear")),
                Expectation = c(expectation(it)),
                Covariance = c(covariance(it)[!upper.tri(covariance(it))]))
    names(ret$Expectation) <- NULL
    names(ret$TestStatistic) <- NULL
    names(ret$p.value) <- NULL
    ret
}

### perl -pe 's/([a-z]*_test)\(/lc(\"\1\",/g;' < coin-Ex.R > tmp.R
lc <- function(FUN, ...) {

    set.seed(29)
    object <- do.call(FUN, list(...))
    blk <- object@statistic@block
    if (nlevels(blk) == 1) blk <- integer(0)
    B <- 1000L
    d <- object@distribution
    if (inherits(d, "AsymptNullDistribution")) B <- 0L
    w  <- as.integer(object@statistic@weights)
    s <- integer(0)
    if (any(w <= 0))
        s <- which(w > 0)
    set.seed(29)
    lev <- LinStatExpCov(X = object@statistic@xtrans, 
                         Y = object@statistic@ytrans, 
                         subset = s,
                         weights = w, 
                         block = blk, B = B)
    ix <- iy <- 1:nrow(object@statistic@xtrans)
    if (max(ix) < 500) {
       attr(ix, "levels") <- attr(iy, "levels") <- 1:max(ix)
        lev2d <- LinStatExpCov(X = rbind(0, object@statistic@xtrans),
                         Y = rbind(0, object@statistic@ytrans),
                         ix = ix, iy = iy,
                         subset = s,
                         weights = w,
                         block = blk, B = ifelse(max(ix) < 50, B, 0L))
        slt <- c("LinearStatistic", "Expectation", "Covariance")
        stopifnot(all.equal(lev[slt], lev2d[slt]))
    }

    n <- sum(w)
    if (FUN == "chisq_test") lev$Covariance <- lev$Covariance * (n - 1) / n
    teststat <- "maximum"
    alternative <- "two.sided"
    if (inherits(object@statistic, "QuadTypeIndependenceTestStatistic")) {
        teststat <- "quadratic"
    } else {
        if (inherits(object@statistic, "ScalarIndependenceTestStatistic"))
            teststat <- "scalar"
        alternative <- object@statistic@alternative
    }
        
    tst <- doTest(lev, teststat = teststat, alternative = alternative)
    tst$LinearStatistic <- lev$LinearStatistic
    tst$Expectation <- lev$Expectation
    tst$Covariance <- lev$Covariance  
    print(all.equal(tst, .testcoin(object), scale = 1))
    object
}
