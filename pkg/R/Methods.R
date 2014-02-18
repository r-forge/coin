### generic method for asymptotic null distributions
setGeneric("AsymptNullDistribution", function(object, ...)
    standardGeneric("AsymptNullDistribution"))

### method for scalar test statistics
setMethod(f = "AsymptNullDistribution",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, ...) {

              new("AsymptNullDistribution",
                  p = function(q) pnorm(q),
                  q = function(p) qnorm(p),
                  d = function(x) dnorm(x),
                  pvalue = function(q)
                      switch(object@alternative,
                          "less"      = pnorm(q),
                          "greater"   = 1 - pnorm(q),
                          "two.sided" = 2 * min(pnorm(q), 1 - pnorm(q))),
                  support = function(p = 1e-5) c(qnorm(p), qnorm(1 - p)),
                  name = "Univariate Normal Distribution")
          }
)

### just a wrapper
pmv <- function(lower, upper, mean, corr, ...) {
    if (length(corr) > 1L)
        pmvnorm(lower = lower, upper = upper, mean = mean, corr = corr, ...)
    else
        pmvnorm(lower = lower, upper = upper, mean = mean, sigma = 1L, ...)
}

### method for max-type test statistics
setMethod(f = "AsymptNullDistribution",
          signature = "MaxTypeIndependenceTestStatistic",
          definition = function(object, ...) {

              corr <- cov2cor(covariance(object))
              pq <- length(expectation(object))
              pmvn <- function(q) {
                  p <- switch(object@alternative,
                           "less"      = pmv(lower = q, upper = Inf,
                                             mean = rep.int(0L, pq),
                                             corr = corr, ...),
                           "greater"   = pmv(lower = -Inf, upper = q,
                                             mean = rep.int(0L, pq),
                                             corr = corr, ...),
                           "two.sided" = pmv(lower = -abs(q), upper = abs(q),
                                             mean = rep.int(0L, pq),
                                             corr = corr, ...))
                  error <- attr(p, "error")
                  attr(p, "error") <- NULL
                  ci <- c(max(0, p - error), min(p + error, 1))
                  attr(ci, "conf.level") <- 0.99
                  attr(p, "conf.int") <- ci
                  class(p) <- "MCp"
                  p
              }
              qmvn <- function(p) {
                  if (length(corr) > 1L)
                      q <- qmvnorm(p, mean = rep.int(0L, pq), corr = corr,
                                   tail = "both.tails", ...)$quantile
                  else
                      q <- qmvnorm(p, mean = rep.int(0L, pq), sigma = 1L,
                                   tail = "both.tails", ...)$quantile
                  attributes(q) <- NULL
                  q
              }

              new("AsymptNullDistribution",
                  p = pmvn,
                  q = qmvn,
                  d = function(x) dmvnorm(x),
                  pvalue = function(q) {
                      p <- 1 - pmvn(q)
                      attr(p, "conf.int") <- 1 - attr(p, "conf.int")[2L:1L]
                      attr(attr(p, "conf.int"), "conf.level") <- 0.99
                      class(p) <- "MCp"
                      p
                  },
                  support = function(p = 1e-5) c(qmvn(p), qmvn(1 - p)),
                  name = "Multivariate Normal Distribution",
                  parameters = list(corr = corr))
          }
)

### method for quad-type test statistics
setMethod(f = "AsymptNullDistribution",
          signature = "QuadTypeIndependenceTestStatistic",
          definition = function(object, ...) {

              new("AsymptNullDistribution",
                  p = function(q) pchisq(q, df = object@df),
                  q = function(p) qchisq(p, df = object@df),
                  d = function(d) dchisq(d, df = object@df),
                  pvalue = function(q) 1 - pchisq(q, df = object@df),
                  support = function(p = 1e-5) c(0, qchisq(1 - p, df = object@df)),
                  name = "Chi-Squared Distribution",
                  parameters = list(df = object@df))
          }
)


### generic method for exact null distributions
setGeneric("ExactNullDistribution", function(object, ...)
    standardGeneric("ExactNullDistribution"))

### method for scalar test statistics
setMethod(f = "ExactNullDistribution",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, algorithm = c("auto", "shift", "split-up"),
                                ...) {

              algorithm <- match.arg(algorithm)

              if (object@paired) {
                  if (algorithm == "split-up")
                      stop("split-up algorithm not implemented for paired samples")
                  int <- is_integer(object@ytrans[, 1L], ...)
                  if (int)
                      SR_shift_1sample(object, fact = attr(int, "fact"))
                  else
                      stop("cannot compute exact distribution with real-valued scores")
              } else if (is_2sample(object)) {
                  if (algorithm == "split-up")
                      vdW_split_up_2sample(object)
                  else {
                      int <- is_integer(object@ytrans[, 1L], ...)
                      if (int)
                          SR_shift_2sample(object, fact = attr(int, "fact"))
                      else if (algorithm == "auto")
                          vdW_split_up_2sample(object)
                      else
                          stop("cannot compute exact distribution with real-valued scores")
                  }
              } else
                  stop(sQuote("object"), " is not a two sample problem")
          }
)

### method for quad-type test statistics
setMethod(f = "ExactNullDistribution",
          signature = "QuadTypeIndependenceTestStatistic",
          definition = function(object, algorithm = c("auto", "shift", "split-up"),
                                ...) {

              algorithm <- match.arg(algorithm)

              if (object@paired) {
                  if (algorithm == "split-up")
                      stop("split-up algorithm not implemented for paired samples")
                  int <- is_integer(object@ytrans[, 1L], ...)
                  if (int)
                      SR_shift_1sample(object, fact = attr(int, "fact"))
                  else
                      stop("cannot compute exact distribution with real-valued scores")
              } else if (is_2sample(object)) {
                  if (algorithm == "split-up")
                      stop("split-up algorithm not implemented for quad-type tests")
                  else {
                      int <- is_integer(object@ytrans[, 1L], ...)
                      if (int)
                          SR_shift_2sample(object, fact = attr(int, "fact"))
                      else if (algorithm == "auto")
                          stop("split-up algorithm not implemented for quad-type tests")
                      else
                          stop("cannot compute exact distribution with real-valued scores")
                  }
              } else
                  stop(sQuote("object"), " is not a two sample problem")
          }
)


### generic method for approximate null distributions
setGeneric("ApproxNullDistribution", function(object, ...)
    standardGeneric("ApproxNullDistribution"))

MCfun <- function(x, y, w, b, B) {
    ### expand observations for non-unit weights
    if (!is_unity(w)) {
        indx <- rep.int(seq_along(w), w)
        x <- x[indx, , drop = FALSE]
        y <- y[indx, , drop = FALSE]
        b <- b[indx]
    }
    .Call("R_MonteCarloIndependenceTest", x, y, as.integer(b), as.integer(B),
          PACKAGE = "coin")
}

### method for scalar test statistics
setMethod(f = "ApproxNullDistribution",
          signature = "ScalarIndependenceTestStatistic",
          definition = function(object, B = 10000, ...) {

              pls <- plsraw <-
                  MCfun(object@xtrans, object@ytrans, object@weights,
                        as.integer(object@block), as.integer(B))

              ### <FIXME> can transform p, q, x instead of those </FIXME>
              pls <- sort((pls - expectation(object)) / sqrt(variance(object)))

              new("ApproxNullDistribution",
                  p = function(q) {
                      p <- mean(LE(pls, q))
                      attr(p, "conf.int") <-
                          binom.test(round(p * B), B, conf.level = 0.99)$conf.int
                      class(p) <- "MCp"
                      p
                  },
                  q = function(p) {
                      quantile(pls, probs = p, names = FALSE, type = 1L)
                  },
                  d = function(x) {
                      tmp <- abs(pls - x)
                      mean(tmp == tmp[which.min(tmp)])
                  },
                  pvalue = function(q) {
                      p <- switch(object@alternative,
                               "less"      = mean(LE(pls, q)),
                               "greater"   = mean(GE(pls, q)),
                               "two.sided" = mean(GE(abs(pls), abs(q))))
                      attr(p, "conf.int") <-
                          binom.test(round(p * B), B, conf.level = 0.99)$conf.int
                      class(p) <- "MCp"
                      p
                  },
                  support = function(raw = FALSE) {
                      if (raw)
                          plsraw
                      else
                          sort(unique(drop(pls)))
                  },
                  name = "Monte Carlo Distribution")
          }
)

### method for max-type test statistics
setMethod(f = "ApproxNullDistribution",
          signature = "MaxTypeIndependenceTestStatistic",
          definition = function(object, B = 10000, ...) {

              pls <- plsraw <-
                  MCfun(object@xtrans, object@ytrans, object@weights,
                        as.integer(object@block), as.integer(B))

              fun <- switch(object@alternative,
                  "less" = min,
                  "greater" = max,
                  "two.sided" = function(x) max(abs(x))
              )

              dcov <- sqrt(variance(object))
              expect <- expectation(object)
              pls <- (pls - expect) / dcov

              ### <FIXME>
              ### pls is a rather large object (potentially)
              ### try not to copy it too often -- abs() kills you
              ### </FIXME>

              pmaxmin <- function() {
                  pls <- switch(object@alternative,
                      "less" = do.call("pmin.int", as.data.frame(t(pls))),
                      "greater" = do.call("pmax.int", as.data.frame(t(pls))),
                      "two.sided" = do.call("pmax.int", as.data.frame(t(abs(pls)))))
                  sort(pls)
              }

              new("ApproxNullDistribution",
                  p = function(q) {
                      p <- switch(object@alternative,
                               "less" = mean(colSums(GE(pls, q)) == nrow(pls)),
                               "greater" = mean(colSums(LE(pls, q)) == nrow(pls)),
                               "two.sided" = mean(colSums(LE(abs(pls), q)) == nrow(pls)))
                      attr(p, "conf.int") <-
                          binom.test(round(p * B), B, conf.level = 0.99)$conf.int
                      class(p) <- "MCp"
                      p
                  },
                  q = function(p) {
                      quantile(pmaxmin(), probs = p, names = FALSE, type = 1L)
                  },
                  d = function(x) {
                      tmp <- abs(pmaxmin() - x)
                      mean(tmp == tmp[which.min(tmp)])
                  },
                  pvalue = function(q) {
                      p <- switch(object@alternative,
                               "less" = mean(colSums(LE(pls, q)) > 0),
                               "greater" = mean(colSums(GE(pls, q)) > 0),
                               "two.sided" = mean(colSums(GE(abs(pls), q)) > 0))
                      attr(p, "conf.int") <-
                          binom.test(round(p * B), B, conf.level = 0.99)$conf.int
                      class(p) <- "MCp"
                      p
                  },
                  support = function(raw = FALSE) {
                      if (raw)
                          plsraw
                      else
                          sort(unique(drop(pmaxmin())))
                  },
                  name = "Monte Carlo Distribution")
          }
)

### method for quad-type test statistics
setMethod(f = "ApproxNullDistribution",
          signature = "QuadTypeIndependenceTestStatistic",
          definition = function(object, B = 10000, ...) {

              pls <- plsraw <-
                  MCfun(object@xtrans, object@ytrans, object@weights,
                        as.integer(object@block), as.integer(B))

              dcov <- object@covarianceplus
              expect <- expectation(object)
              a <- pls - expect
              pls <- sort(rowSums(crossprod(a, dcov) * t(a)))

              new("ApproxNullDistribution",
                  p = function(q) {
                      p <- mean(LE(pls, q))
                      attr(p, "conf.int") <-
                          binom.test(round(p * B), B, conf.level = 0.99)$conf.int
                      class(p) <- "MCp"
                      p
                  },
                  q = function(p) {
                      quantile(pls, probs = p, names = FALSE, type = 1L)
                  },
                  d = function(x) {
                      tmp <- abs(pls - x)
                      mean(tmp == tmp[which.min(tmp)])
                  },
                  pvalue = function(q) {
                      p <- mean(GE(pls, q))
                      attr(p, "conf.int") <-
                          binom.test(round(p * B), B, conf.level = 0.99)$conf.int
                      class(p) <- "MCp"
                      p
                  },
                  support = function(raw = FALSE) {
                      if (raw)
                          plsraw
                      else
                          sort(unique(drop(pls)))
                  },
                  name = "Monte Carlo Distribution")
          }
)


confint.ScalarIndependenceTestConfint <- function(object, parm, level = 0.95,
    ...) {
        if ("level" %in% names(match.call()))
            x <- object@confint(level)
        else
            x <- object@confint(object@conf.level)
        class(x) <- "ci"
        x
}
