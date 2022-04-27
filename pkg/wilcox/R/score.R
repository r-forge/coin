
.wilcox_score <- function(x, y, level = 0.95, exact, ...) {

  z <- c(x, y)
  r <- rank(z)
  n.x <- length(x)
  n.y <- length(y)
  U <- sum(r[seq_along(x)]) - n.x * (n.x + 1) / 2
  if (U < 2 || U > (n.x * n.y - 2))
      stop("samples are (almost) non-overlapping, cannot compute confidence interval")

  ### PI = P(X < Y) -> logOR
  PI <- U / (n.x * n.y)
  PI2logOR <- function(x, PI)
      x + (exp(-x) * (PI + exp(2 * x) * (PI - 1) + exp(x)* (1 - 2 * PI)))
  ilogOR <- uniroot(PI2logOR, PI = PI, interval = 50 * c(-1, 1))$root

  N <- n.x + n.y
  alpha <- (1 - level) / 2

  ### use at most 500 levels in MASS::polr
  uz <- unique(z)
  ct <- sort(uz)
  zc <- cut(z, breaks = c(-Inf, ct, Inf), ordered = TRUE)[, drop = TRUE]
  ix <- rep(0:1, c(n.x, n.y))

  q <- -2000:2000 / 100
  ECDF <- ecdf(z)(sort(z)[-N])

  .zeta <- function(beta) {
      ### match marginal ECDF and conditional distributions,
      ### see doi:10.1002/sim.7890, page 3998 
      marg <- (n.x * plogis(q) + n.y * plogis(q - beta)) / (n.x + n.y)
      ### invert numerically
      s <- spline(x = q, y = marg, method = "hyman")
      return(approx(x = s$y, y = s$x, xout = ECDF, rule = 1)$y)
  }

  if (length(uz) > 500) {

      stopifnot(!exact)

      ### return Wald interval; optimize 
      ### hard-coded profile likelihood here

      ### nonparametric profile likelihood for beta
      prfll <- function(beta) {
          zeta <- .zeta(beta)
          off <- beta * ix
          ### likelihood contributions
          prb <- plogis(c(zeta, Inf)[zc] - off) - plogis(c(-Inf, zeta)[zc] - off)
          -sum(log(pmax(sqrt(.Machine$double.eps), prb)))
      }

      ### this is a one-dimensional problem
      if (length(uz) < 1000) {
          opt <- optimize(prfll, interval = c(ilogOR - 5, ilogOR + 5))
          logOR <- opt$minimum
      } else {
          logOR <- ilogOR
      }
      ### compute Hessian numerically
      idx <- seq(from = logOR - 5, to = logOR + 5, length.out = 1000)
      eval_prfll <- sapply(idx, prfll)
      sd1 <- sqrt(1 / splinefun(idx, eval_prfll)(logOR, deriv = 2))
      ### Wald interval
      Sci <- logOR + sd1 * qnorm(c(alpha, 1 - alpha))
      return(c(Sci, logOR))
  }

  ### good starting values for polr: approximate zeta for initial ilogOR
  izeta <- .zeta(ilogOR)
  parm0 <- c(izeta, ilogOR)

  ### fit proportional odds model
  d <- expand.grid(zc = sort(unique(zc)), grp = factor(0:1))
  grp <- factor(ix)
  d$w <- c(table(zc, grp))

  plr <- polr(zc ~ grp, data = d, weights = d$w, 
              start = parm0[c(length(parm0), 1:(length(parm0) - 1))], 
              Hess = FALSE)
  logOR <- coef(plr)
  Xgrp <- model.matrix(plr)[, names(logOR)]

  ### set-up design matrices
  X <- Matrix(0L, nrow = nrow(d), ncol = nlevels(zc))
  X[cbind(1:nrow(X), unclass(d$zc))] <- 1L
  X <- cbind(X, Xgrp)
  X[,ncol(X)] <- -X[, ncol(X)] 
  offl <- offr <- numeric(nrow(X))
  yi <- unclass(d$zc)
  offr[yi == max(yi)] <- 32
  offl[yi == min(yi)] <- -32

  mmr <- X[, -(ncol(X) - 1)]
  mml <- X[, -1]

  ## score function wrt logOR parameter in proportional odds model 
  ## for some parameters parm
  .sc <- function(parm, w) {
    lpl <- as(mml %*% parm + offl, "vector")
    lpr <- as(mmr %*% parm + offr, "vector")
    mmr[yi == max(yi),] <- 0
    mml[yi == min(yi),] <- 0

    Fr <- plogis(lpr)
    Fl <- plogis(lpl)
    fr <- dlogis(lpr)
    fl <- dlogis(lpl)
    sum(w * (fr - fl) / (Fr - Fl) * Xgrp)
  }

  ## standard deviation for logOR, for some parameters parm
  .sd <- function(parm, w) {
    lpl <- as(mml %*% parm + offl, "vector")
    lpr <- as(mmr %*% parm + offr, "vector")
    mmr[yi == max(yi),] <- 0
    mml[yi == min(yi),] <- 0

    Fr <- plogis(lpr)
    Fl <- plogis(lpl)
    fr <- dlogis(lpr)
    fl <- dlogis(lpl)
    fp <- function(z) dlogis(z) * (1 - 2 * plogis(z))
    dfr <- fp(lpr)
    dfl <- fp(lpl)

    Frl <- Fr - Fl
    w1 <- c(dfr / Frl) * w
    w2 <- c(dfl / Frl) * w
    w3 <- c(fr / Frl) * sqrt(w)
    w4 <- c(fl / Frl) * sqrt(w)
    W3 <- mmr * w3
    W4 <- mml * w4
    H <- -(crossprod(mmr * w1, mmr) - crossprod(mml * w2, mml) - 
           (crossprod(W3) - crossprod(W3, W4) - 
            crossprod(W4, W3) + crossprod(W4)))
    ret <- H[ncol(H), ncol(H)] - 
      H[ncol(H),-ncol(H)] %*% solve(H[-ncol(H), -ncol(H)],
                                    H[-ncol(H), ncol(H)])
    sqrt(1 / as(ret, "vector"))
  }

  parm1 <- c(plr$zeta, logOR)
  sd1 <- .sd(parm1, d$w)
  sd0 <- .sd(parm0, d$w)

  ## compute initial (too wide) Wald interval
  alpha_init <- (1 - level) / 5
  CI <- logOR + sd1 * qnorm(c(alpha_init, 1 - alpha_init))
  grd <- seq(from = CI[1], to = CI[2], length.out = 25)

  ## evaluate score statistic for hypothetial values of b = logOR
  prm <- parm1
  sc <- function(b) {
      ### polr is not needed here; this is a series of logistic
      ### regressions with offset => intercepts can be computed analytically
      m <- polr(zc ~ 1 + offset(Xgrp * b), start = prm[-length(prm)],
                data = d, weights = d$w, Hess = FALSE)
      prm <<- c(m$zeta, b)
      .sc(prm, d$w) * .sd(prm, d$w)
  }

  ### Note: qwilcox gives the distribution of the U statistic, this is
  ### decreasing
  Qa <- .qwilcox(c(1 - alpha, alpha), sizes = c(n.x, n.y), z = z, ...) 
  Qa <- Qa + n.x * (n.x + 1) / 2
  ### note: polr score under H_0 is (1 - 2 * R / N + 1 / N) * sd0
  Qa <- sd0 * (n.x - Qa * 2 / N + n.x / N)

  ## invert score statistic 
  grd_sc <- numeric(length(grd))
  for (i in 1:length(grd)) {
      grd_sc[i] <- sc(grd[i])
      if (!is.finite(grd_sc)[i]) break()
  }

  Sci <- NULL
  if (all(is.finite(grd_sc))) {
      if (all(diff(grd_sc) < 0) || all(diff(grd_sc) > 0)) {
          s <- spline(x = grd, y = grd_sc, method = "hyman")
          Sci <- approx(x = s$y, y = s$x, xout = Qa)$y
      }
  }

  ### use Taylor approximation if inversion failed
  if (is.null(Sci))
     Sci <- logOR + sd1 * Qa

  return(c(Sci, logOR))

}

