
trafo.test <- function(y, x, ...)
    UseMethod("trafo.test")

trafo.test.numeric <- function(y, x, nbins = 0, ...) {

    uy <- unique(y)
    if (nbins && nbins < length(uy)) {
        breaks <- c(-Inf, quantile(y, prob = seq_len(nbins) / (nbins + 1L)), Inf)
    } else {
        breaks <- c(-Inf, sort(uy), Inf)
    }
    r <- cut(y, breaks = breaks)[, drop = TRUE]
    trafo.test(r, x, ...)
}

trafo.test.factor <- function(y, x, type = Wilcoxon(), mu = 0, conf.level = .95, 
                           Wald = FALSE, B = 0, ...)
{

    tol <- sqrt(.Machine$double.eps)
    stopifnot(nlevels(y <- y[,drop = TRUE]) > 1L)
    stopifnot(is.factor(x))
    stopifnot(nlevels(x <- x[,drop = TRUE]) == 2L)

    xrt <- table(y, x)
    xt1 <- xrt[,1]
    xt2 <- xrt[,2]

    ### starting value via AUC
    xt <- table(x)
    W <- sum(rank(y)[which(x == levels(x)[1])])
    U <- prod(xt) + xt[1] * (xt[1] + 1) / 2 - W
    AUC <- U / prod(xt)

    betastart <- type$PI2lp(AUC)
    F <- type$p
    Q <- type$q
    f <- type$d
    fp <- type$dd

    ll <- function(parm) {
        beta <- parm[1L]
        theta <- cumsum(parm[-1L])
        pltheta <- c(0, p <- F(theta))
        putheta <- c(p, 1)
        plxtheta <- c(0, p <- F(theta - mu - beta))
        puxtheta <- c(p, 1)
        ret <- sum(xt1 * log(pmax(tol, putheta - pltheta))) + 
               sum(xt2 * log(pmax(tol, puxtheta - plxtheta)))
        -ret
    }

    profile <- function(beta, parm_start, lwr, upr) {
        bll <- function(parm)
            ll(c(beta, parm))
        bsc <- function(parm)
            sc(c(beta, parm))[-1L]
        ### <TH> try() and convergence? </TH>
        ret <- optim(par = parm_start, fn = bll, gr = bsc, 
                     lower = lwr, upper = upr, 
                     method = "L-BFGS-B", ...)
        ret$par
    }

    resid <- function(beta, theta) {
        ltheta <- c(-Inf, theta)
        utheta <- c(theta, Inf)
        Ful <- pmax(tol, F(utheta) - F(ltheta))
        Fulb <- pmax(tol, F(utheta - beta) - F(ltheta - beta))
        ret <- xt1 * (f(utheta) - f(ltheta)) / Ful + 
               xt2 * (f(utheta - beta) - f(ltheta - beta)) / Fulb
        ret <- ret / (xt1 + xt2)
        ret
    }

    sc <- function(parm) {
        beta <- parm[1L]
        theta <- cumsum(parm[-1L])
        ltheta <- c(-Inf, theta)
        utheta <- c(theta, Inf)
        Ful <- pmax(tol, F(utheta) - F(ltheta))
        Fulb <- pmax(tol, F(utheta - mu - beta) - F(ltheta - mu - beta))
        z <- xt1 * f(utheta) / Ful
        ret <- c(0, rev(cumsum(rev(z)[-1])))
        z <- xt1 * f(ltheta) / Ful
        ret <- ret - c(0, rev(cumsum(rev(z[-1]))))
        z <- xt2 * f(utheta - mu - beta) / Fulb
        ret <- ret + c(-sum(z[-length(z)]), rev(cumsum(rev(z)[-1])))
        z <- xt2 * f(ltheta - mu - beta) / Fulb
        ret <- ret - c(-sum(z), rev(cumsum(rev(z[-1]))))
        -ret
    }

    he <- function(parm) {
        beta <- parm[1L]
        theta <- cumsum(parm[-1L])
        ltheta <- c(-Inf, theta)
        utheta <- c(theta, Inf)

        Ful <- pmax(tol, F(utheta) - F(ltheta))
        Fulb <- pmax(tol, F(utheta - mu - beta) - F(ltheta - mu - beta))

        i1 <- length(utheta)
        i2 <- 1

        fu <- f(utheta)
        fub <- f(utheta - mu - beta)
        fl <- f(ltheta)
        flb <- f(ltheta - mu - beta)
        fpu <- fp(utheta)
        fpl <- fp(ltheta)
        fpub <- fp(utheta - mu - beta)
        fplb <- fp(ltheta - mu - beta)

        b <- -((xt1 * fu * fl / Ful^2)[-c(i2)] +
               (xt2 * fub * flb / Fulb^2)[-c(i2)])
        b <- b[-length(b)]
        x <- ((xt2 * fpub / Fulb)[-i1] - 
              (xt2 * fplb / Fulb)[-i2] - 
              ((xt2 * fub^2 / Fulb^2)[-i1] - 
               (xt2 * fub * flb / Fulb^2)[-i2] -
               (xt2 * fub * flb / Fulb^2)[-i1] +
               (xt2 * flb^2 / Fulb^2)[-i2]
              )
             )
        a <- ((xt1 * fpu / Ful)[-i1] - 
              (xt1 * fpl / Ful)[-i2] - 
              ((xt1 * fu^2 / Ful^2)[-i1] + 
               (xt1 * fl^2 / Ful^2)[-i2]
              )
             )
        a <- a + ((xt2 * fpub / Fulb)[-i1] - 
                  (xt2 * fplb / Fulb)[-i2] - 
                  ((xt2 * fub^2 / Fulb^2)[-i1] +
                   (xt2 * flb^2 / Fulb^2)[-i2]
                  )
                 )
        z <- -sum(xt2 * (fpub / Fulb - 
                         fplb / Fulb -
                         (fub^2 / Fulb^2 - 
                          2 * fub * flb / Fulb^2 +
                          flb^2 / Fulb^2
                         )
                        )
                 )
        c(Schur_symtri(a = -a, b = b, X = x, Z = z))
    }

    ### ECDF
    cs <- cumsum(xt1 + xt2)
    ql <- Q(cs[-length(cs)] / cs[length(cs)])
    parm_start <- c(betastart, ql[1], diff(ql))
    lwr <- c(-Inf, -Inf, rep(tol, length(parm_start) - 2))
    upr <- rep(Inf, length(parm_start))
    ### <TH> check & convergence </TH>
    ret <- optim(par = parm_start, fn = ll, gr = sc, 
                 lower = lwr, upper = upr, 
                 method = "L-BFGS-B", ...)
    cf <- ret$par

    alpha <- (1 - conf.level) / 2
    aW <- alpha
    if (!Wald) aW <- alpha / 4
    WALD <- c(cf[1], cf[1] + sqrt(1 / he(cf)) * qnorm(1 - aW) * c(-1, 1))
    if (Wald) return(WALD)
    
    ### <TH> exact: as part of family-style argument ?</TH>

    if (B) {
        if (mu == 0) {
            res <- resid(0, ql)
            pstart <- parm_start[-1L]
        } else {
            res <- resid(mu, pstart <- profile(0, start = cf[-1L], lwr = lwr[-1L], upr = upr[-1L]))
        }
        rt <- r2dtable(B, r = xt1 + xt2, c = c(sum(xt1), sum(xt2)))
        se0 <- sqrt(1 / he(c(0, pstart)))
        U <- sapply(rt, function(x) sum(x[,1] * res)) * se0
        qz <- quantile(U, probs = c(alpha, 1 - alpha))
        ### <TH> achieved alpha ? </TH>
    } else {
        ### score
        qz <- qnorm(c(alpha, 1 - alpha))
    }
    pstart <- cf[-1L]
    sf <- function(b) {
        bparm <- profile(b, parm_start = pstart, lwr = lwr[-1L], upr = upr[-1L])
        pstart <<- bparm
        sc(c(b, bparm))[1L] * sqrt(1 / he(c(b, bparm)))
    }
    grd <- c(WALD[2], WALD[1])
    lci <- uniroot(function(b) sf(b) - qz[1], interval = grd)$root
    grd <- c(WALD[1], WALD[3])
    uci <- uniroot(function(b) sf(b) - qz[2], interval = grd)$root
    c(cf[1], lci, uci)
}
