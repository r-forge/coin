
.AUC2logOR <- function(PI) {
    f <- function(x, PI) x + (exp(-x) * (PI + exp(2 * x) * 
        (PI - 1) + exp(x) * (1 - 2 * PI)))
    ret <- sapply(PI, function(p) uniroot(f, PI = p, interval = 50 * c(-1, 1))$root)
    return(ret)
}

Lehmann <- function(y, x, ...)
    UseMethod("Lehmann")

Lehmann.numeric <- function(y, x, nbins = 0, ...) {

    if (nbins) {
        breaks <- c(-Inf, quantile(y, prob = 1:nbins / (nbins + 1)), Inf)
    } else {
        breaks <- c(-Inf, sort(unique(y)))
    }
    r <- cut(y, breaks = breaks)
    Lehmann(r, x, ...)
}

Lehmann.factor <- function(y, x, type = c("OddsRatio", "HazardRatio", "Lehmann"), conf.level = .95, ...)
{

    tol <- sqrt(.Machine$double.eps)
    xrt <- table(y, x)
    xt1 <- xrt[,1]
    xt2 <- xrt[,2]

    type <- match.arg(type)

    if (type == "OddsRatio") {
        ### starting value
        xt <- table(x)
        W <- sum(rank(y)[which(x == levels(x)[1])])
        U <- prod(xt) + xt[1] * (xt[1] + 1) / 2 - W
        AUC <- U / prod(xt)
        start <- .AUC2logOR(AUC)

        F <- plogis
        Q <- qlogis
        f <- dlogis
        fp <- function(z) {
            p <- plogis(z)
            return(p * (1 - p)^2 - p^2 * (1 - p))
        }
    } else {
        stop("not yet implemented")
    }

    ll <- function(parm) {
        beta <- parm[1L]
        theta <- cumsum(parm[-1L])
        pltheta <- c(0, p <- F(theta))
        putheta <- c(p, 1)
        plxtheta <- c(0, p <- F(theta - beta))
        puxtheta <- c(p, 1)
        ret <- sum(xt1 * log(pmax(tol, putheta - pltheta))) + 
               sum(xt2 * log(pmax(tol, puxtheta - plxtheta)))
        -ret
    }

    sc <- function(parm) {
        beta <- parm[1L]
        theta <- cumsum(parm[-1L])
        ltheta <- c(-Inf, theta)
        utheta <- c(theta, Inf)
        f <- pmax(tol, F(utheta) - F(ltheta))
        fx <- pmax(tol, F(utheta - beta) - F(ltheta - beta))
        z <- xt1 * f(utheta) / f
        ret <- c(0, rev(cumsum(rev(z)[-1])))
        z <- xt1 * f(ltheta) / f
        ret <- ret - c(0, rev(cumsum(rev(z[-1]))))
        z <- xt2 * f(utheta - beta) / fx
        ret <- ret + c(-sum(z[-length(z)]), rev(cumsum(rev(z)[-1])))
        z <- xt2 * f(ltheta - beta) / fx
        ret <- ret - c(-sum(z), rev(cumsum(rev(z[-1]))))
        -ret
    }

    se <- function(parm) {
        beta <- parm[1L]
        theta <- cumsum(parm[-1L])
        ltheta <- c(-Inf, theta)
        utheta <- c(theta, Inf)

        f <- pmax(tol, F(utheta) - F(ltheta))
        fb <- pmax(tol, F(utheta - beta) - F(ltheta - beta))

        i1 <- length(utheta)
        i2 <- 1

        fu <- f(utheta)
        fub <- f(utheta - beta)
        fl <- f(ltheta)
        flb <- f(ltheta - beta)
        fpu <- fp(utheta)
        fpl <- fp(ltheta)
        fpub <- fp(utheta - beta)
        fplb <- fp(ltheta - beta)

        b <- -((xt1 * fu * fl / f^2)[-c(i2)] +
               (xt2 * fub * flb / fb^2)[-c(i2)])
        b <- b[-length(b)]
        x <- ((xt2 * fpub / fb)[-i1] - 
              (xt2 * fplb / fb)[-i2] - 
              ((xt2 * fub^2 / fb^2)[-i1] - 
               (xt2 * fub * flb / fb^2)[-i2] -
               (xt2 * fub * flb / fb^2)[-i1] +
               (xt2 * flb^2 / fb^2)[-i2]))
        a <- ((xt1 * fpu / f)[-i1] - 
              (xt1 * fpl / f)[-i2] - 
              ((xt1 * fu^2 / f^2)[-i1] + 
               (xt1 * fl^2 / f^2)[-i2]))
        a <- a + ((xt2 * fpub / fb)[-i1] - 
              (xt2 * fplb / fb)[-i2] - 
              ((xt2 * fub^2 / fb^2)[-i1] +
               (xt2 * flb^2 / fb^2)[-i2]))
        z <- -sum(xt2 * (fpub / fb - 
                         fplb / fb -
                         (fub^2 / fb^2 - 
                          2 * fub * flb / fb^2 +
                          flb^2 / fb^2)
                         )
                 )
        c(Schur_symtri(a = -a, b = b, X = x, Z = z))
    }

    ql <- Q(1:(length(xt1) - 1L) / length(xt1))
    theta <- c(start, ql[1], diff(ql))
    lwr <- c(-Inf, -Inf, rep(tol, length(theta) - 2))
    upr <- rep(Inf, length(theta))

    ret <- optim(par = theta, fn = ll, gr = sc, 
                 lower = lwr, upper = upr, 
                 method = "L-BFGS-B", ...)
    cf <- ret$par
    c(cf[1], cf[1] + sqrt(se(cf)) * qnorm(1 - (1 - conf.level) / 2) * c(-1, 1))
}
