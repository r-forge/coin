
trafo.test <- function(y, ...)
    UseMethod("trafo.test")

trafo.test.formula <- function(formula, data, weights, subset, na.action = na.pass, ...)
{
    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    if (length(attr(terms(formula[-2L]), "term.labels")) != 1L)
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ") # works in all cases
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    w <- as.vector(model.weights(mf))
    y <- mf[[response]]
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
    ## Call the default method.
    RVAL <- trafo.test(y = y, x = g, weights = w, ...)
    RVAL$data.name <- DNAME
    RVAL
}
 
trafo.test.numeric <- function(y, x, weights = NULL, nbins = 0, ...) {

    stopifnot(is.null(weights))
    DNAME <- paste(deparse1(substitute(x)), "and",
                   deparse1(substitute(y)))
    uy <- unique(y)
    if (nbins && nbins < length(uy)) {
        nbins <- ceiling(nbins)
        breaks <- c(-Inf, quantile(y, prob = seq_len(nbins) / (nbins + 1L)), Inf)
    } else {
        breaks <- c(-Inf, sort(uy), Inf)
    }
    r <- cut(y, breaks = breaks)[, drop = TRUE]
    RVAL <- trafo.test(r, x, ...)
    RVAL$data.name <- DNAME
    RVAL$method <- gsub("Ordinal|Binary", "Semiparametric", RVAL$method)
    RVAL
}

trafo.test.factor <- function(y, x, weights = 1, ...) {

    DNAME <- paste(deparse1(substitute(x)), "and",
                   deparse1(substitute(y)))

    stopifnot(is.factor(x))
    d <- data.frame(y = y, x = x, w = weights)
    tab <- xtabs(w ~ y + x, data = d)
    RVAL <- trafo.test(tab, ...)
    RVAL$data.name <- DNAME
    RVAL
}

trafo.test.table <- function(y,
                             link = c("logit", "probit", "cloglog", "loglog"), 
                             alternative = c("two.sided", "less", "greater"),
                             inference = c("Wald", "LRatio", "MLScore", "PermScore"),
                             parmscale = c("shift", "AUC/PI", "Overlap"),
                             mu = 0, conf.level = .95, B = 0, delta = NULL, ...)
{

    d <- dim(y)
    stopifnot(length(d) == 2 && d[2] == 2)
    dn <- dimnames(y)
    DNAME <- NULL
    if (!is.null(dn))
        DNAME <- paste(names(dn)[1], "and", names(dn)[2])

    y <- y[rowSums(y) > 0,,drop = FALSE]
    
    xt1 <- y[,1]
    xt2 <- y[,2]
    mm1 <- which(xt1 > 0)
    mm1 <- mm1[c(1, length(mm1))]
    mm2 <- which(xt2 > 0)
    mm2 <- mm2[c(1, length(mm2))]
    if (mm1[1] > mm2[2] || mm2[1] > mm1[2])
        warning("Data fully separated, results instable")
    xt <- c(xt1, xt2)
    stopifnot(sum(xt1 + xt2 > 0) > 1L)
    stopifnot(sum(xt1) > 0 && sum(xt2) > 0)

    tol <- sqrt(.Machine$double.eps)

    inference <- match.arg(inference)
    alternative <- match.arg(alternative)
    parmscale <- match.arg(parmscale)
    if (parmscale != "shift") {
        stopifnot(inference != "Wald")
        stopifnot(mu == 0)
    }

    if (!inherits(link, "linkfun")) {
        link <- match.arg(link)
        link <- do.call(link, list())
    }
    names(mu) <- link$parm

    F <- function(x) .p(link, x)
    Q <- function(p) .q(link, p)
    f <- function(x) .d(link, x)
    fp <- function(x) .dd(link, x)

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
        ret <- Schur_symtri(a = -a, b = b, X = x, Z = z)
#        while (ret < tol) {
#            a <- a - sqrt(tol)
#            z <- z + sqrt(tol)
#            ret <- Schur_symtri(a = -a, b = b, X = x, Z = z)
#        }
        return(c(ret))
    }

    ### log-OR from 2x2 table as starting value
    cs <- cumsum(xt1 + xt2) 
    ymed <- min(c(which.min(abs(cs / cs[length(cs)] - .5)), length(xt1) - 1))
    tab2x2 <- cbind(colSums(y[1:ymed,,drop = FALSE]), colSums(y[-(1:ymed),,drop = FALSE]))
    betastart <- sum(log(diag(tab2x2))) - log(tab2x2[2,1]) - log(tab2x2[1,2])
    if (!is.finite(betastart)) betastart <- 0

    ### ECDF
    cs <- cumsum(xt1 + xt2)
    ql <- Q(cs[-length(cs)] / cs[length(cs)])
    if (is.null(delta)) {
        parm_start <- c(betastart, ql[1], diff(ql))
        lwr <- c(-Inf, -Inf, rep(tol, length(parm_start) - 2))
        upr <- rep(Inf, length(parm_start))
        ret <- try(optim(par = parm_start, fn = ll, gr = sc, 
                         lower = lwr, upper = upr, 
                         method = "L-BFGS-B", hessian = FALSE, ...))
    } else {
        stopifnot(inference == "Wald")
        llf <- function(parm) ll(c(delta, parm))
        scf <- function(parm) sc(c(delta, parm))[-1L]
        parm_start <- c(ql[1], diff(ql))
        lwr <- c(-Inf, rep(tol, length(parm_start) - 1))
        upr <- rep(Inf, length(parm_start) - 1)
        ret <- try(optim(par = parm_start, fn = llf, gr = scf, 
                         lower = lwr, upper = upr, 
                         method = "L-BFGS-B", hessian = FALSE, ...))
    }

    if (inherits(ret, "try-error"))
        stop("optimisation failed")
    if (ret$convergence > 1)
        warning(ret$message)

    cf <- ret$par
    if (!is.null(delta)) cf <- c(delta, cf)
    ESTIMATE <- cf[1]
    names(ESTIMATE) <- paste(link$parm, ifelse(mu == 0, "", paste0("-", mu)))
    HE <- try(he(cf))
    if (inherits(ret, "try-error"))
        stop("Computation of Hessian failed")

    if (inference == "Wald") {
        STATISTIC <- c("Wald Z" = unname(ESTIMATE * sqrt(HE)))
        if (alternative == "less") {
            PVAL <- pnorm(STATISTIC)
            cint <- c(-Inf, ESTIMATE + qnorm(conf.level) / sqrt(HE))
        } else if (alternative == "greater") {
            PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
            cint <- c(ESTIMATE - qnorm(conf.level) / sqrt(HE), Inf)
        } else {
            PVAL <- 2 * pnorm(-abs(STATISTIC))
            cint <- ESTIMATE + c(-1, 1) * qnorm(1 - (1 - conf.level) / 2) / sqrt(HE)
        }
        attr(cint, "conf.level") <- conf.level
        TYPE <- "Wald"
    } else {
        alpha <- (1 - conf.level)
        if (alternative == "two.sided") alpha <- alpha / 2
        aW <- alpha / 10
        WALD <- c(ESTIMATE, ESTIMATE + sqrt(1 / HE) * qnorm(1 - aW) * c(-1, 1))

        if (inference == "LRatio") {
            stopifnot(alternative == "two.sided")
            qc <- qchisq(conf.level, df = 1)
            pstart <- cf[-1L]
            logLRstat <- function(b) {
                bparm <- profile(b, parm_start = pstart, lwr = lwr[-1L], upr = upr[-1L])
                pstart <<- bparm
                -2 * (ret$value - ll(c(b, bparm)))
            }
            STATISTIC <- c("LR Chisq" = unname(logLRstat(0)))
            PVAL <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)

            lci <- -Inf
            grd <- c(WALD[2], WALD[1])
            if (alternative != "greater") {
                grd <- c(WALD[2], WALD[1])
                lci <- uniroot(function(b) logLRstat(b) - qc, interval = grd, extendInt = "yes")$root
            }
            uci <- Inf
            if (alternative != "less") {
                grd <- c(WALD[1], WALD[3])
                uci <- uniroot(function(b) logLRstat(b) - qc, interval = grd, extendInt = "yes")$root
            }
            cint <- c(lci, uci)
            attr(cint, "conf.level") <- conf.level
            TYPE <- "LR"
        } else {

            pstart <- cf[-1L]
            sf <- function(b) {
                bparm <- profile(b, parm_start = pstart, lwr = lwr[-1L], upr = upr[-1L])
                pstart <<- bparm
                sc(c(b, bparm))[1L] * sqrt(1 / he(c(b, bparm)))
            }
            STATISTIC <- c("Score Z" = unname(sf(0)))

            if (inference == "MLScore") {
                ### score
                qz <- qnorm(c(alpha, 1 - alpha))
                if (alternative == "less") {
                    PVAL <- pnorm(STATISTIC)
                } else if (alternative == "greater") {
                    PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
                } else {
                    PVAL <- 2 * pnorm(-abs(STATISTIC))
                }
                TYPE <- "score"
            } else { ### PermScore
                if (mu == 0) {
                    res <- resid(0, ql)
                    pstart <- parm_start[-1L]
                } else {
                    res <- resid(mu, cumsum(pstart <- profile(0, parm_start = cf[-1L], lwr = lwr[-1L], upr = upr[-1L])))
                }
                se0 <- sqrt(1 / he(c(0, pstart)))
                sp <- statpvalPerm(res = res * se0, xt = y,
                                   alternative = alternative, B = B)
                STATISTIC <- sp$STATISTIC
                PVAL <- sp$PVAL
                qz <- qPerm(p = c(alpha, 1 - alpha), res = res * se0, xt = y, B = B)
                ### <TH> achieved alpha ? </TH>
                TYPE <- paste(sp$TYPE, "permutation")
           }
            lci <- -Inf
            grd <- c(WALD[2], WALD[1])
            if (alternative != "greater") {
                grd <- c(WALD[2], WALD[1])
                lci <- uniroot(function(b) sf(b) - qz[1], interval = grd, extendInt = "yes")$root
            }
            uci <- Inf
            if (alternative != "less") {
                grd <- c(WALD[1], WALD[3])
                uci <- uniroot(function(b) sf(b) - qz[2], interval = grd, extendInt = "yes")$root
           }
           cint <- c(lci, uci)
           attr(cint, "conf.level") <- conf.level
        }
    }
    METHOD <- paste(ifelse(nlevels(y) > 2, "Ordinal", "Binary"), 
                    "two-sample", TYPE, "inference for",
                    link$model, "models")
    if (parmscale != "shift") {
        if (parmscale == "AUC/PI") {
            ESTIMATE <- link$parm2PI(ESTIMATE)
            names(ESTIMATE) <- "AUC/PI"
            cint <- link$parm2PI(cint)
        }
        if (parmscale == "Overlap") {
            ESTIMATE <- link$parm2OVL(ESTIMATE)
            names(ESTIMATE) <- "Overlap coefficient"
            cint <- link$parm2OVL(sort(abs(cint), decreasing = TRUE))
        }
    }
    RVAL <- list(statistic = STATISTIC, parameter = NULL, p.value = as.numeric(PVAL), 
                 null.value = mu, alternative = alternative, method = METHOD, 
                 data.name = DNAME, conf.int = cint, estimate = ESTIMATE,
                 logLik = ret$value, hessian = HE)
    RVAL$link <- link
    class(RVAL) <- c("trafo.test", "htest")
    return(RVAL)
}

power.trafo.test <- function(n = NULL, prob = NULL, aratio = 1, delta = NULL, sig.level = .05, power = NULL,
                             link = c("logit", "probit", "cloglog", "loglog"),
                             alternative = c("two.sided", "less", "greater"), ...) 
{

    if (sum(vapply(list(n, delta, power, sig.level), is.null, 
        NA)) != 1) 
        stop("exactly one of 'n', 'delta', 'power', and 'sig.level' must be NULL")
    stats:::assert_NULL_or_prob(sig.level)
    stats:::assert_NULL_or_prob(power)

    tol <- sqrt(.Machine$double.eps)

    if (is.null(n)) 
        n <- uniroot(function(n)
             power.trafo.test(n = n, prob = prob, aratio = aratio, delta = delta, 
                             sig.level = sig.level, link = link, alternative = alternative, 
                             ...)$power
            - power, c(5, 1e+03), 
            tol = tol, extendInt = "upX")$root
    else if (is.null(delta)) 
        delta <- uniroot(function(delta) 
            power.trafo.test(n = n, prob = prob, aratio = aratio, delta = delta, 
                             sig.level = sig.level, link = link, alternative = alternative, 
                             ...)$power
            - power, 
    ### <TH> interval depending on alternative, symmetry? </TH>
            c(0, 10), tol = tol, extendInt = "upX")$root
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) 
            power.trafo.test(n = n, prob = prob, aratio = aratio, delta = delta, 
                             sig.level = sig.level, link = link, alternative = alternative, 
                             ...)$power
            - power, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
    else if (is.null(power)) {

        if (!inherits(link, "linkfun")) {
            link <- match.arg(link)
            link <- do.call(link, list())
        }

        if (is.null(prob)) prob <- rep(1 / n, n)
        p0 <- cumsum(prob)
        h0 <- .q(link, p0)
        h1 <- h0 - delta
        p1 <- .p(link, h1)
        se <- rep(NA, 100)

        for (i in 1:length(se)) {
            n0 <- rmultinom(1, size = n, prob = prob)
            n1 <- rmultinom(1, size = aratio * n, prob = c(p1[1], diff(p1)))
            y <- as.table(cbind(n0, n1))
            suppressWarnings(HE <- try(trafo.test(y, delta = delta, link = link,
                                 ...)$hessian, silent = TRUE))
            if (!inherits(HE, "try-error"))
                se[i] <- 1 / sqrt(HE)
        }
        se <- mean(se, na.rm = TRUE)
        if (is.na(se)) stop("approximating hessian failed")
        alternative <- match.arg(alternative)
        power  <- switch(alternative, 
            "two.sided" = pnorm(qnorm(sig.level / 2) + delta / se) + 
                          pnorm(qnorm(sig.level / 2) - delta / se),
            "less" = pnorm(qnorm(sig.level) - delta / se),
            "greater" = pnorm(qnorm(sig.level) + delta / se)
        )
    }
    list(power = power, n = n, delta = delta, sig.level = sig.level)
}
