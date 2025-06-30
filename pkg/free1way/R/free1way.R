
# ML estimation

.free1wayML <- function(x, link, mu = 0, start = NULL, fix = NULL, 
                        residuals = TRUE, score = TRUE, hessian = TRUE, 
                        tol = sqrt(.Machine$double.eps), ...) {

    stopifnot(is.table(x))
    dx <- dim(x)
    dn <- dimnames(x)
    if (length(dx) == 2L) {
        x <- as.table(array(c(x), dim = dx <- c(dx, 1L)))
        dimnames(x) <- dn <- c(dn, list(A = "A"))
    }

    ms <- c(list(x), lapply(seq_along(dx), function(j) marginSums(x, j) > 0))
    ms$drop <- FALSE
    x <- do.call("[", ms)

    dx <- dim(x)
    dn <- dimnames(x)
    stopifnot(length(dx) >= 3L)
    stopifnot(dx[1L] > 1L)
    K <- dx[2L]
    stopifnot(K > 1L)

    F <- function(q) .p(link, q = q)
    Q <- function(p) .q(link, p = p)
    f <- function(q) .d(link, x = q)
    fp <- function(q) .dd(link, x = q)

    # setup
    
    # table2list

    .table2list <- function(x) {

        dx <- dim(x)
        if (length(dx) == 1L)
            stop("")
        if (length(dx) == 2L)
            x <- as.table(array(x, dim = c(dx, 1)))
        ms <- c(list(x), lapply(seq_along(dx), function(j) marginSums(x, j) > 0))
        ms$drop <- FALSE
        x <- do.call("[", ms)
        dx <- dim(x)
        stopifnot(length(dx) >= 3L)
        K <- dim(x)[2L]
        B <- dim(x)[3L]
        stopifnot(dx[1L] > 1L)
        stopifnot(K > 1L)

        xrc <- NULL
        if (length(dx) == 4L) {
            if (dx[4] == 2L) {
                xrc <- array(x[,,,"FALSE", drop = TRUE], dim = dx[1:3])
                x <- array(x[,,,"TRUE", drop = TRUE], dim = dx[1:3])
            } else {
                stop("")
            }
        }

        xlist <- xrclist <- vector(mode = "list", length = B)

        lwr <- rep(-Inf, times = K - 1)
        for (b in seq_len(B)) {
            xb <- matrix(x[,,b, drop = TRUE], ncol = K)
            xw <- rowSums(abs(xb)) > 0
            ### do not remove last parameter if there are corresponding
            ### right-censored observations
            if (!is.null(xrc) && any(xrc[dx[1],,b,drop = TRUE] > 0))
                xw[length(xw)] <- TRUE
            if (sum(xw) > 1L) {
                xlist[[b]] <- xb[xw,,drop = FALSE]
                attr(xlist[[b]], "idx") <- xw
                if (!is.null(xrc)) {
                    xrclist[[b]] <- matrix(xrc[xw,,b,drop = TRUE], ncol = K)
                    attr(xrclist[[b]], "idx") <- xw
                }
            }
        }
        nn <- !sapply(xlist, is.null)
        ret <- list(xlist = xlist[nn])
        if (!is.null(xrc))
            ret$xrclist <- xrclist[nn]
        ret$strata <- nn
        ret
    }
    
    C <- dim(x)[1L]
    K <- dim(x)[2L]
    B <- dim(x)[3L]
    xl <- .table2list(x)
    xlist <- xl$xlist
    xrclist <- xl$xrclist
    if (NS <- is.null(start))
        start <- rep.int(0, K - 1)
    lwr <- rep(-Inf, times = K - 1)
    for (b in seq_len(length(xlist))) {
        lwr <- c(lwr, -Inf, rep.int(tol, times = nrow(xlist[[b]]) - 2L))
        if (NS) {
            ecdf0 <- cumsum(rowSums(xlist[[b]]))
            ecdf0 <- ecdf0[-length(ecdf0)] / ecdf0[length(ecdf0)]
            Qecdf <- Q(ecdf0)
            start <- c(start, Qecdf[1], diff(Qecdf))
            start[!is.finite(start)] <- 0
        }
    }
    
    # cumsumrev
    
    .rcr <- function(z)
        # Reduce('+', z, accumulate = TRUE, right = TRUE)
        rev.default(cumsum(rev.default(z)))
    
    # negative logLik
    
    .nll <- function(parm, x, mu = 0, rightcensored = FALSE) {
        # parm to prob
        
        bidx <- seq_len(ncol(x) - 1L)
        beta <- c(0, mu + parm[bidx])
        intercepts <- c(-Inf, cumsum(parm[- bidx]), Inf)
        tmb <- intercepts - matrix(beta, nrow = length(intercepts),  
                                         ncol = ncol(x),
                                         byrow = TRUE)
        Ftmb <- F(tmb)
        if (rightcensored) {
            prb <- pmax(1 - Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
        } else {
            prb <- pmax(Ftmb[- 1L, , drop = FALSE] - 
                        Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
        } 
        
        return(- sum(x * log(prb)))
    }
    
    # negative score
    
    .nsc <- function(parm, x, mu = 0, rightcensored = FALSE) {
        # parm to prob
        
        bidx <- seq_len(ncol(x) - 1L)
        beta <- c(0, mu + parm[bidx])
        intercepts <- c(-Inf, cumsum(parm[- bidx]), Inf)
        tmb <- intercepts - matrix(beta, nrow = length(intercepts),  
                                         ncol = ncol(x),
                                         byrow = TRUE)
        Ftmb <- F(tmb)
        if (rightcensored) {
            prb <- pmax(1 - Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
        } else {
            prb <- pmax(Ftmb[- 1L, , drop = FALSE] - 
                        Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
        } 
        
        # density prob ratio
        
        ftmb <- f(tmb)
        zu <- x * ftmb[- 1, , drop = FALSE] / prb
        if (rightcensored) zu[] <- 0 ### derivative of a constant
        zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
        

        ret <- numeric(length(parm))
        ret[bidx] <- colSums(zl)[-1L] -
                     colSums(zu[-nrow(zu),,drop = FALSE])[-1L]
        ret[-bidx] <- Reduce("+", 
                             lapply(1:ncol(x), 
                                 function(j) {
                                     .rcr(zu[-nrow(zu),j]) - 
                                     .rcr(zl[-1,j])
                                 })
                             )
        - ret
    }
    
    # negative score residuals
    
    .nsr <- function(parm, x, mu = 0, rightcensored = FALSE) {
        # parm to prob
        
        bidx <- seq_len(ncol(x) - 1L)
        beta <- c(0, mu + parm[bidx])
        intercepts <- c(-Inf, cumsum(parm[- bidx]), Inf)
        tmb <- intercepts - matrix(beta, nrow = length(intercepts),  
                                         ncol = ncol(x),
                                         byrow = TRUE)
        Ftmb <- F(tmb)
        if (rightcensored) {
            prb <- pmax(1 - Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
        } else {
            prb <- pmax(Ftmb[- 1L, , drop = FALSE] - 
                        Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
        } 
        
        # density prob ratio
        
        ftmb <- f(tmb)
        zu <- x * ftmb[- 1, , drop = FALSE] / prb
        if (rightcensored) zu[] <- 0 ### derivative of a constant
        zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
        

        ret <- rowSums(zu - zl) / rowSums(x)
        ret[!is.finite(ret)] <- 0
        ret
    }
    
    # Hessian
    
    .hes <- function(parm, x, mu = 0, rightcensored = FALSE) {
        # parm to prob
        
        bidx <- seq_len(ncol(x) - 1L)
        beta <- c(0, mu + parm[bidx])
        intercepts <- c(-Inf, cumsum(parm[- bidx]), Inf)
        tmb <- intercepts - matrix(beta, nrow = length(intercepts),  
                                         ncol = ncol(x),
                                         byrow = TRUE)
        Ftmb <- F(tmb)
        if (rightcensored) {
            prb <- pmax(1 - Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
        } else {
            prb <- pmax(Ftmb[- 1L, , drop = FALSE] - 
                        Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
        } 
        

        # Hessian prep
        
        ftmb <- f(tmb)
        fptmb <- fp(tmb)

        dl <- ftmb[- nrow(ftmb), , drop = FALSE]
        du <- ftmb[- 1, , drop = FALSE]
        if (rightcensored) du[] <- 0
        dpl <- fptmb[- nrow(ftmb), , drop = FALSE]
        dpu <- fptmb[- 1, , drop = FALSE]
        if (rightcensored) dpu[] <- 0
        dlm1 <- dl[,-1L, drop = FALSE]
        dum1 <- du[,-1L, drop = FALSE]
        dplm1 <- dpl[,-1L, drop = FALSE]
        dpum1 <- dpu[,-1L, drop = FALSE]
        prbm1 <- prb[,-1L, drop = FALSE]

        i1 <- length(intercepts) - 1L
        i2 <- 1L
        

        # off-diagonal elements for Hessian of intercepts
        
        Aoffdiag <- -rowSums(x * du * dl / prb^2)[-i2]
        Aoffdiag <- Aoffdiag[-length(Aoffdiag)]
        
        # diagonal elements for Hessian of intercepts
        
        Adiag <- -rowSums((x * dpu / prb)[-i1,,drop = FALSE] - 
                          (x * dpl / prb)[-i2,,drop = FALSE] - 
                          ((x * du^2 / prb^2)[-i1,,drop = FALSE] + 
                           (x * dl^2 / prb^2)[-i2,,drop = FALSE]
                          )
                         )
                          
        
        # intercept / shift contributions to Hessian
        
        xm1 <- x[,-1L,drop = FALSE] 
        X <- ((xm1 * dpum1 / prbm1)[-i1,,drop = FALSE] - 
              (xm1 * dplm1 / prbm1)[-i2,,drop = FALSE] - 
              ((xm1 * dum1^2 / prbm1^2)[-i1,,drop = FALSE] - 
               (xm1 * dum1 * dlm1 / prbm1^2)[-i2,,drop = FALSE] -
               (xm1 * dum1 * dlm1 / prbm1^2)[-i1,,drop = FALSE] +
               (xm1 * dlm1^2 / prbm1^2)[-i2,,drop = FALSE]
              )
             )

        Z <- -colSums(xm1 * (dpum1 / prbm1 - 
                             dplm1 / prbm1 -
                             (dum1^2 / prbm1^2 - 
                              2 * dum1 * dlm1 / prbm1^2 +
                              dlm1^2 / prbm1^2
                             )
                            )
                     )
        if (length(Z) > 1L) Z <- diag(Z)
        

        if (length(Adiag) > 1L) {
            if(is.null(tryCatch(loadNamespace("Matrix"), error = function(e)NULL)))
                    stop(gettextf("%s needs package 'Matrix' correctly installed",
                                  "free1way.test"),
                         domain = NA)
            A <- Matrix::bandSparse(length(Adiag), k = 0:1, diagonals = list(Adiag, Aoffdiag), 
                                    symmetric = TRUE)
        } else {
            A <- matrix(Adiag)
        }
        return(list(A = A, X = X, Z = Z))
    }
    
    # stratified negative logLik
    
    .snll <- function(parm, x, mu = 0, rightcensored = FALSE) {
        # stratum prep
        
        C <- sapply(x, NROW) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        beta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        
        ret <- 0
        for (b in seq_len(B))
            ret <- ret + .nll(c(beta, intercepts[[b]]), x[[b]], mu = mu,
                              rightcensored = rightcensored)
        return(ret)
    }
    
    # stratified negative score
    
    .snsc <- function(parm, x, mu = 0, rightcensored = FALSE) {
        # stratum prep
        
        C <- sapply(x, NROW) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        beta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        
        ret <- numeric(length(bidx))
        for (b in seq_len(B)) {
            nsc <- .nsc(c(beta, intercepts[[b]]), x[[b]], mu = mu,
                        rightcensored = rightcensored)
            ret[bidx] <- ret[bidx] + nsc[bidx]
            ret <- c(ret, nsc[-bidx])
        }
        return(ret)
    }
    
    # stratified Hessian
    
    .shes <- function(parm, x, mu = 0, xrc = NULL) {
        # stratum prep
        
        C <- sapply(x, NROW) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        beta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        
        ret <- matrix(0, nrow = length(bidx), ncol = length(bidx))
        for (b in seq_len(B)) {
            H <- .hes(c(beta, intercepts[[b]]), x[[b]], mu = mu)
            if (!is.null(xrc)) {
                Hrc <- .hes(c(beta, intercepts[[b]]), xrc[[b]], mu = mu, 
                            rightcensored = TRUE)
                H$X <- H$X + Hrc$X
                H$A <- H$A + Hrc$A
                H$Z <- H$Z + Hrc$Z
            }
            ret <- ret + (H$Z - crossprod(H$X, solve(H$A, H$X)))
        }
        ret
    }
    
    # stratified negative score residual
    
    .snsr <- function(parm, x, mu = 0, rightcensored = FALSE) {
        # stratum prep
        
        C <- sapply(x, NROW) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        beta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        
        ret <- c()
        for (b in seq_len(B)) {
            idx <- attr(x[[b]], "idx")
            sr <- numeric(length(idx))
            sr[idx] <- .nsr(c(beta, intercepts[[b]]), x[[b]], mu = mu,
                            rightcensored = rightcensored)
            ret <- c(ret, sr)
        }
        return(ret)
    }
    
    # profile
    
    .profile <- function(start, fix = seq_len(K - 1)) {
        stopifnot(all(fix %in% seq_len(K - 1)))
        beta <- start[fix]
        ret <- optim(par = start[-fix], fn = function(par) {
                         p <- numeric(length(par) + length(fix))
                         p[fix] <- beta
                         p[-fix] <- par
                         ret <- .snll(p, x = xlist, mu = mu)
                         if (!is.null(xrclist))
                             ret <- ret + .snll(p, x = xrclist, mu = mu, rightcensored = TRUE)
                         ret
                     },
                     gr = function(par) {
                         p <- numeric(length(par) + length(fix))
                         p[fix] <- beta
                         p[-fix] <- par
                         ret <- .snsc(p, x = xlist, mu = mu)[-fix]
                         if (!is.null(xrclist))
                             ret <- ret + .snsc(p, x = xrclist, mu = mu, rightcensored = TRUE)[-fix]
                         ret
                     },
                     lower = lwr[-fix], method = "L-BFGS-B", 
                     hessian = FALSE, ...)
        p <- numeric(length(start))
        p[fix] <- beta
        p[-fix] <- ret$par
        ret$par <- p
        ret
    }
    
    # optim
    
    if (!length(fix)) {
        ret <- optim(par = start, 
                     fn = function(par) {
                         ret <- .snll(par, x = xlist, mu = mu)
                         if (!is.null(xrclist))
                             ret <- ret + .snll(par, x = xrclist, mu = mu, rightcensored = TRUE)
                         ret
                     },
                     gr = function(par) {
                         ret <- .snsc(par, x = xlist, mu = mu)
                         if (!is.null(xrclist))
                             ret <- ret + .snsc(par, x = xrclist, mu = mu, rightcensored = TRUE)
                         ret
                     },
                     lower = lwr, method = "L-BFGS-B", 
                     hessian = FALSE, ...)
    } else if (length(fix) == length(start)) {
        fn <- function(par) {
                         ret <- .snll(par, x = xlist, mu = mu)
                         if (!is.null(xrclist))
                             ret <- ret + .snll(par, x = xrclist, mu = mu, rightcensored = TRUE)
                         ret
                     }
        ret <- list(par = start, 
                    value = fn(start))
    } else {
        ret <- .profile(start, fix = fix)
    }
     
    # post processing
    
    if (is.null(fix) || (length(fix) == length(start)))
        parm <- seq_len(K - 1)
    else 
        parm <- fix
    if (any(parm >= K)) return(ret)

    ret$coefficients <- ret$par[parm]
    dn2 <- dimnames(x)[2L]
    names(ret$coefficients) <- cnames <- paste0(names(dn2), dn2[[1L]][1L + parm])

    if (score)
        ret$negscore <- .snsc(ret$par, x = xlist, mu = mu)[parm]
        if (!is.null(xrclist))
            ret$negscore <- ret$negscore + .snsc(ret$par, x = xrclist, mu = mu, rightcensored = TRUE)[parm]
    if (hessian) {
        ret$hessian <- .shes(ret$par, x = xlist, mu = mu, xrc = xrclist)
        if (length(parm) != nrow(ret$hessian))
           ret$hessian <- solve(ret$vcov <- solve(ret$hessian)[parm,parm])
        ret$vcov <- solve(ret$hessian)
        rownames(ret$vcov) <- colnames(ret$vcov) <- rownames(ret$hessian) <-
            colnames(ret$hessian) <-  cnames
    }
    if (residuals) {
        ret$residuals <- .snsr(ret$par, x = xlist, mu = mu)
        if (!is.null(xrclist)) {
            rcr <- .snsr(ret$par, x = xrclist, mu = mu, rightcensored = TRUE)
            ret$residuals <- c(rbind(matrix(ret$residuals, nrow = C),
                  matrix(rcr, nrow = C)))
         }
    }
    ret$profile <- function(start, fix)
        .free1wayML(x, link = link, mu = mu, start = start, fix = fix, tol = tol, ...) 

    ret$table <- x
    ret$mu <- mu
    ret$strata <- xl$strata
    names(ret$mu) <- link$parm
    

    class(ret) <- "free1wayML"
    ret
}

# free1way

free1way.test <- function(y, ...)
    UseMethod("free1way.test")
free1way.test.table <- function(y, link = c("logit", "probit", "cloglog", "loglog"), mu = 0, B = 0, ...)
{

    cl <- match.call()

    d <- dim(y)
    dn <- dimnames(y)
    DNAME <- NULL
    if (!is.null(dn)) {
        DNAME <- paste(names(dn)[1], "by", names(dn)[2], paste0("(", paste0(dn[2], collapse = ", "), ")"))
        if (length(dn) == 3L)
            DNAME <- paste(DNAME, "\n\t stratified by", names(dn)[3])
    }

    # link2fun
    
    if (!inherits(link, "linkfun")) {
        link <- match.arg(link)
        link <- do.call(link, list())
    }
    

    ret <- .free1wayML(y, link = link, mu = mu, ...)
    ret$link <- link
    ret$data.name <- DNAME
    ret$call <- cl

    alias <- link$alias
    if (length(link$alias) == 2L) alias <- alias[1L + (d[2] > 2L)]
    ret$method <- paste(ifelse(length(d) == 3L, "Stratified", ""), 
                        paste0(d[2], "-sample"), alias, 
                        "test against", link$model, "alternatives")

    cf <- ret$par
    cf[idx <- seq_len(d[2L] - 1L)] <- 0
    pr <- ret$profile(cf, idx)
    if (d[2L] == 2L)
        res <- pr$residuals / sqrt(c(pr$hessian))
    else
        res <- pr$residuals

    # Strasser Weber
    
    .SW <- function(res, xt) {

        if (length(dim(xt)) == 3L) {
            res <- matrix(res, nrow = dim(xt)[1L], ncol = dim(xt)[3])
            STAT <-  Exp <- Cov <- 0
            for (b in seq_len(dim(xt)[3L])) {
                sw <- .SW(res[,b, drop = TRUE], xt[,,b, drop = TRUE])
                STAT <- STAT + sw$Statistic
                Exp <- Exp + sw$Expectation
                Cov <- Cov + sw$Covariance
            }
            return(list(Statistic = STAT, Expectation = as.vector(Exp),
                        Covariance = Cov))
        }

        Y <- matrix(res, ncol = 1, nrow = length(xt))
        weights <- c(xt)
        x <- gl(ncol(xt), nrow(xt))
        X <- model.matrix(~ x, data = data.frame(x = x))[,-1L,drop = FALSE]

        w. <- sum(weights)
        wX <- weights * X
        wY <- weights * Y
        ExpX <- colSums(wX)
        ExpY <- colSums(wY) / w.
        CovX <- crossprod(X, wX)
        Yc <- t(t(Y) - ExpY)
        CovY <- crossprod(Yc, weights * Yc) / w.
        Exp <- kronecker(ExpY, ExpX)
        Cov <- w. / (w. - 1) * kronecker(CovY, CovX) -
               1 / (w. - 1) * kronecker(CovY, tcrossprod(ExpX))
        STAT <- crossprod(X, wY)
        list(Statistic = STAT, Expectation = as.vector(Exp),
             Covariance = Cov)
    }
    
    # resampling
    
    .resample <- function(res, xt, B = 10000) {

        if (length(dim(xt)) == 2L)
            xt <- as.table(array(xt, dim = c(dim(xt), 1)))

        res <- matrix(res, nrow = dim(xt)[1L], ncol = dim(xt)[3L])
        stat <- 0
        ret <- .SW(res, xt)
        if (dim(xt)[2L] == 2L) {
            ret$testStat <- c((ret$Statistic - ret$Expectation) / sqrt(c(ret$Covariance)))
        } else {
            ES <- t(ret$Statistic - ret$Expectation)
            ret$testStat <- sum(ES %*% solve(ret$Covariance) * ES)
        }
        ret$DF <- dim(xt)[2L] - 1L

        if (B) {
            for (j in 1:dim(xt)[3L]) {
               rt <- r2dtable(B, r = rowSums(xt[,,j]), c = colSums(xt[,,j]))
               stat <- stat + sapply(rt, function(x) colSums(x[,-1L, drop = FALSE] * res[,j]))
            }
            if (dim(xt)[2L] == 2L) {
                 ret$permStat <- (stat - ret$Expectation) / sqrt(c(ret$Covariance))
            } else {
                ES <- t(matrix(stat, ncol = B) - ret$Expectation)
                ret$permStat <- rowSums(ES %*% solve(ret$Covariance) * ES)
            }
        }
        ret
    }
    

    if (length(dim(y)) == 3L) y <- y[,,ret$strata, drop = FALSE]
    if (length(dim(y)) == 4L) {
        y <- y[,,ret$strata,, drop = FALSE]
        dy <- dim(y)
        dy[1] <- dy[1] * 2
        y <- apply(y, 3, function(x) rbind(x[,,2], x[,,1]))
        y <- array(y, dim = dy[1:3])
    }
    ret$perm <- .resample(res, y, B = B)

    if (!is.null(names(dn))) {
        fm <- as.formula(paste(names(dn)[1:2], collapse = "~"))
        ret$terms <- terms(fm, data = as.data.frame(y))
    }

    class(ret) <- "free1way"
    return(ret)
}

# free1way methods

coef.free1way <- function(object, what = c("shift", "PI", "AUC", "OVL"), ...)
{
    what <- match.arg(what)
    cf <- object$coefficients
    return(switch(what, "shift" = cf,
                        "PI" = object$link$parm2PI(cf),
                        "AUC" = object$link$parm2PI(cf),        ### same as PI
                        "OVL" = object$link$parm2OVL(cf)))
}
vcov.free1way <- function(object, ...)
    object$vcov
logLik.free1way <- function(object, ...)
    -object$value
### the next two could go into multcomp
model.frame.free1way <- function(formula, ...)
    as.data.frame(formula$table)
model.matrix.free1way <- function (object, ...) 
{
    mm <- model.matrix(delete.response(terms(object)), data = model.frame(object))
    at <- attributes(mm)
    mm <- mm[, -1]
    at$dim[2] <- at$dim[2] - 1
    at$dimnames[[2]] <- at$dimnames[[2]][-1]
    at$assign <- at$assign[-1]
    attributes(mm) <- at
    mm
}

# free1way summary

.print.free1way <- function(x, test = c("Permutation", "Wald", "LRT", "Rao"), 
                           alternative = c("two.sided", "less", "greater"), 
                           tol = .Machine$double.eps, ...)
{

    test <- match.arg(test)
    alternative <- match.arg(alternative)

    ### global
    cf <- coef(x)
    if ((length(cf) > 1L || test == "LRT") && alternative != "two.sided") 
        stop("Cannot compute one-sided p-values")

    DF <- NULL
    parm <- seq_along(cf)
    value <- 0

    # statistics
    
    if (test == "Wald") {
        # Wald statistic
        
        if (alternative == "two.sided") {
            STATISTIC <- c("Wald chi-squared" = c(crossprod(cf, x$hessian %*% cf)))
            DF <- c("df" = length(parm))
            PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
        } else {
            STATISTIC <- c("Wald Z" = c(cf * sqrt(c(x$hessian))))
            PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
        }
        
    } else if (test == "LRT") {
        # LRT
        
        par <- x$par
        par[parm] <- value
        unll <- x$value ### neg logLik
        rnll <- x$profile(par, parm)$value ### neg logLik
        STATISTIC <- c("logLR chi-squared" = - 2 * (unll - rnll))
        DF <- c("df" = length(parm))
        PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
        
    } else if (test == "Rao") {
        # Rao
        
        par <- x$par
        par[parm] <- value
        ret <- x$profile(par, parm)
        if (alternative == "two.sided") {
            STATISTIC <- c("Rao chi-squared" = c(crossprod(ret$negscore, ret$vcov %*% ret$negscore)))
            DF <- c("df" = length(parm))
            PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
        } else {
            STATISTIC <- c("Rao Z" = -ret$negscore * sqrt(c(ret$vcov)))
            PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
        }
        
    } else if (test == "Permutation") {
        # Permutation
        
        par <- x$par
        par[parm] <- value
        ret <- x$profile(par, parm)
        sc <- -ret$negscore
        if (length(cf) == 1L)
           sc <- sc / sqrt(c(ret$hessian))
        Esc <- sc - x$perm$Expectation
        if (alternative == "two.sided" && length(cf) > 1L) {
            STATISTIC <- c("Perm chi-squared" = sum(Esc %*% solve(x$perm$Covariance) * Esc))
            ps <- x$perm$permStat
            if (!is.null(x$perm$permStat))
                PVAL <- mean(ps > STATISTIC + tol)
            else {
                DF <- c("df" = x$perm$DF)
                PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
            }
        } else {
            STATISTIC <- c("Perm Z" = Esc / sqrt(c(x$perm$Covariance)))
            if (!is.null(x$perm$permStat)) {
                if (alternative == "two.sided")
                    PVAL <- mean(abs(x$perm$permStat) > abs(STATISTIC) + tol)
                else if (alternative == "less")
                    PVAL <- mean(x$perm$permStat < STATISTIC - tol)
                else
                    PVAL <- mean(x$perm$permStat > STATISTIC + tol)
            } else {
                if (alternative == "two.sided")
                    PVAL <- pchisq(STATISTIC^2, df = 1, lower.tail = FALSE)
                else
                    PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
            }
        }
        
    }
    

    RVAL <- list(statistic = STATISTIC, parameter = DF, p.value = PVAL, 
        null.value = x$mu, alternative = alternative, method = x$method, 
        data.name = x$data.name)
    class(RVAL) <- "htest"
    return(RVAL)
}

print.free1way <- function(x, ...) {
    print(ret <- .print.free1way(x))
    return(invisible(x))
}

summary.free1way <- function(object, test, alternative = c("two.sided", "less", "greater"), 
                             tol = .Machine$double.eps, ...)
{

    if (!missing(test))
        return(.print.free1way(object, test = test, alternative = alternative, tol = tol))
   
    alternative <- match.arg(alternative)

    ESTIMATE <- coef(object)
    SE <- sqrt(diag(vcov(object)))
    STATISTIC <- unname(ESTIMATE / SE)
    if (alternative == "less") {
        PVAL <- pnorm(STATISTIC)
    } else if (alternative == "greater") {
        PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
    } else {
        PVAL <- 2 * pnorm(-abs(STATISTIC))
    }
    cfmat <- cbind(ESTIMATE, SE, STATISTIC, PVAL)
    colnames(cfmat) <- c(object$link$parm, "Std. Error", "z value",
                         switch(alternative, "two.sided" = "P(>|z|)",
                                             "less" = "P(<z)",
                                             "greater" = "P(>z)"))
    ret <- list(call = object$call, coefficients = cfmat)
    class(ret) <- "summary.free1way"
    return(ret)
}
print.summary.free1way <- function(x, ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Coefficients:\n")
    printCoefmat(x$coefficients)
}

# free1way confint

confint.free1way <- function(object, parm,
    level = .95, test = c("Permutation", "Wald", "LRT", "Rao"), 
    what = c("shift", "PI", "AUC", "OVL"), ...)
{

    test <- match.arg(test)
    conf.level <- 1 - (1 - level) / 2

    cf <- coef(object)
    if (missing(parm)) 
        parm <- seq_along(cf)

    CINT <- confint.default(object, level = level)
    if (test == "Wald")
        return(CINT)
    wlevel <- level
    wlevel <- 1 - (1 - level) / 10
    CINT[] <- confint.default(object, level = wlevel)

    sfun <- function(value, parm, quantile) {
        x <- object
        alternative <- "two.sided"
        tol <- .Machine$double.eps
        # statistics
        
        if (test == "Wald") {
            # Wald statistic
            
            if (alternative == "two.sided") {
                STATISTIC <- c("Wald chi-squared" = c(crossprod(cf, x$hessian %*% cf)))
                DF <- c("df" = length(parm))
                PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
            } else {
                STATISTIC <- c("Wald Z" = c(cf * sqrt(c(x$hessian))))
                PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
            }
            
        } else if (test == "LRT") {
            # LRT
            
            par <- x$par
            par[parm] <- value
            unll <- x$value ### neg logLik
            rnll <- x$profile(par, parm)$value ### neg logLik
            STATISTIC <- c("logLR chi-squared" = - 2 * (unll - rnll))
            DF <- c("df" = length(parm))
            PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
            
        } else if (test == "Rao") {
            # Rao
            
            par <- x$par
            par[parm] <- value
            ret <- x$profile(par, parm)
            if (alternative == "two.sided") {
                STATISTIC <- c("Rao chi-squared" = c(crossprod(ret$negscore, ret$vcov %*% ret$negscore)))
                DF <- c("df" = length(parm))
                PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
            } else {
                STATISTIC <- c("Rao Z" = -ret$negscore * sqrt(c(ret$vcov)))
                PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
            }
            
        } else if (test == "Permutation") {
            # Permutation
            
            par <- x$par
            par[parm] <- value
            ret <- x$profile(par, parm)
            sc <- -ret$negscore
            if (length(cf) == 1L)
               sc <- sc / sqrt(c(ret$hessian))
            Esc <- sc - x$perm$Expectation
            if (alternative == "two.sided" && length(cf) > 1L) {
                STATISTIC <- c("Perm chi-squared" = sum(Esc %*% solve(x$perm$Covariance) * Esc))
                ps <- x$perm$permStat
                if (!is.null(x$perm$permStat))
                    PVAL <- mean(ps > STATISTIC + tol)
                else {
                    DF <- c("df" = x$perm$DF)
                    PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
                }
            } else {
                STATISTIC <- c("Perm Z" = Esc / sqrt(c(x$perm$Covariance)))
                if (!is.null(x$perm$permStat)) {
                    if (alternative == "two.sided")
                        PVAL <- mean(abs(x$perm$permStat) > abs(STATISTIC) + tol)
                    else if (alternative == "less")
                        PVAL <- mean(x$perm$permStat < STATISTIC - tol)
                    else
                        PVAL <- mean(x$perm$permStat > STATISTIC + tol)
                } else {
                    if (alternative == "two.sided")
                        PVAL <- pchisq(STATISTIC^2, df = 1, lower.tail = FALSE)
                    else
                        PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
                }
            }
            
        }
        
        return(STATISTIC - quantile)
    }

    if (test == "Permutation") {
        stopifnot(length(cf) == 1L)
        if (is.null(object$perm$permStat)) {
            qu <- qnorm(conf.level) * c(-1, 1)
        } else {
            qu <- quantile(object$perm$permStat, probs = c(1 - conf.level, conf.level))
            att.level <- mean(object$perm$permStat > qu[1] & object$perm$permStat < qu[2])
            attr(CINT, "Attained level") <- att.level
        }
    } else {
        qu <- rep.int(qchisq(level, df = 1), 2) ### always two.sided
    }

    for (p in parm) {
        CINT[p, 1] <- uniroot(sfun, interval = c(CINT[p,1], cf[p]), parm = p, quantile = qu[2])$root
        CINT[p, 2] <- uniroot(sfun, interval = c(cf[p], CINT[p, 2]), parm = p, quantile = qu[1])$root
    }

    what <- match.arg(what)
    CINT <- switch(what, "shift" = CINT,
                         "PI" = object$link$parm2PI(CINT),
                         "AUC" = object$link$parm2PI(CINT), ### same as PI 
                         "OVL" = object$link$parm2OVL(CINT))
    return(CINT)
}

# free1way formula

free1way.test.formula <- function(formula, data, weights, subset, na.action = na.pass, ...)
{

    cl <- match.call()

    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")

    strata <- function(object) object
    formula <- terms(formula, specials = "strata")

    stratum <- attr(formula, "specials")$strata
    if (is.null(stratum)) stratum <- 0L
    
    if (length(attr(formula, "term.labels")) > 1L + stratum)
        stop("'formula' missing or incorrect")
    group <- attr(formula, "term.labels") 
    if (stratum) group <- group[-(stratum - 1L)]

    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    DNAME <- paste(vn <- c(names(mf)[response], group), collapse = " by ") # works in all cases
    w <- as.vector(model.weights(mf))
    y <- mf[[response]]
    g <- mf[[group]]
    stopifnot(is.factor(g))
    lev <- levels(g)
    DNAME <- paste(DNAME, paste0("(", paste0(lev, collapse = ", "), ")"))
    if (nlevels(g) < 2L)
        stop("grouping factor must have at least 2 levels")
    if (stratum) {
        st <- factor(mf[[stratum]], levels = )
        if (nlevels(st) < 2L)
            stop("at least two strata must be present")
        vn <- c(vn, names(mf)[stratum])
        RVAL <- free1way.test(y = y, x = g, z = st, weights = w, 
                              varnames = vn, ...)
        DNAME <- paste(DNAME, paste("\n\t stratified by", names(mf)[stratum]))
    } else {
        ## Call the default method.
        RVAL <- free1way.test(y = y, x = g, weights = w, varnames = vn, ...)
    }
    RVAL$data.name <- DNAME
    RVAL$call <- cl
    RVAL
}

# free1way numeric

free1way.test.numeric <- function(y, x, z = NULL, event = NULL, weights = NULL, nbins = 0, 
    varnames = c(deparse1(substitute(y)), 
                 deparse1(substitute(x)), 
                 deparse1(substitute(z))), ...) {

    cl <- match.call()
    DNAME <- paste(varnames[1], "by", varnames[2])
    DNAME <- paste(DNAME, paste0("(", paste0(levels(x), collapse = ", "), ")"))

    if (!is.null(z))
        DNAME <- paste(DNAME, "\n\t stratified by", varnames[3])
    varnames <- varnames[varnames != "NULL"]

    if (!is.null(event)) {
        stopifnot(is.logical(event))
        uy <- sort(unique(y[event]))
        if (all(y[!event] < uy[length(uy)]))
            uy <- uy[-length(uy)]
    } else {
        uy <- sort(unique(y))
    }
    if (nbins && nbins < length(uy)) {
        nbins <- ceiling(nbins)
        breaks <- c(-Inf, quantile(y, prob = seq_len(nbins) / (nbins + 1L)), Inf)
    } else {
        breaks <- c(-Inf, uy, Inf)
    }
    r <- cut(y, breaks = breaks, ordered_result = TRUE)[, drop = TRUE]
    RVAL <- free1way.test(y = r, x = x, z = z, event = event, weights = weights, 
                          varnames = varnames, ...)
    RVAL$data.name <- DNAME
    RVAL$call <- cl
    RVAL
}

# free1way factor

free1way.test.factor <- function(y, x, z = NULL, event = NULL, weights = NULL, 
    varnames = c(deparse1(substitute(y)), 
                 deparse1(substitute(x)), 
                 deparse1(substitute(z))), ...) {

    cl <- match.call()
    DNAME <- paste(varnames[1], "by", varnames[2])
    DNAME <- paste(DNAME, paste0("(", paste0(levels(x), collapse = ", "), ")"))

    if (!is.null(z))
        DNAME <- paste(DNAME, "\n\t stratified by", varnames[3])
    varnames <- varnames[varnames != "NULL"]

    stopifnot(is.factor(x))
    if (nlevels(y) > 2L)
        stopifnot(is.ordered(y))
    d <- data.frame(w = 1, y = y, x = x)
    if (!is.null(weights)) d$w <- weights
    if (is.null(z)) z <- gl(1, nrow(d))
    d$z <- z 
    if (!is.null(event)) {
        stopifnot(is.logical(event))
        d$event <- event
    }
    tab <- xtabs(w ~ ., data = d)
    dn <- dimnames(tab)
    names(dn)[seq_along(varnames)] <- varnames
    dimnames(tab) <- dn
    RVAL <- free1way.test(tab, ...)
    RVAL$data.name <- DNAME
    RVAL$call <- cl
    RVAL
}

# ppplot

ppplot <- function(x, y, plot.it = TRUE,
            xlab = deparse1(substitute(x)),
            ylab = deparse1(substitute(y)), 
            interpolate = FALSE, ...,
            conf.level = NULL, conf.args = list(type = "Wald", col = NA, border = NULL)) {

    force(xlab)
    force(ylab)
    if (xlab == ylab) {
        xlab <- paste0("x = ", xlab)
        ylab <- paste0("y = ", ylab)
    }

    ex <- ecdf(x)
    if (interpolate) {
        vals <- sort(unique(x))
        ex <- splinefun(vals, ex(vals), method = "hyman")
    }
    sy <- sort(unique(y))
    py <- ecdf(y)(sy)
    px <- ex(sy)
    ret <- list(x = px, y = py)
    if (!plot.it)
        return(ret)

    plot(px, py, xlim = c(0, 1), ylim = c(0, 1), 
         xlab = xlab, ylab = ylab, type = "n", ...)

    # ROC bands
    
     if (!is.null(conf.level)) {
        prb <- seq_len(1000) / 1001
        res <- c(x, y)
        grp <- gl(2, 1, labels = c(xlab, ylab))
        grp <- grp[rep(1:2, c(length(x), length(y)))]
        args <- conf.args
        args$y <- res
        args$x <- grp
        args$border <- args$col <- args$type <- NULL
        f1w <- do.call("free1way.test", args)

        ci <- confint(f1w, level = conf.level, type = args$type)
        lwr <- .p(f1w$link, .q(f1w$link, prb) - ci[1,1])
        upr <- .p(f1w$link, .q(f1w$link, prb) - ci[1,2])
        x <- c(prb, rev(prb))
        y <- c(lwr, rev(upr))
        xn <- c(x[1L], rep(x[-1L], each = 2))
        yn <- c(rep(y[-length(y)], each = 2), y[length(y)])
        polygon(x = xn, y = yn, col = conf.args$col, border = conf.args$border)
        lines(prb, .p(f1w$link, .q(f1w$link, prb) - coef(f1w)))
    }
    

    points(px, py, ...)
    return(invisible(ret)) 
}

# r2dsim

r2dsim <- function(n, r, c, delta = 0,
                   link = c("logit", "probit", "cloglog", "loglog")) 
{

    if (length(n <- as.integer(n)) == 0L || (n < 0) || is.na(n)) 
        stop("invalid argument 'n'")
    colsums <- c
    if (length(colsums[] <- as.integer(c)) <= 1L || 
        any(colsums < 0) || anyNA(colsums)) 
        stop("invalid argument 'c'")

    prob <- r
    if (length(prob[] <- as.double(r / sum(r))) <= 1L || 
        any(prob < 0) || anyNA(prob)) 
        stop("invalid argument 'r'")

    if (is.null(names(prob))) 
        names(prob) <- paste0("i", seq_along(prob))
    
    K <- length(colsums)
    if (is.null(names(colsums))) 
        names(colsums) <- LETTERS[seq_len(K)]
    delta <- rep_len(delta, K - 1L)

    # link2fun
    
    if (!inherits(link, "linkfun")) {
        link <- match.arg(link)
        link <- do.call(link, list())
    }
    

    p0 <- cumsum(prob)
    h0 <- .q(link, p0)

    h1 <- h0 - matrix(delta, nrow = length(prob), ncol = K - 1, byrow = TRUE)
    p1 <- .p(link, h1)
    p <- cbind(p0, p1)
    ret <- vector(mode = "list", length = n)

    for (i in seq_len(n)) {
        tab <- sapply(seq_len(K), function(k) 
            rmultinom(1L, size = colsums[k], 
                      prob = c(p[1,k], diff(p[,k]))))
        ret[[i]] <- as.table(array(tab, dim = c(length(prob), K), 
                          dimnames = list(names(prob), 
                                          names(colsums))))
    }
    return(ret)
}

# power

power.free1way.test <- function(n = NULL, prob = rep.int(1 / n, n), 
                                alloc_ratio = 1, strata_ratio = 1, 
                                delta = NULL, mu = 0, sig.level = .05, power = NULL,
                                link = c("logit", "probit", "cloglog", "loglog"),
                                alternative = c("two.sided", "less", "greater"), 
                                nsim = 100, seed = NULL, tol = .Machine$double.eps^0.25) 
{

    if (sum(vapply(list(n, delta, power, sig.level), is.null, 
        NA)) != 1) 
        stop("exactly one of 'n', 'delta', 'power', and 'sig.level' must be NULL")
    stats:::assert_NULL_or_prob(sig.level)
    stats:::assert_NULL_or_prob(power)

    # random seed
    
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        seed <- RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    

    if (is.null(n)) 
        n <- ceiling(uniroot(function(n) {
                 # power call
                 
                 power.free1way.test(n = n, prob = prob, alloc_ratio = alloc_ratio, 
                                     strata_ratio = strata_ratio, delta = delta, mu = mu,
                                     sig.level = sig.level, link = link, 
                                     alternative = alternative, 
                                     nsim = nsim, seed = seed, tol = tol)$power - power
                 
             }, interval = c(5, 1e+03), tol = tol, extendInt = "upX")$root)
    else if (is.null(delta)) {
        ### 2-sample only
        stopifnot(K == 2L)
        delta <- uniroot(function(delta) {
                 # power call
                 
                 power.free1way.test(n = n, prob = prob, alloc_ratio = alloc_ratio, 
                                     strata_ratio = strata_ratio, delta = delta, mu = mu,
                                     sig.level = sig.level, link = link, 
                                     alternative = alternative, 
                                     nsim = nsim, seed = seed, tol = tol)$power - power
                 
    ### <TH> interval depending on alternative, symmetry? </TH>
            }, interval = c(0, 10), tol = tol, extendInt = "upX")$root
        }
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) {
                # power call
                
                power.free1way.test(n = n, prob = prob, alloc_ratio = alloc_ratio, 
                                    strata_ratio = strata_ratio, delta = delta, mu = mu,
                                    sig.level = sig.level, link = link, 
                                    alternative = alternative, 
                                    nsim = nsim, seed = seed, tol = tol)$power - power
                
            }, interval = c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
    
    # power setup
    

    # link2fun

    if (!inherits(link, "linkfun")) {
        link <- match.arg(link)
        link <- do.call(link, list())
    }
    

    ### matrix means control distributions in different strata
    if (!is.matrix(prob))
        prob <- matrix(prob, nrow = NROW(prob))
    prob <- prop.table(prob, margin = 2L)
    C <- nrow(prob)
    K <- length(delta) + 1L
    B <- ncol(prob)
    if (is.null(colnames(prob))) 
        colnames(prob) <- paste0("stratum", seq_len(B))
    if (is.null(names(delta))) 
        names(delta) <- LETTERS[seq_len(K)[-1]]
    p0 <- apply(prob, 2, cumsum)
    h0 <- .q(link, p0)
    if (length(alloc_ratio) == 1L) 
        alloc_ratio <- rep_len(alloc_ratio, K - 1)
    stopifnot(length(alloc_ratio) == K - 1)
    if (length(strata_ratio) == 1L) 
        strata_ratio <- rep_len(strata_ratio, B - 1)
    stopifnot(length(strata_ratio) == B - 1)
    ### sample size per group (columns) and stratum (rows)
    N <- n * matrix(c(1, alloc_ratio), nrow = B, ncol = K, byrow = TRUE) * 
             matrix(c(1, strata_ratio), nrow = B, ncol = K)
    rownames(N) <- colnames(prob)
    ctrl <- "Control"
    dn <- dimnames(prob)
    if (!is.null(names(dn)[1L]))
        ctrl <- names(dn)[1L]
    colnames(N) <- c(ctrl, names(delta))
    
    # estimate Fisher information
    
    he <- 0
    deltamu <- delta - mu
    for (i in seq_len(nsim)) {
        parm <- deltamu
        x <- as.table(array(0, dim = c(C, K, B)))
        for (b in seq_len(B)) {
            x[,,b] <- r2dsim(1L, r = prob[, b], c = N[b,], delta = delta, link = link)[[1L]]
            rs <- rowSums(x[,,b]) > 0
            h <- h0[rs, b]
            theta <- c(h[1], diff(h[-length(h)]))
            parm <- c(parm, theta)
        }
        ### evaluate observed hessian for true parameters parm and x data
        he <- he + .free1wayML(x, link = link, mu = mu, start = parm, fix = seq_along(parm))$hessian
    }
    ### estimate expected Fisher information
    he <- he / nsim
    


    alternative <- match.arg(alternative)
    if (K == 2L) {
        se <- 1 / sqrt(c(he))
        power  <- switch(alternative, 
            "two.sided" = pnorm(qnorm(sig.level / 2) + deltamu / se) + 
                          pnorm(qnorm(sig.level / 2) - deltamu / se),
            "less" = pnorm(qnorm(sig.level) - deltamu / se),
            "greater" = pnorm(qnorm(sig.level) + deltamu / se)
        )
    } else {
        stopifnot(alternative == "two.sided")
        ncp <- sum((chol(he) %*% deltamu)^2)
        qsig <- qchisq(sig.level, df = K - 1L, lower.tail = FALSE)
        power <- pchisq(qsig, df = K - 1L, ncp = ncp, lower.tail = FALSE)
    }

    # power htest output
    
    ss <- paste(colSums(N), paste0("(", colnames(N), ")"), collapse = " + ")
    ret <- list(n = n, 
                "Total sample size" = paste(ss, "=", sum(N)),
                power = power, 
                sig.level = sig.level)
    if (mu != 0) ret$mu <- mu
    ret[[link$parm]] <- delta
    ret$note <- "'n' is sample size in control group"
    if (B > 1) ret$note <- paste(ret$note, "of first stratum")
    alias <- link$alias
    if (length(link$alias) == 2L) alias <- alias[1L + (K > 2L)]
    ret$method <- paste(ifelse(B > 1L, "Stratified", ""), 
                        paste0(K, "-sample"), alias, 
                        "test against", link$model, "alternatives")
    class(ret) <- "power.htest"
    

    ret
}

