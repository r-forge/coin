
is.numeric.Surv <- function(x, ...)
    return(FALSE)

BDR <- function(object, nmax = 20, ...)
    UseMethod("BDR")

BDR.default <- function(object, nmax = 20, ...)
    stop("cannot handle objects of class", " ", sQuote(class(object)))

BDR.Surv <- function(object, nmax = 20, ...) {
    x <- BDR(as.data.frame(unclass(object)), nmax = nmax, total = TRUE)
    lev <- as.matrix(attr(x, "levels"))
    atr <- attributes(object)
    atr$dim <- dim(lev)
    attributes(lev) <- atr
    attr(x, "levels") <- lev
    x
}

BDR.data.frame <- function(object, nmax = 20, ignore = NULL, total = FALSE, ...) {

    if (total) {
        bdr <- BDR(object, nmax = nmax, ignore = ignore, total = FALSE) 
        ret <- do.call("interaction", bdr)
        lDF <- do.call("expand.grid", 
            lapply(bdr, function(x) {
                ret <- attr(x, "levels")
                if (!is.null(dim(ret)))
                    ret <- 1:NROW(ret)
                if (any(x < 1)) {
                    ret <- ret[c(1, 1:length(ret))]
                    ret[1] <- NA
                }
                ret
            }))
        for (j in 1:NCOL(object)) {
            if (!is.null(dim(object[[j]])))
                lDF[[j]] <- attr(bdr[[j]], "levels")[lDF[[j]],,drop = FALSE]
        }
        sDF <- lDF[table(ret) > 0,,drop = FALSE]
        rownames(sDF) <- NULL
        ret <- unclass(ret[, drop = TRUE])
        attr(ret, "levels") <- sDF
        class(ret) <- "BDRtotal"
        return(ret)
    }

    ret <- vector(mode = "list", length = ncol(object))
    names(ret) <- cn <- colnames(object)

    if (!is.null(ignore)) {
        if (is.integer(ignore)) cn <- cn[-ignore]
        if (is.character(ignore)) cn <- cn[cn != ignore]
    }

    for (v in cn) {
        x <- object[[v]]
        if (is.logical(x)) {
            ix <- x + 1L
            levels(ix) <- ux <- c(FALSE, TRUE)
            X <- rbind(0, diag(2))
        } else if (is.numeric(x)) {
            ux <- sort(unique(x))
            if (length(ux) <= nmax) {
                ix <- cut.default(x, breaks = c(-Inf, ux, Inf), 
                                  right = TRUE, labels = FALSE)
            } else {
                ux <- unique(quantile(x, prob = 1:nmax / nmax, na.rm = TRUE))
                ix <- cut.default(x, breaks = c(-Inf, ux, Inf), 
                                  right = TRUE, labels = FALSE)
            }
            levels(ix) <- ux
            X <- rbind(0, matrix(ux, ncol = 1L))
        } else if (is.factor(x) && !is.ordered(x)) {
            ix <- unclass(x)
            attr(ix, "levels") <- factor(levels(x), levels = levels(x), 
                                         labels = levels(x))
            X <- rbind(0, diag(nlevels(x)))
        } else if (is.ordered(x)) {
            ix <- unclass(x)
            attr(ix, "levels") <- ordered(levels(x), 
                                          levels = levels(x), 
                                          labels = levels(x))
            sc <- attr(x, "scores")
            if (is.null(sc)) sc <- 1:nlevels(x)
            X <- rbind(0, matrix(sc, ncol = 1L))
        } else if (is.data.frame(x)) {
            ix <- BDR(x, nmax = nmax, ignore = ignore, total = TRUE)
            X <- NA
        } else {
            ix <- BDR(x, nmax = nmax, ...)
            X <- NA
        }
        ix[is.na(ix)] <- 0L
        storage.mode(X) <- "double"
        attr(ix, "X") <- X
        ret[[v]] <- ix
    }
    class(ret) <- "BDR"
    ret
}

as.data.frame.BDR <- function(x, ...) {
    ret <- lapply(x, function(x) {
        lev <- attr(x, "levels")
        d <- dim(lev)
        if (is.null(d) || inherits(lev, "Surv")) {
            lev <- lev[c(1, 1:NROW(lev))]
            lev[1] <- NA
            return(lev[x + 1])
        }
        if (length(d) > 2) stop("cannot handle arrays")
        x[x < 1] <- NA
        return(lev[x, ,drop = FALSE])
    })
    class(ret) <- "data.frame"
    attr(ret, "row.names") <- 1:NROW(ret[[1]])
    ret
}

as.data.frame.BDRtotal <- function(x, ...) {
    lev <- attr(x, "levels")
    x[x < 1] <- NA
    return(lev[x, , drop = FALSE])
}

