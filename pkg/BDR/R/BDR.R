
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
    atr$dimnames <- dimnames(lev)
    attributes(lev) <- atr
    attr(x, "levels") <- lev
    x
}

BDR.data.frame <- function(object, nmax = 20, ignore = NULL, total = FALSE, weights = NULL, ...) {

    if (total) {
        bdr <- BDR(object, nmax = nmax, ignore = ignore, total = FALSE) 
        ret <- do.call("interaction", bdr)
        if (!is.null(weights)) {
            tab <- xtabs(weights ~ ret)
        } else {
            tab <- table(ret)
        }
        tab0 <- which(tab > 0)

        sDF <- vector(mode = "list", length = length(bdr))
        len <- sapply(bdr, function(x) NROW(attr(x, "levels")))
        ### do.call("expand.grid", bdr), essentially
        for (j in 1:length(len)) {
            ix <- 1:len[j]
            if (j > 1)
                ix <- rep(ix, each = prod(len[1:(j - 1)]))
            idx <- rep(ix, length.out = prod(len))[tab0]
            lev <- attr(bdr[[j]], "levels")
            idx[idx < 1] <- NA
            if (is.null(dim(lev)))
                sDF[[j]] <- lev[idx, drop = FALSE]
            else
                sDF[[j]] <- lev[idx,,drop = FALSE]
        }
        sDF <- as.data.frame(sDF)
        sDF[["(weights)"]] <- tab[tab0]
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
            ### <FIXME> use findInterval instead of cut ?
            if (length(ux) <= nmax) {
                ix <- cut.default(x, breaks = c(-Inf, ux, Inf), 
                                  right = TRUE, labels = FALSE)
            } else {
                ux <- unique(quantile(x, prob = 1:nmax / nmax, na.rm = TRUE))
                ix <- cut.default(x, breaks = c(-Inf, ux, Inf), 
                                  right = TRUE, labels = FALSE)
            }
            ### </FIXME>
            levels(ix) <- ux
        } else if (is.factor(x) && !is.ordered(x)) {
            ix <- unclass(x)
            attr(ix, "levels") <- factor(levels(x), levels = levels(x), 
                                         labels = levels(x))
        } else if (is.ordered(x)) {
            ix <- unclass(x)
            attr(ix, "levels") <- ordered(levels(x), 
                                          levels = levels(x), 
                                          labels = levels(x))
        } else if (is.data.frame(x)) {
            ix <- BDR(x, nmax = nmax, ignore = ignore, total = TRUE)
        } else {
            ix <- BDR(x, nmax = nmax, ...)
        }
        ix[is.na(ix)] <- 0L
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

weights.BDRtotal <- function(object, ...)
    attr(object, "levels")[["(weights)"]]

