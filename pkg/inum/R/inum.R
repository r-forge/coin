
inum <- function(object, nmax = 20, ...)
    UseMethod("inum")

inum.default <- function(object, nmax = 20, ...)
    stop("cannot handle objects of class", " ", sQuote(class(object)))

inum.data.frame <- function(object, nmax = 20, ignore = NULL, total = FALSE, 
                           weights = NULL, as.interval = "", 
                           complete.cases.only = FALSE, ...) {

    if (total) {
        bdr <- inum(object, nmax = nmax, ignore = ignore, 
                   total = FALSE, as.interval = as.interval) 
        bdr2 <- lapply(bdr, function(x)
            factor(x, levels = 0:nlevels(x)))
        ret <- do.call("interaction", bdr2)
        if (!is.null(weights)) {
            tab <- xtabs(weights ~ ret)
        } else {
            tab <- table(ret)
        }
        tab0 <- which(tab > 0)

        sDF <- vector(mode = "list", length = length(bdr))
        len <- sapply(bdr2, nlevels)
        ### do.call("expand.grid", bdr), essentially
        for (j in 1:length(len)) {
            ix <- 1:len[j]
            if (j > 1)
                ix <- rep(ix, each = prod(len[1:(j - 1)]))
            idx <- rep(ix, length.out = prod(len))[tab0]
            if (inherits(bdr[[j]], "interval")) {
                sDF[[j]] <- (0:nlevels(bdr[[j]]))[idx]
                attr(sDF[[j]], "levels") <- attr(bdr[[j]], "levels")
                class(sDF[[j]]) <- class(bdr[[j]])
            } else {
                lev <- attr(bdr[[j]], "levels")
                lev <- lev[c(1, 1:length(lev))]
                lev[1] <- NA
                sDF[[j]] <- lev[idx, drop = FALSE]
            }
        }
        ### note: sDF contains missings and 
        ### ret is always > 0 (is, no missings)
        ### this is different from enum/integer types
        ### should we handle this here?
        sDF <- as.data.frame(sDF)
        colnames(sDF) <- names(bdr)
        sDF[["(weights)"]] <- as.numeric(tab[tab0])
        rownames(sDF) <- NULL
        ret <- unclass(ret[, drop = TRUE])

        if (complete.cases.only) {
            cc <- rowSums(sapply(sDF[colnames(sDF) != "(weights)"], 
                                 function(x) is.na(x))) == 0
            cc[is.na(cc)] <- TRUE
            if (any(!cc)) {
                sDF <- sDF[cc,,drop = FALSE]
                ret[!cc] <- 0L
                ret <- unclass(factor(ret)) - 1L
            }
        }  

        attr(ret, "levels") <- sDF
        class(ret) <- "inumtotal"
        return(ret)
    }

    ret <- vector(mode = "list", length = ncol(object))
    names(ret) <- cn <- colnames(object)

    if (!is.null(ignore)) {
        if (is.integer(ignore)) cn <- cn[-ignore]
        if (is.character(ignore)) cn <- cn[cn != ignore]
    }

    if (as.interval != "") {
        if (!is.character(as.interval))
            stop(sQuote("as.interval"), " ", "is not a character")
    }

    for (v in cn) {
        x <- object[[v]]
        if (is.logical(x) || is.factor(x) || is.integer(x)) {
            ix <- enum(x)
        } else if (is.numeric(x)) {
            ux <- sort(unique(x))
            xmin <- ux[1]
            xmax <- ux[length(ux)] 
            if (length(ux) > nmax)
                ux <- unique(quantile(x, prob = 1:(nmax - 1L) / nmax, 
                                      na.rm = TRUE))
            ux <- ux[ux < xmax]
            tol <- min(diff(sort(ux))) ### sqrt(.Machine$double.eps)
            ix <- interval(x, breaks = c(xmin - tol, ux, xmax))
            if (all(as.interval != v)) {
                ### <FIXME> this minimises distances to original
                ### measurements but leads to incorrect cutpoints
                ### (where c(ux, xmax) would be OK)
                # nux <- c(xmin, ux) + diff(c(xmin, ux, xmax)) / 2
                nux <- c(ux, xmax)
                ### </FIXME>
                attr(ix, "levels") <- as.double(nux)
                class(ix) <- c("enum", "integer")
             }
        } else if (is.data.frame(x)) {
            ix <- inum(x, nmax = nmax, ignore = ignore, total = TRUE,
                      as.interval = as.interval)
        } else {
            ix <- inum(x, nmax = nmax, ...) ### nothing as of now
        }
        ret[[v]] <- ix
    }
    class(ret) <- "inum"
    ret
}

### only useful for checks
as.data.frame.inum <- function(x, ...) {
    ret <- lapply(x, function(x) {
        if (inherits(x, "interval")) return(x)
        lev <- attr(x, "levels")
        lev <- lev[c(1, 1:length(lev))]
        lev[1] <- NA
        return(lev[x + 1])
    })
    class(ret) <- "data.frame"
    attr(ret, "row.names") <- 1:NROW(ret[[1]])
    ret
}

as.data.frame.inumtotal <- function(x, ...) 
    attr(x, "levels")

weights.inumtotal <- function(object, ...)
    attr(object, "levels")[["(weights)"]]

### does not make sense
# is.numeric.Surv <- function(x, ...)
#    return(FALSE)
# inum.Surv <- function(object, nmax = 20, ...) {
#     x <- inum(as.data.frame(unclass(object)), nmax = nmax, total = TRUE)
#     lev <- as.matrix(attr(x, "levels"))
#     atr <- attributes(object)
#     atr$dim <- dim(lev)
#     atr$dimnames <- dimnames(lev)
#     attributes(lev) <- atr
#     attr(x, "levels") <- lev
#     x
# }
