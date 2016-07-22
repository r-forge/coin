
BDR <- function(object, nmax = 20, ignore = NULL) {

    stopifnot(is.data.frame(object))
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
                ix <- match(x, ux)
            } else {
                ux <- unique(quantile(x, prob = 1:(nmax - 1) / nmax ))
                ix <- unclass(cut(x, breaks = c(-Inf, ux, Inf)))
            }
            levels(ix) <- ux
            X <- rbind(0, matrix(ux, ncol = 1L))
        } else if (is.factor(x) && !is.ordered(x)) {
            ix <- unclass(x)
            X <- rbind(0, diag(nlevels(x)))
        } else if (is.ordered(x)) {
            ix <- unclass(x)
            X <- rbind(0, matrix(1:nlevels(x), ncol = 1L)) ### scores?
        } else {
            stop("cannot handle class", class(x))
        }
        ix[is.na(ix)] <- 0L
        storage.mode(X) <- "double"
        attr(ix, "X") <- X
        ret[[v]] <- ix
    }
    class(ret) <- "BDR"
    ret
}
