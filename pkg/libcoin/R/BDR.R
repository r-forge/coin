
df2BDR <- function(object, max = 10, ignore = NULL) {

    ret <- vector(mode = "list", length = ncol(object))
    names(ret) <- cn <- colnames(object)
    if (!is.null(ignore)) {
        if (is.integer(ignore)) cn <- cn[-ignore]
        if (is.character(ignore)) cn <- cn[cn != ignore]
    }

    for (v in colnames(object)) {
        x <- object[[v]]
        if (!is.factor(x)) {
            q <- unique(quantile(x, prob = 1:max / (max + 1)))
            ind <- cut(x, breaks = c(-Inf, q, Inf))
            ret[[v]] <- list(ind = unclass(ind),
                             val = rbind(NA, cbind(c(-Inf, q), c(q, Inf))))
        } else {
            ret[[v]] <- list(ind = unclass(x),
                             val = factor(c(NA, levels(x)), ordered = is.ordered(x)))
        }
        ret[[v]]$ind[is.na(x)] <- 0
        attributes(ret[[v]]$ind) <- NULL
    }
    ret
}

IRIS <- iris[sample(1:nrow(iris), 1000000, replace = TRUE),]

IRIS2 <- df2BDR(IRIS, 20)

object.size(IRIS)
object.size(IRIS2)

