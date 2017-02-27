	
interval <- function(x, ...)
    UseMethod("interval")

interval.default <- function(x, ...)
    stop("no interval method for class", " ", sQuote(class(x)), 
         " ", "found")

interval.numeric <- function(x, breaks = 50, ...) {

    ### from cut.default()
    if (length(breaks) == 1L) {
        if (is.na(breaks) || breaks < 2L) 
            stop("invalid number of intervals")
        nb <- as.integer(breaks + 1)
        dx <- diff(rx <- range(x, na.rm = TRUE))
        if (dx == 0) {
            dx <- abs(rx[1L])
            breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, 
                length.out = nb)
        }
        else {
            breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
            breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + 
                dx/1000)
        }
    } else {
        breaks <- sort(as.double(breaks))
    }
    if (anyDuplicated(breaks)) 
        stop("'breaks' are not unique")

    ret <- cut.default(x, breaks = breaks, labels = FALSE)
    ret[is.na(x)] <- 0L
    attr(ret, "levels") <- breaks
    class(ret) <- c("interval", "integer")
    ret
}


levels.interval <- function(x) {
    breaks <- attr(x, "levels")
    return(paste("(", breaks[-length(breaks)], ",", 
                 breaks[-1], "]", sep = ""))
}

nlevels.interval <- function(x)
    length(attr(x, "levels")) - 1L

print.interval <- function(x, quote = FALSE, max.levels = NULL, 
                           width = getOption("width"), ...) 
{

    print(c("<NA>", levels(x))[x + 1L], quote = quote)
    maxl <- if (is.null(max.levels)) 
        TRUE
    else max.levels
    if (maxl) {
        n <- length(lev <- encodeString(levels(x), quote = ifelse(quote, 
            "\"", "")))
        colsep <- " < "
        T0 <- "Intervals: "
        if (is.logical(maxl)) 
            maxl <- {
                width <- width - (nchar(T0, "w") + 3L + 1L + 
                  3L)
                lenl <- cumsum(nchar(lev, "w") + nchar(colsep, 
                  "w"))
                if (n <= 1L || lenl[n] <= width) 
                  n
                else max(1L, which.max(lenl > width) - 1L)
            }
        drop <- n > maxl
        cat(if (drop) 
            paste(format(n), ""), T0, paste(if (drop) 
            c(lev[1L:max(1, maxl - 1)], "...", if (maxl > 1) lev[n])
        else lev, collapse = colsep), "\n", sep = "")
    }
    return(invisible(x))
}

"[.interval" <- function(x, i, ..., drop = FALSE) {
    ix <- unclass(x)
    ret <- ix[i]
    lev <- attr(x, "levels")
    if (drop)
        stop(sQuote("drop = TRUE"), " ", "not implemented in", " ", 
             sQuote("[.interval"))
    attr(ret, "levels") <- lev
    class(ret) <- class(x)
    ret
}

format.interval <- function(x, ...)
    c("<NA>", levels(x))[x + 1L]

is.na.interval <- function(x)
    unclass(x) == 0L
