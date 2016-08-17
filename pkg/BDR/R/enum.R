
enum <- function(x)
    UseMethod("enum")

enum.default <- function(x)
    stop("no enum method for class", " ", sQuote(class(x)),
         " ", "found")

enum.factor <- function(x) {
    ret <- unclass(x)
    attr(ret, "levels") <- factor(levels(x), levels = levels(x),
                                  labels = levels(x), ordered = is.ordered(x))
    ret[is.na(x)] <- 0L
    class(ret) <- c("enum", "integer")
    ret
}

enum.logical <- function(x) {
    ret <- x + 1L
    attr(ret, "levels") <- c(FALSE, TRUE)
    ret[is.na(x)] <- 0L
    class(ret) <- c("enum", "integer")
    ret
}

enum.integer <- function(x) {
    breaks <- sort(unique(x))
    ret <- match(x, breaks)
    ret[is.na(x)] <- 0L
    attr(ret, "levels") <- breaks
    class(ret) <- c("enum", "integer")
    ret
}

enum.numeric <- function(x)
    return(enum.integer(x))


levels.enum <- function(x)
    attr(x, "levels")

nlevels.enum <- function(x)
    length(levels(x))

print.enum <- function(x, quote = FALSE, max.levels = NULL,
                       width = getOption("width"), ...)
{

    print(c("<NA>", as.character(levels(x)))[x + 1L], quote = quote)
    maxl <- if (is.null(max.levels))
        TRUE
    else max.levels
    if (maxl) {
        n <- length(lev <- encodeString(as.character(levels(x)), 
            quote = ifelse(quote, "\"", "")))
        colsep <- if (is.ordered(levels(x))) " < " else " "
        T0 <- "Levels: "
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
