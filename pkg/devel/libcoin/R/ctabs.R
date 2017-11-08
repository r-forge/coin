
# R Header

###
### TO NOT EDIT THIS FILE
### 
### Edit `libcoin.w' and run `nuweb -r libcoin.w'
###

ctabs <- function(ix, iy = integer(0), block = integer(0), weights = integer(0),
                    subset = integer(0), checkNAs = TRUE)
{

    stopifnot(is.integer(ix) || is.factor(ix))
    N <- length(ix)

    # Check ix
    
    if (is.null(attr(ix, "levels"))) {
        rg <- range(ix)
        if (any(is.na(rg)))
            stop("no missing values allowed in ix") 
        stopifnot(rg[1] >= 0)
        attr(ix, "levels") <- 1:rg[2]
    } else {
        if (checkNAs) stopifnot(all(!is.na(ix)))
    }
    

    if (length(iy) > 0) {
        stopifnot(length(iy) == N)
        stopifnot(is.integer(iy) || is.factor(iy))
        # Check iy
        
        if (is.null(attr(iy, "levels"))) {
            rg <- range(iy)
            if (any(is.na(rg)))
                stop("no missing values allowed in iy") 
            stopifnot(rg[1] >= 0)
            attr(iy, "levels") <- 1:rg[2]
        } else {
            if (checkNAs) stopifnot(all(!is.na(ix)))
        }
        
    }

    # Check weights, subset, block
    

    if (is.null(weights)) weights <- integer(0)

    if (length(weights) > 0) {
        if (!((N == length(weights)) && all(weights >= 0)))
            stop("incorrect weights")
        if (checkNAs) stopifnot(all(!is.na(weights)))
    }

    if (is.null(subset)) subset <- integer(0)

    if (length(subset) > 0) {
        rs <- range(subset)
        if (any(is.na(rs))) stop("no missing values allowed in subset")
        if (!((rs[2] <= N) && (rs[1] >= 1L)))
            stop("incorrect subset")
    }

    if (is.null(block)) block <- integer(0)

    if (length(block) > 0) {
        if (!((N == length(block)) && is.factor(block)))
            stop("incorrect block")
        if (checkNAs) stopifnot(all(!is.na(block)))
    }
    

    if (length(iy) == 0 && length(block) == 0)
        return(.Call(R_OneTableSums, ix, weights, subset))
    if (length(block) == 0)
        return(.Call(R_TwoTableSums, ix, iy, weights, subset))
    if (length(iy) == 0)
        return(.Call(R_TwoTableSums, ix, block, weights, subset)[,-1,drop = FALSE])
    return(.Call(R_ThreeTableSums, ix, iy, block, weights, subset))
}
