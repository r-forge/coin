
Schur_symtri <- function(a, b, X, Z) {

    n <- length(a)
    stopifnot(length(b) == n - 1)
    p <- NROW(Z)
    stopifnot(NCOL(X) == p)
    stopifnot(NCOL(Z) == p)
    stopifnot(NROW(X) == n)
    if (!is.matrix(X)) X <- matrix(X, nrow = n)
    storage.mode(X) <- "double"

    idx <- which(abs(b) < .Machine$double.eps)
    if (length(idx) > 0) {
        qf <- 0
        idx <- c(0, idx, n)
        for (i in 1:(length(idx) - 1)) {
            ix <- (idx[i] + 1):idx[i + 1]
            qf <- qf + .Call("R_symtrisolve_quadform", as.double(a)[ix], 
                                                       as.double(b)[ix[-length(ix)]], 
                                                       X[ix,,drop = FALSE])
        }
    } else {
        qf <- .Call("R_symtrisolve_quadform", as.double(a), as.double(b), X)
    }

    solve(Z - qf)
}
