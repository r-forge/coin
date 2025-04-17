
Schur_symtri <- function(a, b, X, Z) {

    n <- length(a)
    stopifnot(length(b) == n - 1)
    p <- NROW(Z)
    stopifnot(NCOL(X) == p)
    stopifnot(NCOL(Z) == p)
    stopifnot(NROW(X) == n)

    storage.mode(X) <- "double"

    ret <- Z - .Call("R_symtrisolve_quadform", as.double(a), as.double(b), X)
    solve(ret)
}
