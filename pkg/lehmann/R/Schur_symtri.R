
###
### Schur complement Z - t(X) %*% solve(A) %*% X of a block matrix 
###
###      | A     X |
###      | t(X)  Z |
###
### when A is tridiagonal symmetric with diagonal a
### and off-diagonal b
### 
Schur_symtri <- function(a, b, X, Z, tol = sqrt(.Machine$double.eps)) {

    N <- length(a)
    stopifnot(length(b) == N - 1)
    P <- NROW(Z)
    stopifnot(NCOL(X) == P)
    stopifnot(NCOL(Z) == P)
    stopifnot(NROW(X) == N)
    if (!is.matrix(X)) X <- matrix(X, nrow = N)
    storage.mode(X) <- "double"

    idx <- which(abs(b) < .Machine$double.eps)
    ### is problem reducible?
    if (length(idx) > 0) {
        qf <- 0
        idx <- c(0, idx, N)
        for (i in 1:(length(idx) - 1)) {
            ix <- (idx[i] + 1L):idx[i + 1L]
            qf <- qf + .Call("R_symtrisolve_quadform", 
                             as.double(a)[ix], 
                             as.double(b)[ix[-length(ix)]], 
                             X[ix,,drop = FALSE], 
                             as.double(tol))
        }
    } else {
        qf <- .Call("R_symtrisolve_quadform", 
                    as.double(a), 
                    as.double(b), X, 
                    as.double(tol))
    }

    return(Z - qf)
}
