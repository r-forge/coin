
#include "lehmann.h"

R_xlen_t NROW
(
    SEXP x
) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(XLENGTH(x));
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[0]);
    return((R_xlen_t) INTEGER(a)[0]);
}

/* NCOL */

R_xlen_t NCOL
(
    SEXP x
) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[1]);
    return((R_xlen_t) INTEGER(a)[1]);
}


/* 

   Compute inverse of a symmetric N x N tridiagonal matrix T with 
   diagonals `a' (length N) and off-diagonals `b' (length N - 1)
   
       T^{-1}_ij = u_i v_j for i <= j
   
   as described in DOI:10.1137/0613045
   
   The function computes a real N x 2 matrix `ans' with 1st column `u'
   and second column `v'
   
*/

void C_symtrisolve (double *a, double *b, R_xlen_t n, double tol, double *ans)
{

    SEXP Rd;
    double *d, *delta, *u, *v, prodb, det;
    R_xlen_t i;

    /* output vectors */
    u = ans;
    v = ans + n + 1;
    
    /* n = N - 1 */
    PROTECT(Rd = allocVector(REALSXP, n + 1));
    d = REAL(Rd);
    
    d[n] = a[n];
    det = d[n];
    for (i = n - 1; i >= 0; i--) {
        d[i] = a[i] - pow(b[i], 2) / d[i + 1];
        /* DOI:10.1137/0613045 page 710: T = U D^-1 U^t with
           diag(U) = d (upper triangular),
           diag(D) = d (diagonal) => det(T) = prod(d)
        */
        det *= d[i];
    }

    if (fabs(det) < tol) {
        error("Matrix not invertible");
    } else {
        u[0] = 1 / d[0];
        prodb = 1.0;
        for (i = 1; i <= n; i++) {
            prodb *= -b[i - 1] / d[i - 1];
            u[i] = prodb / d[i];
        }
    
        delta = d;
        delta[0] = a[0];
        for (i = 1; i <= n; i++)
            delta[i] = a[i] - pow(b[i - 1], 2) / delta[i - 1];

        v[n] = 1 / (ans[n] * delta[n]);
        v[0] = 1.0;
        prodb = 1.0;
        for (i = 1; i < n; i++) {
            prodb *= -b[n - i] / delta[n - i];
            v[n - i] = prodb * v[n];
        }
    }
    
    UNPROTECT(1);
}

/* R interface */

SEXP R_symtrisolve (SEXP a, SEXP b, SEXP tol)
{

    SEXP ans;
    R_xlen_t N; 

    N = XLENGTH(a);
    if (XLENGTH(b) != N - 1)
        error("incorrect length of b");
        
    if (!isReal(a))
        error("incorrect type of a");
        
    if (!isReal(b))
        error("incorrect type of b");

    if (!isReal(tol))
        error("incorrect type of tol");
        
    PROTECT(ans = allocMatrix(REALSXP, N, 2));
    C_symtrisolve(REAL(a), REAL(b), N - 1, REAL(tol)[0], REAL(ans));
    UNPROTECT(1);

    return(ans);
}

/*
    Compute t(X) %*% solve(T) * X
    for a symmetric N x N tridiagonal matrix T with 
    diagonals `a' and off-diagonals `b' and N x P matrix `X'
*/

SEXP R_symtrisolve_quadform (SEXP a, SEXP b, SEXP X, SEXP tol)
{

    SEXP ans, vu, cumsumvux;
    double *dans, *dx, dxA, dxA1, dxA2, *dvu, *dcs;
    R_xlen_t N, i, j;
    int p, pp, P;
    
    N = XLENGTH(a);

    if (NROW(X) != N)
        error("incorrect number of rows in X");
    if (!isReal(X))
        error("incorrect type of X");
    dx = REAL(X);
    P = (int) NCOL(X);

    PROTECT(ans = allocMatrix(REALSXP, P, P));
    dans = REAL(ans);

    if (XLENGTH(b) != N - 1)
        error("incorrect length of b");

    if (!isReal(a))
        error("incorrect type of a");

    if (!isReal(b))
        error("incorrect type of b");

    if (!isReal(tol))
        error("incorrect type of tol");

    PROTECT(vu = allocMatrix(REALSXP, N, 2));
    dvu = REAL(vu);
    C_symtrisolve(REAL(a), REAL(b), N - 1, REAL(tol)[0], dvu);
    PROTECT(cumsumvux = allocMatrix(REALSXP, N, 2));
    dcs = REAL(cumsumvux);

    for (p = 0; p < P * P; p++)
        dans[p] = 0.0;

    for (p = 0; p < P; p++) {
        i = 0;
        dcs[i] = dx[p * N + i] * dvu[i];
        dcs[N + i] = dx[p * N + i] * dvu[N + i];
        for (i = 1; i < N; i++) {
            dcs[i] = dcs[i - 1] + dx[p * N + i] * dvu[i];
            dcs[N + i] = dcs[N + i - 1] + dx[p * N + i] * dvu[N + i];
        }
        for (j = 0; j < N; j++) {
            dxA = 0.0;
            /*
            for (i = 0; i < N; i++) {
                if (j > i) {
                    dxA += dx[p * N + i] * dvu[N + i] * dvu[j];
                } else {
                    dxA += dx[p * N + i] * dvu[i] * dvu[N + j];
                }
            }
            */
            dxA1 = dcs[N + j];
            dxA2 = dcs[N - 1] - dcs[j];
            dxA = dxA1 * dvu[j] + dxA2 * dvu[N + j];
            for (pp = p; pp < P; pp++)
                dans[p * P + pp] += dxA * dx[pp * N + j];
        }
    }
    
    /* dans = t(X) %*% solve(T) * X is symmetric */
    for (p = 0; p < P; p++) {
        for (pp = p + 1; pp < P; pp++)
            dans[pp * P + p] = dans[p * P + pp];
    }
    UNPROTECT(3);
    return(ans);
}

SEXP R_cumsumrev (SEXP x) {

    double sum = 0.;
    SEXP ret;
    double *rx = REAL(x), *rs;

    PROTECT(ret = allocVector(REALSXP, XLENGTH(x)));
    rs = REAL(ret);

    for (R_xlen_t i = XLENGTH(x) - 1; i >= 0 ; i--) {
        sum += rx[i]; /* NA and NaN propagated */
        rs[i] = sum;
    }
    UNPROTECT(1);
    return(ret);
}
