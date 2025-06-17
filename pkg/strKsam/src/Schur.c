
#define STRICT_R_HEADERS
#include <R.h>
#include <Rinternals.h>
/* C_symtrisolve */

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

    /* d vec */
    
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
        /* u vec */
        
        u[0] = 1 / d[0];
        prodb = 1.0;
        for (i = 1; i <= n; i++) {
            prodb *= -b[i - 1] / d[i - 1];
            u[i] = prodb / d[i];
        }
        
        delta = d;
        /* delta vec */
        
        delta[0] = a[0];
        for (i = 1; i <= n; i++)
            delta[i] = a[i] - pow(b[i - 1], 2) / delta[i - 1];
        
        /* v vec */
        
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

/* NROW */

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

/* wcrossprod */

SEXP R_wcrossprod (SEXP a, SEXP b, SEXP X, SEXP tol)
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

    /* lower wcrossprod */
    
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
            dxA1 = dcs[N + j];
            dxA2 = dcs[N - 1] - dcs[j];
            dxA = dxA1 * dvu[j] + dxA2 * dvu[N + j];
            for (pp = p; pp < P; pp++)
                dans[p * P + pp] += dxA * dx[pp * N + j];
        }
    }
    
    /* upper wcrossprod */
    
    for (p = 0; p < P; p++) {
        for (pp = p + 1; pp < P; pp++)
            dans[pp * P + p] = dans[p * P + pp];
    }
    

    UNPROTECT(3);
    return(ans);
}

