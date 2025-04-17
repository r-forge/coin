
#include "lehmann.h"

int NROW
(
    SEXP x
) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(XLENGTH(x));
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[0]);
    return(INTEGER(a)[0]);
}

/* NCOL */

int NCOL
(
    SEXP x
) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[1]);
    return(INTEGER(a)[1]);
}


/* 
   Compute inverse of a symmetric n x n tridiagonal matrix A with 
   diagonals `a' and off-diagonals `b'
   
   A^{-1}_ij = u_i v_j for i < j
   
   as described in DOI:10.1137/0613045
   
   The function returns a real n x 2 matrix with 1st column u
   and second column v
   
*/

void C_symtrisolve (double *a, double *b, int n, double *ans)
{

    SEXP Rd;
    double *d, prodb;
    int i;
    
    PROTECT(Rd = allocVector(REALSXP, n));
    d = REAL(Rd);
    
    d[n] = a[n];
    for (i = n - 1; i >= 0; i--)
        d[i] = a[i] - b[i] * b[i] / d[i + 1];
        
    ans[0] = 1 / d[0];
    prodb = 1.0;
    for (i = 1; i <= n; i++) {
        prodb *= -b[i - 1] / d[i - 1];
        ans[i] = prodb;
        ans[i] /= d[i];
    }
    
    d[0] = a[0];
    for (i = 1; i <= n; i++)
        d[i] = a[i] - b[i - 1] * b[i - 1] / d[i - 1];

    ans[2 * n + 1] = 1 / (ans[n] * d[n]);
    ans[n + 1] = 1.0;
    prodb = 1.0;
    for (i = 1; i < n; i++) {
        prodb *= -b[n - i] / d[n - i];
        ans[n + (n - i) + 1] = prodb;
        ans[n + (n - i) + 1] /= (ans[n] * d[n]);
    }
    
    UNPROTECT(1);
}

/* R interface */

SEXP R_symtrisolve (SEXP a, SEXP b)
{

    SEXP ans;
    int N;

    N = LENGTH(a);
    if (LENGTH(b) != N - 1)
        error("incorrect length of b");
        
    if (!isReal(a))
        error("incorrect type of a");
        
    if (!isReal(b))
        error("incorrect type of b");
        
    PROTECT(ans = allocMatrix(REALSXP, N, 2));
    C_symtrisolve(REAL(a), REAL(b), N - 1, REAL(ans));
    UNPROTECT(1);

    return(ans);
}

/*
    Compute t(X) %*% solve(A) * X
    for a symmetric n x n tridiagonal matrix A with 
    diagonals `a' and off-diagonals `b'
*/

SEXP R_symtrisolve_quadform (SEXP a, SEXP b, SEXP X)
{

    SEXP ans, vu, xA;
    double *dans, *dx, *dxA, *dvu;
    int N, i, j, p, pp, P;
    
    N = LENGTH(a);
    
    if (NROW(X) != N)
        error("incorrect number of rows in X");
    if (!isReal(X))
        error("incorrect type of X");
    dx = REAL(X);
    P = NCOL(X);

    PROTECT(ans = allocMatrix(REALSXP, P, P));
    dans = REAL(ans);

    if (LENGTH(b) != N - 1)
        error("incorrect length of b");

    if (!isReal(a))
        error("incorrect type of a");

    if (!isReal(b))
        error("incorrect type of b");

    PROTECT(vu = allocMatrix(REALSXP, N, 2));
    dvu = REAL(vu);
    C_symtrisolve(REAL(a), REAL(b), N - 1, dvu);

    for (p = 0; p < P * P; p++)
        dans[p] = 0.0;

    PROTECT(xA = allocVector(REALSXP, N - 1));
    dxA = REAL(xA);
        
    for (p = 0; p < P; p++) {
        for (j = 0; j < N; j++) {
            dxA[j] = 0.0;
            for (i = 0; i < N; i++) {
                if (j > i) {
                    dxA[j] += dx[p * N + i] * dvu[N + i] * dvu[j];
                } else {
                    dxA[j] += dx[p * N + i] * dvu[i] * dvu[N + j];
                }
            }
            for (pp = p; pp < P; pp++)
                dans[p * P + pp] += dxA[j] * dx[pp * N + j];
        }
    }
    
    /* dans is symmetric */
    for (p = 0; p < P; p++) {
        for (pp = p + 1; pp < P; pp++)
            dans[pp * P + p] = dans[p * P + pp];
    }
    UNPROTECT(3);
    return(ans);
}
