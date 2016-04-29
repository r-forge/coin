
#include "libcoin.h"
#include "helpers.h"

/* MP inv of symmetric matrix in lower triangular packed form */

SEXP R_MPinv_sym (SEXP x, SEXP tol) {

    SEXP ans, rank, MP;
    double *val, *vec, *dMP, dtol, *rx, *work, valinv;
    int n, valzero = 0, info = 0, kn;

    n = (int) (.5 + sqrt(.25 + 2 * LENGTH(x))) - 1;
    
    rx = Calloc(LENGTH(x), double);
    Memcpy(rx, REAL(x), LENGTH(x));
    work = Calloc(3 * n, double);
    val = Calloc(n, double);
    vec = Calloc(n * n, double);
    
    F77_CALL(dspev)("V", "L", &n, rx, val, vec, &n, work,
                    &info);
                                            
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, MP = allocVector(REALSXP, n * (n + 1) / 2));
    SET_VECTOR_ELT(ans, 1, rank = allocVector(INTSXP, 1));
    dMP = REAL(MP);
    
    dtol = val[n - 1] * REAL(tol)[0];

    for (int k = 0; k < n; k++)
        valzero += (val[k] < dtol); 
    INTEGER(rank)[0] = n - valzero;

    for (int i = 0; i < n * (n + 1) / 2; i++) dMP[i] = 0.0;
    
    for (int k = valzero; k < n; k++) {
        valinv = 1 / val[k];
        kn = k * n;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                /* MP is symmetric */
                dMP[S(i, j, n)] += valinv * vec[kn + i] * vec[kn + j];
            }
        }
    }
    Free(rx); Free(work); Free(val); Free(vec);
    UNPROTECT(1);
    return(ans);
}
