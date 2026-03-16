/**
    Some additional functionality for package 'coin'

    *\file Helpers.c
*/

#include "coin_common.h"
#include <R_ext/Rdynload.h>
#include <libcoinAPI.h>

SEXP R_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset,
                                      SEXP block, SEXP varonly, SEXP tol) {
    return(libcoin_R_ExpectationCovarianceStatistic(x, y, weights, subset,
                                                    block, varonly, tol));
}

SEXP R_PermutedLinearStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset,
                               SEXP block, SEXP nresample) {
    return(libcoin_R_PermutedLinearStatistic(x, y, weights, subset,
                                             block, nresample));
}

SEXP R_quadform(SEXP linstat, SEXP expect, SEXP MPinv_sym) {
    return(libcoin_R_quadform(linstat, expect, MPinv_sym));
}

SEXP R_kronecker(SEXP A, SEXP B) {
    return(libcoin_R_kronecker(A, B));
}

SEXP R_MPinv_sym(SEXP x, SEXP n, SEXP tol) {
    return(libcoin_R_MPinv_sym(x, n, tol));
}

SEXP R_unpack_sym(SEXP x, SEXP names, SEXP diagonly) {
    return(libcoin_R_unpack_sym(x, names, diagonly));
}

SEXP R_maxstattrafo(SEXP x, SEXP cutpoints) {

    int n, nc, jn;
    SEXP ans;
    double *dans, *dx, *dcutpoints, cj;

    if (!isReal(x) || !isReal(cutpoints))
        error("x or cutpoints are not of type REALSXP");

    n = LENGTH(x);
    nc = LENGTH(cutpoints);
    PROTECT(ans = allocMatrix(REALSXP, n, nc));
    dans = REAL(ans);
    dx = REAL(x);
    dcutpoints = REAL(cutpoints);

    for (int j = 0; j < nc; j++) {
        jn = j * n;
        cj = dcutpoints[j];
        for (int i = 0; i < n; i++) {
            if (ISNAN(dx[i])) {
                dans[jn + i] = dx[i];
            } else if (dx[i] > cj) {
                dans[jn + i] = 0.0;
            } else {
                dans[jn + i] = 1.0;
            }
        }
    }
    UNPROTECT(1);
    return(ans);
}
