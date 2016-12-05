
#include <R_ext/Rdynload.h>
#include <libcoin.h>

extern SEXP libcoin_R_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset,
                                                     SEXP block, SEXP varonly) {

    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP))
            R_GetCCallable("libcoin", "R_ExpectationCovarianceStatistic");
    return fun(x, y, weights, subset, block, varonly);
}

extern SEXP libcoin_R_PermutedLinearStatistic(SEXP LEV, SEXP x, SEXP y, SEXP weights,
                               SEXP subset, SEXP block, SEXP B, SEXP standardise, SEXP tol) {

    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP))
            R_GetCCallable("libcoin", "R_PermutedLinearStatistic");
    return fun(LEV, x, y, weights, subset, block, B, standardise, tol);
}

extern SEXP libcoin_R_ExpectationCovarianceStatistic_2d(SEXP x, SEXP ix, SEXP y, SEXP iy,
                                                 SEXP weights, SEXP subset, SEXP block,
                                                 SEXP varonly) {

    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP))
            R_GetCCallable("libcoin", "R_ExpectationCovarianceStatistic_2d");
    return fun(x, ix, y, iy, weights, subset, block, varonly);
}

extern SEXP libcoin_R_PermutedLinearStatistic_2d(SEXP LEV, SEXP x, SEXP ix, SEXP y, SEXP iy,
                                          SEXP block, SEXP B, SEXP standardise, SEXP tol) {

    static SEXP(*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = NULL;
    if(fun == NULL)
        fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP SEXP))
            R_GetCCallable("libcoin", "R_PermutedLinearStatistic_2d");
    return fun(LEV, x, ix, y, iy, block, B, standardise, tol);
}
