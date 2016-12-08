
#include <R_ext/Rdynload.h>  /* required by R */
#include <libcoinAPI.h>

SEXP R_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset,
                                      SEXP block, SEXP varonly, SEXP tol)
{
    return(libcoin_R_ExpectationCovarianceStatistic(x, y, weights, subset,
                                                    block, varonly, tol));
}
