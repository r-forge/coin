
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <R_ext/Rdynload.h>  /* required by R */

#include <libcoinAPI.h>

SEXP myR_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset,
                                        SEXP block, SEXP varonly) {

   return(libcoin_R_ExpectationCovarianceStatistic(x, y, weights, subset,
                                           block, varonly));
}
