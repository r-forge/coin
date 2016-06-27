
#include "libcoin_internal.h"

SEXP R_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset, 
                                      SEXP block, SEXP varonly);
SEXP R_PermutedLinearStatistic(SEXP LEV, SEXP x, SEXP y, SEXP weights, 
                               SEXP subset, SEXP block, SEXP B);
SEXP R_ExpectationCovarianceStatistic_2d(SEXP x, SEXP ix, SEXP y, SEXP iy,
                                         SEXP weights, SEXP subset, SEXP block,
                                         SEXP varonly);
SEXP R_PermutedLinearStatistic_2d(SEXP LEV, SEXP x, SEXP ix, SEXP y, SEXP iy,
                                  SEXP block, SEXP B);
                                  