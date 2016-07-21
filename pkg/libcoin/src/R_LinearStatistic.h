
#include "libcoin_internal.h"

extern SEXP R_ExpectationCovarianceStatistic
(
    const SEXP x, 
    const SEXP y, 
    const SEXP weights, 
    const SEXP subset, 
    const SEXP block, 
    const SEXP varonly
);

extern SEXP R_PermutedLinearStatistic
(
    const SEXP LEV, 
    const SEXP x, 
    const SEXP y, 
    const SEXP weights, 
    const SEXP subset, 
    const SEXP block, 
    const SEXP B
);

extern SEXP R_ExpectationCovarianceStatistic_2d
(
    const SEXP x,
    const SEXP ix, 
    const SEXP y, 
    const SEXP iy,
    const SEXP weights, 
    const SEXP subset, 
    const SEXP block,
    const SEXP varonly
);

extern SEXP R_PermutedLinearStatistic_2d
(
    const SEXP LEV, 
    const SEXP x, 
    const SEXP ix,
    const SEXP y, 
    const SEXP iy,
    const SEXP block, 
    const SEXP B
);
                                  