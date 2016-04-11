
#include <R.h>
#include <Rinternals.h>
#include "LinearStatistic.h"
#include "helpers.h"

SEXP R_LinearStatistic(SEXP x, SEXP y) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(x)));
    C_LinearStatistic(REAL(x), NROW(x), NCOL(x), REAL(y), NCOL(y), REAL(ans));
    UNPROTECT(1);
    return(ans);
}
