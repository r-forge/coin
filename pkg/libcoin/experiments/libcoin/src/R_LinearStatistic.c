
#include <R.h>
#include <Rinternals.h>
#include "LinearStatistic.h"
#include "helpers.h"
#include "Tables.h"

SEXP R_LinearStatistic(SEXP x, SEXP y) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(x)));
    C_LinearStatistic(REAL(x), NROW(x), NCOL(x), REAL(y), NCOL(y), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_LinearStatistic_weights(SEXP x, SEXP y, SEXP weights) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(x)));
    C_LinearStatistic_weights(REAL(x), NROW(x), NCOL(x), REAL(y), NCOL(y), 
                              INTEGER(weights), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_LinearStatistic_subset(SEXP x, SEXP y, SEXP subset) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(x)));
    C_LinearStatistic_subset(REAL(x), NROW(x), NCOL(x), REAL(y), NCOL(y), 
                             INTEGER(subset), LENGTH(subset), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_LinearStatistic_weights_subset(SEXP x, SEXP y, SEXP weights, SEXP subset) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(x)));
    C_LinearStatistic_weights_subset(REAL(x), NROW(x), NCOL(x), 
                                     REAL(y), NCOL(y), INTEGER(weights), 
                                     INTEGER(subset), LENGTH(subset), 
                                     REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_LinearStatistic_2d(SEXP x, SEXP y, SEXP fx, SEXP fy) {

    SEXP table, ans;
    
    PROTECT(table = R_2dtable(fy, fx));
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(x)));
    C_LinearStatistic_2d(REAL(x), NROW(x), NCOL(x), REAL(y), NROW(y), NCOL(y), 
                         INTEGER(table), REAL(ans));
    UNPROTECT(2);
    return(ans);
}
