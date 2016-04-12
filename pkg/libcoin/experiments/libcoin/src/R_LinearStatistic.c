
#include <R.h>
#include <Rinternals.h>
#include "LinearStatistic.h"
#include "helpers.h"
#include "Tables.h"
#include "Sums.h"

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

SEXP R_ExpectationInfluence(SEXP y) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_ExpectationInfluence(REAL(y), NROW(y), NCOL(y), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationInfluence_weights(SEXP y, SEXP weights) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_ExpectationInfluence_weights(REAL(y), NROW(y), NCOL(y), 
                                   INTEGER(weights), 
                                   C_sum(INTEGER(weights), LENGTH(weights)), 
                                   REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationInfluence_subset(SEXP y, SEXP subset) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_ExpectationInfluence_subset(REAL(y), NROW(y), NCOL(y), 
                                  INTEGER(subset), LENGTH(subset), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationInfluence_weights_subset(SEXP y, SEXP weights, SEXP subset) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_ExpectationInfluence_weights_subset(REAL(y), NROW(y), NCOL(y), 
                                          INTEGER(weights), 
                                          C_sum_subset(INTEGER(weights), 
                                                       LENGTH(weights), 
                                                       INTEGER(subset), 
                                                       LENGTH(subset)),
                                          INTEGER(subset), LENGTH(subset), 
                                          REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceInfluence(SEXP y) {

    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence(y));
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(y)));
    C_CovarianceInfluence(REAL(y), NROW(y), NCOL(y), REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_CovarianceInfluence_weights(SEXP y, SEXP weights) {

    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_weights(y, weights));
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(y)));
    C_CovarianceInfluence_weights(REAL(y), NROW(y), NCOL(y), 
                                   INTEGER(weights), 
                                   C_sum(INTEGER(weights), LENGTH(weights)), 
                                   REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_CovarianceInfluence_subset(SEXP y, SEXP subset) {

    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_subset(y, subset));
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(y)));
    C_CovarianceInfluence_subset(REAL(y), NROW(y), NCOL(y), 
                                 INTEGER(subset), LENGTH(subset), 
                                 REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_CovarianceInfluence_weights_subset(SEXP y, SEXP weights, SEXP subset) {

    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_weights_subset(y, weights, subset));
    PROTECT(ans = allocMatrix(REALSXP, NCOL(y), NCOL(y)));
    C_CovarianceInfluence_weights_subset(REAL(y), NROW(y), NCOL(y), 
                                          INTEGER(weights), 
                                          C_sum_subset(INTEGER(weights), 
                                                       LENGTH(weights), 
                                                       INTEGER(subset), 
                                                       LENGTH(subset)),
                                          INTEGER(subset), LENGTH(subset), 
                                          REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_ExpectationX(SEXP x) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_ExpectationX(REAL(x), NROW(x), NCOL(x), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationX_weights(SEXP x, SEXP weights) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_ExpectationX_weights(REAL(x), NROW(x), NCOL(x), 
                                   INTEGER(weights), 
                                   C_sum(INTEGER(weights), LENGTH(weights)), 
                                   REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationX_subset(SEXP x, SEXP subset) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_ExpectationX_subset(REAL(x), NROW(x), NCOL(x), 
                                  INTEGER(subset), LENGTH(subset), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationX_weights_subset(SEXP x, SEXP weights, SEXP subset) {

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_ExpectationX_weights_subset(REAL(x), NROW(x), NCOL(x), 
                                          INTEGER(weights), 
                                          C_sum_subset(INTEGER(weights), 
                                                       LENGTH(weights), 
                                                       INTEGER(subset), 
                                                       LENGTH(subset)),
                                          INTEGER(subset), LENGTH(subset), 
                                          REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceX(SEXP x) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(x)));
    C_CovarianceX(REAL(x), NROW(x), NCOL(x), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceX_weights(SEXP x, SEXP weights) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(x)));
    C_CovarianceX_weights(REAL(x), NROW(x), NCOL(x), 
                                   INTEGER(weights), 
                                   REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceX_subset(SEXP x, SEXP subset) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(x)));
    C_CovarianceX_subset(REAL(x), NROW(x), NCOL(x), 
                                 INTEGER(subset), LENGTH(subset), 
                                 REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceX_weights_subset(SEXP x, SEXP weights, SEXP subset) {

    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(x)));
    C_CovarianceX_weights_subset(REAL(x), NROW(x), NCOL(x), 
                                          INTEGER(weights), 
                                          INTEGER(subset), LENGTH(subset), 
                                          REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationLinearStatistic(SEXP x, SEXP y) {

    SEXP ans, ExpInf, ExpX;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    PROTECT(ExpInf = R_ExpectationInfluence(y));
    PROTECT(ExpX = R_ExpectationX(x));    
    C_ExpectationLinearStatistic(NCOL(x), NCOL(y), REAL(ExpInf), REAL(ExpX), REAL(ans));
    UNPROTECT(3);
    return(ans);
}

SEXP R_ExpectationLinearStatistic_weights(SEXP x, SEXP y, SEXP weights) {

    SEXP ans, ExpInf, ExpX;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    PROTECT(ExpInf = R_ExpectationInfluence_weights(y, weights));
    PROTECT(ExpX = R_ExpectationX_weights(x, weights));    
    C_ExpectationLinearStatistic(NCOL(x), NCOL(y), REAL(ExpInf), REAL(ExpX), REAL(ans));
    UNPROTECT(3);
    return(ans);
}

SEXP R_ExpectationLinearStatistic_subset(SEXP x, SEXP y, SEXP subset) {

    SEXP ans, ExpInf, ExpX;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    PROTECT(ExpInf = R_ExpectationInfluence_subset(y, subset));
    PROTECT(ExpX = R_ExpectationX_subset(x, subset));    
    C_ExpectationLinearStatistic(NCOL(x), NCOL(y), REAL(ExpInf), REAL(ExpX), REAL(ans));
    UNPROTECT(3);
    return(ans);
}

SEXP R_ExpectationLinearStatistic_weights_subset(SEXP x, SEXP y, SEXP weights, SEXP subset) {

    SEXP ans, ExpInf, ExpX;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    PROTECT(ExpInf = R_ExpectationInfluence_weights_subset(y, weights, subset));
    PROTECT(ExpX = R_ExpectationX_weights_subset(x, weights, subset));    
    C_ExpectationLinearStatistic(NCOL(x), NCOL(y), REAL(ExpInf), REAL(ExpX), REAL(ans));
    UNPROTECT(3);
    return(ans);
}

SEXP R_CovarianceLinearStatistic(SEXP CovInf, SEXP ExpX, SEXP CovX, SEXP sumweights) {

    SEXP ans, PQPQ_tmp, PQQ_tmp;
    int P, Q, PQ;
    
    P = LENGTH(ExpX);
    Q = NCOL(CovInf);
    PQ = P * Q;
    PROTECT(ans = allocMatrix(REALSXP, PQ, PQ));
    PROTECT(PQPQ_tmp = allocMatrix(REALSXP, PQ, PQ));
    PROTECT(PQQ_tmp = allocMatrix(REALSXP, PQ, Q));
    C_CovarianceLinearStatistic(P, Q, REAL(CovInf), REAL(ExpX), REAL(CovX), 
                                INTEGER(sumweights)[0], REAL(PQPQ_tmp), REAL(PQQ_tmp), REAL(ans));
    UNPROTECT(3);
    return(ans);
}
