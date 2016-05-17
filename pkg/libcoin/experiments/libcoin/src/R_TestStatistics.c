
#include "libcoin_internal.h"
#include "TestStatistics.h"
#include "R_LinearStatistic.h"
#include "Utils.h"
#include "Contrasts.h"

SEXP R_quadform(SEXP LinearStatistic, SEXP Expectation, SEXP MPinv)
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = CR_quadform(LinearStatistic, Expectation, MPinv);
    UNPROTECT(1);
    return(ans);
}

SEXP R_maxtype(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP type, SEXP tol)
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    if (INTEGER(type)[0] == ALTERNATIVE_twosided)
        REAL(ans)[0] = CR_maxabsstat(LinearStatistic, Expectation, CoVariance, tol);
    if (INTEGER(type)[0] == ALTERNATIVE_less)
        REAL(ans)[0] = CR_maxstat(LinearStatistic, Expectation, CoVariance, tol);
    if (INTEGER(type)[0] == ALTERNATIVE_greater)
        REAL(ans)[0] = CR_minstat(LinearStatistic, Expectation, CoVariance, tol);
    UNPROTECT(1);
    return(ans);
}

SEXP R_ordered(SEXP x, SEXP y, SEXP weights,
               SEXP subset, SEXP block, SEXP varonly, SEXP teststat, SEXP tol)
{
    SEXP LEV, ans, wmax, maxstat;
    double *contrasts;
    int P, Q;

    PROTECT(LEV = R_ExpectationCovarianceStatistic(x, y, weights, subset,  block, varonly));
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, wmax = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 1, maxstat = allocVector(REALSXP, 1));
    
    P = NLEVELS(x);
    Q = (int) LENGTH(VECTOR_ELT(LEV, 0)) / P;

    contrasts = Calloc(P * (P - 1), double);
    
    for (int p = 0; p < P * (P - 1); p++) contrasts[p] = 0.0;
    for (int p = 0; p < (P - 1); p++) {
        for (int pp = 0; pp <= p; pp++)
            contrasts[pp + p * P] = 1.0;
    }

    if (INTEGER(teststat)[0] == 1) {                            
        C_contrasts_marginal_maxabsstat(REAL(VECTOR_ELT(LEV, 0)),
                                        REAL(VECTOR_ELT(LEV, 1)),
                                        REAL(VECTOR_ELT(LEV, 2)),
                                        contrasts, P, Q, 
                                        P - 1,
                                        REAL(tol)[0], 
                                        INTEGER(wmax), REAL(maxstat));
    } else {
        C_contrasts_marginal_quadform(REAL(VECTOR_ELT(LEV, 0)),
                                      REAL(VECTOR_ELT(LEV, 1)),
                                      REAL(VECTOR_ELT(LEV, 2)),
                                      contrasts, P, Q, 
                                      P - 1,
                                      REAL(tol)[0], 
                                      INTEGER(wmax), REAL(maxstat));
    }
    Free(contrasts);
    UNPROTECT(2);
    return(ans);
}                                      
